#ifndef PARALLEL_TEMPERING_C
#define PARALLEL_TEMPERING_C

#include "../include/macro.h"
#include "../include/gauge_conf.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/memalign.h"
#include "../include/random.h"


// one group of n_pairs swaps: pairs (a, a+1) with a = first, first + 2, ...
//
// Swapping replicas a and b = a+1 exchanges the replica indices of A = GC[replica_to_gc[a]] and B = GC[replica_to_gc[b]]
// The action difference
//   swapped - unswapped = beta * (D_A - D_B),
//   D_X = sum_{r in swap_rect, j != mu} [K_{c_a}(r,mu,j) - K_{c_b}(r,mu,j)] ReTr(Z plaq_X(r,mu,j))/N_c,
// is computed as two independent sums, each one on a single configuration.
static void swap_group(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                       Rectangle const *const swap_rectangle, Acc_Utils *acc_counters,
                       int const first, int const n_pairs, double *const D)
	{
	int const Nr = param->d_N_replica_pt;
	int const mu = param->d_defect_dir; // swaps change only the plaquettes (mu, j), j != mu, around the defect
	long const vol_rect = swap_rectangle->d_vol_rect;
	int const *const replica_to_gc = acc_counters->replica_to_gc;
	double *const aux_p = acc_counters->metro_swap_prob; // alias for the Metropolis probabilities of the pairs

	for(int a = 0; a < Nr; a++)
		D[a] = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:D[:Nr])
	#endif
	for(long s = 0; s < 2 * n_pairs * vol_rect; s++)
		{
		// s = (2 * k + side) * vol_rect + n
		long const n = s % vol_rect;
		long const t = (s - n) / vol_rect;
		int const k = (int) (t / 2);                   // index of the pair in the group
		int const a = 2 * k + first;                   // index of the first replica in the pair
		int const b = a + 1;                           // index of the second replica in the pair
		int const x = (t % 2 == 0) ? a : b;            // replica index of the configuration summed
		long const r = swap_rectangle->rect_sites[n];  // site index of the configuration summed

		double const c_a = param->d_pt_bound_cond_coeff[a];
		double const c_b = param->d_pt_bound_cond_coeff[b];
		Gauge_Conf const *const GC_x = &(GC[replica_to_gc[x]]);

		for(int i = 0; i < STDIM - 1; i++)
			{
			int const j = geo->d_orth_dir[mu][i];
			D[x] += swap_action_site(GC_x, geo, param, r, mu, j, c_a, c_b);
			}
		}

	// Metropolis tests and swaps (serial with Philox: the Metropolis tests use the single host generator)
	#if defined(OPENMP_MODE) && !defined(RNG_PHILOX_MODE)
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int k = 0; k < n_pairs; k++)
		{
		int const a = 2 * k + first;
		int const b = a + 1;

		aux_p[a] = param->d_beta * (D[a] - D[b]); // swapped action - unswapped action
		#ifdef MULTICANONICAL_MODE
		aux_p[a] += delta_topo_potential_swap(&(GC[replica_to_gc[a]]), &(GC[replica_to_gc[b]]), param);
		#endif
		aux_p[a] = exp(-aux_p[a]);                                // Metropolis swap probability
		metropolis_single_swap(GC, a, b, aux_p[a], acc_counters); // Metropolis step
		}
	}


// swaps of all neighbouring replica pairs, in two groups (even and odd pairs) in random order
void swap(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
          Rectangle const *const swap_rectangle, Acc_Utils *acc_counters)
	{
	// N_replica_pt - 1 is the total number of swaps
	int const num_swaps = param->d_N_replica_pt - 1;

	int const is_even = num_swaps % 2;              // check if num_swaps is even or not
	int const num_even = (num_swaps + is_even) / 2; // number of swaps for even replica
	int const num_odd = (num_swaps - is_even) / 2;  // number of swaps for odd replica

	// to be sure detailed balance is satisfied, choose randomly whether to swap first odd or even copies
	int is_even_first, num_swaps_1, num_swaps_2;
	if(casuale(rng_host()) < 0.5) // first swap all even copies, then all odd copies
		{
		is_even_first = 0;
		num_swaps_1 = num_even;
		num_swaps_2 = num_odd;
		}
	else // first swap all odd copies, then all even copies
		{
		is_even_first = 1;
		num_swaps_1 = num_odd;
		num_swaps_2 = num_even;
		}

	// local sums of the action differences, one per replica index
	double *D;
	allocate_array_double(&D, param->d_N_replica_pt, __FILE__, __LINE__);

	swap_group(GC, geo, param, swap_rectangle, acc_counters, is_even_first, num_swaps_1, D);     // first group
	swap_group(GC, geo, param, swap_rectangle, acc_counters, 1 - is_even_first, num_swaps_2, D); // second group: reverse parity

	free(D);
	}


// contribution of the site r, plane (i,j) of the configuration GC to the swap action difference of replicas with
// defect coefficients c_a and c_b: [K_{c_a}(r,i,j) - K_{c_b}(r,i,j)] ReTr(Z plaq(r,i,j))/N_c
ACC_ROUTINE_SEQ
double swap_action_site(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param,
                        long const r, int const i, int const j, double const c_a, double const c_b)
	{
	long const rpi = nnp(geo, r, i);
	long const rpj = nnp(geo, r, j);
	double const delta_K = plaquette_defect_coeff(geo, param, r, rpi, rpj, i, j, c_a) - plaquette_defect_coeff(geo, param, r, rpi, rpj, i, j, c_b);

	return delta_K * plaquettep(GC, geo, param, r, i, j); // plaquettep includes the twist factors
	}

// Metropolis step to swap replicas a and b with probability p: the configurations with replica index a and b
// exchange their replica indices (configurations, twist factors, charges and translations stay in their position in GC)
void metropolis_single_swap(Gauge_Conf *const GC, int const a, int const b, double const p, Acc_Utils *acc_counters)
	{
	// Metropolis test: p<1 => acc=1 with probability p, p>=1 => acc=1 (already assigned)
	int acc = 1;
	if(p < 1 && casuale(rng_host()) > p) acc = 0;

	if(acc == 1)
		{
		int const gc_a = acc_counters->replica_to_gc[a];
		int const gc_b = acc_counters->replica_to_gc[b];

		GC[gc_a].replica_index = b;
		GC[gc_b].replica_index = a;
		acc_counters->replica_to_gc[a] = gc_b;
		acc_counters->replica_to_gc[b] = gc_a;

		// increase counter of successful swaps for replicas (a, a+1)
		acc_counters->num_accepted_swap[a]++;
		}

	// increase counter of tried swaps for replicas (a, a+1) and update acceptance ratio
	acc_counters->num_swap[a]++;
	acc_counters->acc_swap[a] = (double) acc_counters->num_accepted_swap[a] / (double) acc_counters->num_swap[a];
	}


// translation of one lattice spacing of the configuration, including the Z factors
// direction is chosen randomly, verse is always positive
void conf_translation(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	// extract random direction
	int const dir = (int) (STDIM * casuale(rng_host()));
	GC->translation[dir] += 1;

	// translation in direction +dir, including the Z factors and cold lattice
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < (param->d_n_planes + 1) * (param->d_volume); s++)
		{
		// s = j * volume + r
		long const r = s % (param->d_volume);
		int const j = (int) ((s - r) / (param->d_volume));
		if(j < STDIM)
			{
			equal(&(GC->lattice[r][j]), &(GC->lattice_copy[nnm(geo, r, dir)][j]));
			#ifdef MULTICANONICAL_MODE
			equal(&(GC->lattice_cold[r][j]), &(GC->lattice_copy_cold[nnm(geo, r, dir)][j]));
			#endif
			}
		GC->Z[r][j] = GC->Z_copy[nnm(geo, r, dir)][j];
		}

	// swap translated cold lattice with previous cold lattice
	#ifdef MULTICANONICAL_MODE
	GAUGE_GROUP **aux;
	aux = GC->lattice_cold;
	GC->lattice_cold = GC->lattice_copy_cold;
	GC->lattice_copy_cold = aux;
	#endif

	// update the auxiliary confs
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < (param->d_n_planes + 1) * (param->d_volume); s++)
		{
		// s = j * volume + r
		long const r = s % (param->d_volume);
		int const j = (int) ((s - r) / (param->d_volume));
		if(j < STDIM)
			{
			equal(&(GC->lattice_copy[r][j]), &(GC->lattice[r][j]));
			}
		GC->Z_copy[r][j] = GC->Z[r][j];
		}
	}


void init_swap_track_file(FILE **swaptrackfilep, GParam const *const param)
	{
	if(param->d_N_replica_pt == 1)
		{
		*swaptrackfilep = NULL;
		return;
		}
	int write_header = 0;
	if(param->d_start == 2) // starting run from saved conf
		{
		*swaptrackfilep = fopen(param->d_swap_tracking_file, "r");
		if(*swaptrackfilep != NULL) // file exists -> close it and re-open it in append mode
			{
			fclose(*swaptrackfilep);
			*swaptrackfilep = fopen(param->d_swap_tracking_file, "a");
			}
		else // file does not exist -> create it
			{
			*swaptrackfilep = fopen(param->d_swap_tracking_file, "w");
			write_header = 1;
			}
		}
	else // starting run from scratch
		{
		*swaptrackfilep = fopen(param->d_swap_tracking_file, "w");
		write_header = 1;
		}
	REQUIRE(*swaptrackfilep != NULL, "failed to open the swap tracking file %s", param->d_swap_tracking_file);
	if(write_header == 1)
		{
		fprintf(*swaptrackfilep, "# MC_step    conf_labels");
		#ifdef MULTICANONICAL_MODE
		fprintf(*swaptrackfilep, "    conf_charges");
		#endif
		fprintf(*swaptrackfilep, "\n");
		fflush(*swaptrackfilep);
		}
	}


void init_acc_utils(Acc_Utils *acc_counters, GParam const *const param)
	{
	// map replica index -> position in GC, initially the identity (GC[i].replica_index = i)
	allocate_array_int(&(acc_counters->replica_to_gc), param->d_N_replica_pt, __FILE__, __LINE__);
	for(int a = 0; a < param->d_N_replica_pt; a++)
		acc_counters->replica_to_gc[a] = a;

	if(param->d_N_replica_pt == 1)
		{
		acc_counters->num_accepted_swap = NULL;
		acc_counters->num_swap = NULL;
		acc_counters->acc_swap = NULL;
		acc_counters->metro_swap_prob = NULL;
		acc_counters->swaptrackfilep = NULL;
		}
	else
		{
		allocate_array_long(&(acc_counters->num_accepted_swap), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		allocate_array_long(&(acc_counters->num_swap), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		allocate_array_double(&(acc_counters->acc_swap), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		allocate_array_double(&(acc_counters->metro_swap_prob), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		for(int i = 0; i < (param->d_N_replica_pt - 1); i++)
			{
			acc_counters->num_accepted_swap[i] = 0;
			acc_counters->num_swap[i] = 0;
			acc_counters->acc_swap[i] = 0.0;
			}
		init_swap_track_file(&(acc_counters->swaptrackfilep), param);
		}
	#ifdef MULTICANONICAL_MODE
	init_multicanonic_acc_utils(acc_counters, param);
	#endif
	}


void free_acc_utils(Acc_Utils *acc_counters, GParam const *const param)
	{
	free(acc_counters->replica_to_gc);

	if(param->d_N_replica_pt > 1)
		{
		free(acc_counters->num_accepted_swap);
		free(acc_counters->num_swap);
		free(acc_counters->acc_swap);
		free(acc_counters->metro_swap_prob);
		fclose(acc_counters->swaptrackfilep);
		}
	else
		{
		(void) acc_counters; // to suppress compiler warning of unused variable
		}
	#ifdef MULTICANONICAL_MODE
	free_multicanonic_acc_utils(acc_counters);
	#endif
	}


void print_acceptances(Acc_Utils const *const acc_counters, GParam const *const param)
	{
	if(param->d_N_replica_pt == 1)
		{
		(void) acc_counters; // to suppress compiler warning of unused variable
		return;
		}
	FILE *fp = fopen(param->d_swap_acc_file, "w");
	REQUIRE(fp != NULL, "failed to open swap acceptances file %s", param->d_swap_acc_file);

	fprintf(fp, "# %-6s %-6s %-12s %-12s %-12s %-12s %-12s %-12s\n", "from", "to", "c_1", "c_2", "acc(%)", "err_acc(%)", "accepted", "tried");
	double acc, err_acc;
	for(int k = 0; k < param->d_N_replica_pt - 1; k++)
		{
		#ifdef DEBUG
		ASSERT(acc_counters->num_accepted_swap[k] <= acc_counters->num_swap[k], "invalid swap counters for %d-th replica pair", k);
		#endif
		if(acc_counters->num_swap[k] > 1)
			{
			acc = acc_counters->acc_swap[k];
			err_acc = sqrt(acc * (1.0 - acc) / ((double) acc_counters->num_swap[k] - 1.0));
			}
		else if(acc_counters->num_swap[k] == 1)
			{
			acc = (double) acc_counters->num_accepted_swap[k];
			err_acc = 1.0;
			}
		else
			{
			acc = 0.0;
			err_acc = 0.0;
			}
		fprintf(fp,
		        "  %-6d %-6d %12.6f %12.6f %12.3f %12.3f %12ld %12ld\n",
		        k,
		        k + 1,
		        param->d_pt_bound_cond_coeff[k],
		        param->d_pt_bound_cond_coeff[k + 1],
		        100.0 * acc,
		        100.0 * err_acc,
		        acc_counters->num_accepted_swap[k],
		        acc_counters->num_swap[k]);
		}
	fclose(fp);
	}


void print_conf_labels(Gauge_Conf const *const GC, GParam const *const param, Acc_Utils const *const acc_counters)
	{
	if(param->d_N_replica_pt == 1)
		{
		// to suppress compiler warning of unused variable
		(void) GC;
		(void) acc_counters;
		return;
		}
	// for every replica index, label of the configuration currently at that index
	// (= its position in GC because configurations never move)
	FILE *fp = acc_counters->swaptrackfilep;
	fprintf(fp, "%9ld ", GC[0].update_index);
	for(int a = 0; a < param->d_N_replica_pt; a++)
		{
		fprintf(fp, "%4d ", acc_counters->replica_to_gc[a]);
		#ifdef MULTICANONICAL_MODE
		fprintf(fp, "% 9.6f ", GC[acc_counters->replica_to_gc[a]].stored_topcharge);
		#endif
		}
	fprintf(fp, "\n");
	fflush(fp);
	}


#endif
