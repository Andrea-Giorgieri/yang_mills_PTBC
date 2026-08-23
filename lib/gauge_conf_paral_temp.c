#ifndef PARALLEL_TEMPERING_C
#define PARALLEL_TEMPERING_C

#include"../include/macro.h"

#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"../include/memalign.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/random.h"
#include"../include/function_pointers.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"


// swaps are parallelized, evaluation of swap probabilities is parallelized
void swap(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
          Rectangle const *const swap_rectangle, Acc_Utils *acc_counters)
	{
	// Just an alias to be used in reduction clause for OpenMP. icc gives error during
	// optimization if reduction(+:acc_counters->metro_swap_prob[:num_swaps]) is used
	double *aux_p = acc_counters->metro_swap_prob;

	// N_replica_pt - 1 is the total number of swaps
	int const num_swaps = ((param->d_N_replica_pt) - 1);

	// set all probabilities to 0
	for(int a = 0; a < num_swaps; a++)
		aux_p[a] = 0.0;

	int const is_even = num_swaps % 2;              // check if num_swaps is even or not
	int const num_even = (num_swaps + is_even) / 2; // number of swaps for even replica
	int const num_odd = (num_swaps - is_even) / 2;  // number of swaps for odd replica

	// to be sure detailed balance is satisfied, choose randomly whether to swap first odd or even copies
	int is_even_first, num_swaps_1, num_swaps_2;
	if(casuale() < 0.5) // first swap all even copies, then all odd copies
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

	// swaps are done for all couples (mu,j) where mu=defect_dir and j != mu => three couples
	int const mu = param->d_defect_dir;

	// first group of swaps

	// compute action differences (multicanonical contribution in the next loop)
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:aux_p[:num_swaps])
	#endif
	for(long s = 0; s < num_swaps_1 * (swap_rectangle->d_vol_rect); s++)
		{
		long const n = s % (swap_rectangle->d_vol_rect);              // action changes only in the first neighborhood of the defect,
		long const r = swap_rectangle->rect_sites[n];                 // having also swapped the twist factors
		int const k = (int) ((s - n) / (swap_rectangle->d_vol_rect)); // index of replica pair to swap
		int const a = 2 * k + is_even_first;                          // label of 1st replica
		int const b = a + 1;                                          // label of 2nd replica

		for(int i = 0; i < STDIM - 1; i++)
			{
			int const j = geo->d_orth_dir[mu][i];
			// contribution to action difference between replicas a and b of site r on plane (mu,j)
			aux_p[a] += delta_action_swap(GC, geo, param, r, mu, j, a, b);
			}
		}

	// do the swaps
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int k = 0; k < num_swaps_1; k++)
		{
		int const a = 2 * k + is_even_first;
		int const b = a + 1;

		// multicanonical contribution to the action difference
		#ifdef MULTICANONICAL_MODE
		aux_p[a] += delta_topo_potential_swap(GC, a, b, param);
		#endif
		aux_p[a] = exp(-aux_p[a]);                                // metropolis swap probability = exp{ - (swapped action - unswapped action) }
		metropolis_single_swap(GC, a, b, aux_p[a], acc_counters); // metropolis step
		}


	// second group of swaps: reverse parity of replica pairs to swap
	is_even_first = 1 - is_even_first; //

	// compute action differences
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:aux_p[:num_swaps])
	#endif
	for(long s = 0; s < num_swaps_2 * (swap_rectangle->d_vol_rect); s++)
		{
		long const n = s % (swap_rectangle->d_vol_rect);
		long const r = swap_rectangle->rect_sites[n];
		int const k = (int) ((s - n) / (swap_rectangle->d_vol_rect));
		int const a = 2 * k + is_even_first;
		int const b = a + 1;

		for(int i = 0; i < STDIM - 1; i++)
			{
			int const j = geo->d_orth_dir[mu][i];
			aux_p[a] += delta_action_swap(GC, geo, param, r, mu, j, a, b);
			}
		}

	// do the swaps
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int k = 0; k < num_swaps_2; k++)
		{
		int const a = 2 * k + is_even_first;
		int const b = a + 1;

		#ifdef MULTICANONICAL_MODE
		aux_p[a] += delta_topo_potential_swap(GC, a, b, param);
		#endif
		aux_p[a] = exp(-aux_p[a]);
		metropolis_single_swap(GC, a, b, aux_p[a], acc_counters);
		}
	}


double delta_action_swap(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param,
                         long const r, int const i, int const j, int const a, int const b)
	{
	// plaquettes, including twist factors in function plaquettep
	double const re_tr_plaq_a = plaquettep(&(GC[a]), geo, param, r, i, j); // (Re Tr plaq_a(r,i,j) )/N_c , replica label = a
	double const re_tr_plaq_b = plaquettep(&(GC[b]), geo, param, r, i, j); // (Re Tr plaq_b(r,i,j) )/N_c , replica label = b

	// boundary conditions
	long const rpi = nnp(geo, r, i);
	long const rpj = nnp(geo, r, j);
	double const K_a = (GC[a].C[r][i]) * (GC[a].C[rpi][j]) * (GC[a].C[rpj][i]) * (GC[a].C[r][j]);
	double const K_b = (GC[b].C[r][i]) * (GC[b].C[rpi][j]) * (GC[b].C[rpj][i]) * (GC[b].C[r][j]);

	// (swapped action - unswapped action) = beta * delta_K * delta_plaq (twist factors swapped, otherwise delta_K -> K_a*Z_a-K_b*Z_b )
	return param->d_beta * (K_a - K_b) * (re_tr_plaq_a - re_tr_plaq_b);
	}

// swaps are serial, evaluation of swap probability is parallelized (use this version of 'swap' if gcc_version < 6.0 or icc_version < 14.0)
/*
void swap(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
				 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)
  {
	int aux_i, i, j, num_swaps, a, b;
	long n;
	double metro_swap_prob = 0.0;

	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };

	// N_replica_pt - 1 is the total number of swaps
	num_swaps = ((param->d_N_replica_pt)-1);

	// swaps are done for all couples (i,j) where i=defect_dir and j !=i => three couples
	i=param->d_defect_dir;

	for(a=0;a<num_swaps;a++)
		{
		b=a+1;
		metro_swap_prob = 0.0;
		// compute action difference between replica a and b
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) reduction(+:metro_swap_prob) private(n,aux_i,j)
		#endif
		for(n=0;n<(swap_rectangle->d_vol_rect);n++)
			{
			long r = swap_rectangle->rect_sites[n]; // action changes only in the first neighborhood of the defect
			for(aux_i=0; aux_i<STDIM-1; aux_i++)
				{
				j = perp_dir[param->d_defect_dir][aux_i];
				metro_swap_prob += delta_action_swap(GC, geo, param, r, i, j, a, b);
				}
			}

		// do the swap
		metro_swap_prob=exp(-metro_swap_prob); // metropolis swap probability
		metropolis_single_swap(GC,a,b,metro_swap_prob,acc_counters);
		}
	}
*/


// metropolis step to swap replica a and b with probability p, including the twist factors
void metropolis_single_swap(Gauge_Conf *const GC, int const a, int const b, double const p, Acc_Utils *acc_counters)
	{
	// increase counter of tried swaps for replicas (a, a+1)
	acc_counters->num_swap[a]++;

	// Metropolis test: p<1 => acc=1 with probability p, p>=1 => acc=1 (already assigned)
	int acc = 1;
	if(p < 1 && casuale() > p) acc = 0;

	// if Metropolis is accepted, swap replicas, including the twist factors and the stored charges for multicanonic
	if(acc == 1)
		{
		// swap of configurations
		GAUGE_GROUP **aux = GC[a].lattice;
		GC[a].lattice = GC[b].lattice;
		GC[b].lattice = aux;

		double complex **aux_Z = GC[a].Z;
		GC[a].Z = GC[b].Z;
		GC[b].Z = aux_Z;

		// swap of auxiliary configurations
		aux = GC[a].lattice_copy;
		GC[a].lattice_copy = GC[b].lattice_copy;
		GC[b].lattice_copy = aux;

		aux_Z = GC[a].Z_copy;
		GC[a].Z_copy = GC[b].Z_copy;
		GC[b].Z_copy = aux_Z;

		// swap of multicanonic utils
		#ifdef MULTICANONICAL_MODE
		double aux_charge;

		aux_charge = GC[a].stored_topcharge;
		GC[a].stored_topcharge = GC[b].stored_topcharge;
		GC[b].stored_topcharge = aux_charge;

		aux = GC[a].lattice_cold;
		GC[a].lattice_cold = GC[b].lattice_cold;
		GC[b].lattice_cold = aux;

		aux = GC[a].lattice_copy_cold;
		GC[a].lattice_copy_cold = GC[b].lattice_copy_cold;
		GC[b].lattice_copy_cold = aux;
		#endif

		// swap of labels
		int const aux_label = GC[a].conf_label;
		GC[a].conf_label = GC[b].conf_label;
		GC[b].conf_label = aux_label;

		// swap of translation tracking
		for(int i = 0; i < STDIM; i++)
			{
			long const aux_translation = GC[a].translation[i];
			GC[a].translation[i] = GC[b].translation[i];
			GC[b].translation[i] = aux_translation;
			}

		// increase counter of successful swaps for replicas (a, a+1)
		acc_counters->num_accepted_swap[a]++;
		}
	}


// translation of one lattice spacing of the configuration, including the Z factors
// direction is chosen randomly, verse is always positive
void conf_translation(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	// extract random direction
	int const dir = (int) (STDIM * casuale());
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


void init_acc_utils(Acc_Utils *acc_counters, GParam const *const param)
	{
	if(param->d_N_replica_pt == 1)
		{
		acc_counters->num_accepted_swap = NULL;
		acc_counters->num_swap = NULL;
		acc_counters->metro_swap_prob = NULL;
		}
	else
		{
		allocate_array_long(&(acc_counters->num_accepted_swap), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		allocate_array_long(&(acc_counters->num_swap), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		allocate_array_double(&(acc_counters->metro_swap_prob), param->d_N_replica_pt - 1, __FILE__, __LINE__);
		for(int i = 0; i < (param->d_N_replica_pt - 1); i++)
			{
			acc_counters->num_accepted_swap[i] = 0;
			acc_counters->num_swap[i] = 0;
			}
		}
	#ifdef MULTICANONICAL_MODE
	init_multicanonic_acc_utils(acc_counters, param);
	#endif
	}


void free_acc_utils(Acc_Utils *acc_counters, GParam const *const param)
	{
	if(param->d_N_replica_pt > 1)
		{
		free(acc_counters->num_accepted_swap);
		free(acc_counters->num_swap);
		free(acc_counters->metro_swap_prob);
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
			acc = (double) acc_counters->num_accepted_swap[k] / (double) acc_counters->num_swap[k];
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


void print_conf_labels(FILE *fp, Gauge_Conf const *const GC, GParam const *const param)
	{
	if(param->d_N_replica_pt == 1)
		{
		(void) fp; // to suppress compiler warning of unused variable
		(void) GC; // to suppress compiler warning of unused variable
		return;
		}
	fprintf(fp, "%9ld ", GC[0].update_index);
	for(int r = 0; r < param->d_N_replica_pt; r++)
		{
		fprintf(fp, "%4d ", GC[r].conf_label);
		#ifdef MULTICANONICAL_MODE
		fprintf(fp, "% 9.6f ", GC[r].stored_topcharge);
		#endif
		}
	fprintf(fp, "\n");
	fflush(fp);
	}


#endif
