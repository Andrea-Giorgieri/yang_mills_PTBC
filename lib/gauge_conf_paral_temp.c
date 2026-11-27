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

#include "../include/device.h"
#include "../include/device_data.h"
#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/memalign.h"
#include "../include/random.h"

// conf_translation copies the links as arrays of doubles
_Static_assert(sizeof(GAUGE_GROUP) % sizeof(double) == 0, "GAUGE_GROUP must be made of doubles");

// sums P1 and P2 of the configuration GC entering the swap action differences (see swap_group), on the device of GC:
//   P_n = sum_{r in swap_rect, j != mu, n links of the plaquette (r,mu,j) on the defect} ReTr(Z plaq(r,mu,j))/N_c
// The defect coefficient of a plaquette with n links on the defect is K_c = c^n (n = 0, 1, 2), so
//   sum_{r in swap_rect, j != mu} [K_{c_a}(r,mu,j) - K_{c_b}(r,mu,j)] ReTr(Z plaq(r,mu,j))/N_c = (c_a - c_b) P1 + (c_a^2 - c_b^2) P2
// for every pair of coefficients: P1 and P2 depend only on the configuration, not on its replica index
void swap_plaquette_sums(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param,
                         Rectangle const *const swap_rectangle, double *const P1, double *const P2)
	{
	int const mu = param->d_defect_dir; // swaps change only the plaquettes (mu, j), j != mu, around the defect
	long const vol_rect = swap_rectangle->d_vol_rect;
	double p1 = 0.0, p2 = 0.0;
	int const dev = gc_on_device(GC);

	// reduction on the host scalars: queue of the device (after the update), then wait in this function
	#ifdef OPENACC_MODE
	#pragma acc parallel loop gang vector vector_length(ACC_VL_MEAS) if(dev) reduction(+:p1, p2) present(GC[0:1], geo[0:1], param[0:1], swap_rectangle[0:1]) async(QUEUE_DEVICE)
	#endif
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:p1, p2)
	#endif
	for(long n = 0; n < vol_rect; n++)
		{
		long const r = swap_rectangle->rect_sites[n];
		long const rpmu = nnp(geo, r, mu);
		for(int i = 0; i < STDIM - 1; i++)
			{
			int const j = geo->d_orth_dir[mu][i];
			int const n_def = plaquette_defect_links(geo, param, r, rpmu, nnp(geo, r, j), mu, j);
			if(n_def == 1)
				p1 += plaquettep(GC, geo, param, r, mu, j);
			else if(n_def == 2)
				p2 += plaquettep(GC, geo, param, r, mu, j);
			}
		}
	wait_device_queue(dev);

	*P1 = p1;
	*P2 = p2;
	}


// sites of the swap rectangle summed by one gang in swap_plaquette_sums_replica
#define SWAP_SUM_CHUNK 1024

// sums P1[g], P2[g] of swap_plaquette_sums for all the configurations GC[g], 0 <= g < N_replica_pt.
// On the devices a single reduction per configuration would cost a launch and a wait each (the host thread drives all
// the devices in turn): instead every device runs one kernel over all its slots, in which each gang sums a chunk of
// SWAP_SUM_CHUNK sites of one slot into acc_counters->swap_partial, copied back in the same queue; the kernels and the
// copies of all the devices are queued before waiting once per device, then the chunks are added on the host in a
// fixed order (same result in every run, equal to the single reduction up to rounding).
// On the host (no device data): swap_plaquette_sums for every configuration.
static void swap_plaquette_sums_replica(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param,
                                        Rectangle const *const swap_rectangle, Acc_Utils *const acc_counters,
                                        double *const P1, double *const P2)
	{
	int const Nr = param->d_N_replica_pt;

	#ifdef OPENACC_MODE
	if(gc_on_device(&(GC[0])))
		{
		int const mu = param->d_defect_dir; // swaps change only the plaquettes (mu, j), j != mu, around the defect
		long const vol_rect = swap_rectangle->d_vol_rect;
		long const n_chunk = (vol_rect + SWAP_SUM_CHUNK - 1) / SWAP_SUM_CHUNK;
		long const n_part = 2 * (long) Nr * n_chunk; // part[2 * q], part[2 * q + 1]: P1, P2 of chunk q = g * n_chunk + c

		if(acc_counters->swap_partial_size != n_part)
			{
			free_double_on_devices(acc_counters->swap_partial, acc_counters->swap_partial_size);
			acc_counters->swap_partial = alloc_double_on_devices(n_part);
			acc_counters->swap_partial_size = n_part;
			}
		double *const part = acc_counters->swap_partial;

		FOR_EACH_DEVICE(d)
			{
			select_device(d);
			long const first = device_slot_begin(d);
			long const n_sl = device_slot_end(d) - first;
			long const q_begin = first * n_chunk;
			long const n_q = n_sl * n_chunk;
			if(n_q == 0) continue;

			// one gang per chunk, the vector lanes of the gang reduce its sites
			#pragma acc parallel loop gang vector_length(ACC_VL_MEAS) present(GC[first:n_sl], geo[0:1], param[0:1], swap_rectangle[0:1], part[2 * q_begin : 2 * n_q]) async(QUEUE_DEVICE)
			for(long q = q_begin; q < q_begin + n_q; q++)
				{
				long const g = q / n_chunk;
				long const begin = (q - g * n_chunk) * SWAP_SUM_CHUNK;
				long const end = (begin + SWAP_SUM_CHUNK < vol_rect) ? begin + SWAP_SUM_CHUNK : vol_rect;
				double p1 = 0.0, p2 = 0.0;

				#pragma acc loop vector reduction(+:p1, p2)
				for(long n = begin; n < end; n++)
					{
					long const r = swap_rectangle->rect_sites[n];
					long const rpmu = nnp(geo, r, mu);
					for(int i = 0; i < STDIM - 1; i++)
						{
						int const j = geo->d_orth_dir[mu][i];
						int const n_def = plaquette_defect_links(geo, param, r, rpmu, nnp(geo, r, j), mu, j);
						if(n_def == 1)
							p1 += plaquettep(&(GC[g]), geo, param, r, mu, j);
						else if(n_def == 2)
							p2 += plaquettep(&(GC[g]), geo, param, r, mu, j);
						}
					}
				part[2 * q] = p1;
				part[2 * q + 1] = p2;
				}

			#pragma acc update self(part[2 * q_begin : 2 * n_q]) async(QUEUE_DEVICE)
			}

		FOR_EACH_DEVICE(d)
			{
			select_device(d);
			wait_device_queue(1);
			}

		for(int g = 0; g < Nr; g++)
			{
			P1[g] = 0.0;
			P2[g] = 0.0;
			for(long c = 0; c < n_chunk; c++)
				{
				P1[g] += part[2 * (g * n_chunk + c)];
				P2[g] += part[2 * (g * n_chunk + c) + 1];
				}
			}
		return; // sums done: the loop below is only for configurations without device data
		}
	#else
	(void) acc_counters;
	#endif

	// reached without OPENACC_MODE, or with OPENACC_MODE when the configurations have no device data (init_device_data
	// not called, or force_host_compute): every swap_plaquette_sums then runs on the host
	for(int g = 0; g < Nr; g++)
		swap_plaquette_sums(&(GC[g]), geo, param, swap_rectangle, &(P1[g]), &(P2[g]));
	}


// one group of n_pairs swaps: pairs (a, a+1) with a = first, first + 2, ...
//
// Swapping replicas a and b = a+1 exchanges the replica indices of A = GC[replica_to_gc[a]] and B = GC[replica_to_gc[b]]
// The action difference
//   swapped - unswapped = beta * (D_A - D_B),
//   D_X = sum_{r in swap_rect, j != mu} [K_{c_a}(r,mu,j) - K_{c_b}(r,mu,j)] ReTr(Z plaq_X(r,mu,j))/N_c
//       = (c_a - c_b) P1_X + (c_a^2 - c_b^2) P2_X,
// with P1[g], P2[g] the sums of the configuration in slot g (swap_plaquette_sums): the configurations do not change
// during the swaps, so the same sums serve both groups
static void swap_group(Gauge_Conf *const GC, GParam const *const param, Acc_Utils *acc_counters,
                       int const first, int const n_pairs, double const *const P1, double const *const P2)
	{
	int const *const replica_to_gc = acc_counters->replica_to_gc;
	double *const aux_p = acc_counters->metro_swap_prob; // alias for the Metropolis probabilities of the pairs

	// Metropolis tests and swaps (serial with Philox: the Metropolis tests use the single host generator)
	#if defined(OPENMP_MODE) && !defined(RNG_PHILOX_MODE)
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int k = 0; k < n_pairs; k++)
		{
		int const a = 2 * k + first;
		int const b = a + 1;
		int const gc_a = replica_to_gc[a]; // configuration A, with replica index a
		int const gc_b = replica_to_gc[b]; // configuration B, with replica index b

		double const c_a = param->d_pt_bound_cond_coeff[a];
		double const c_b = param->d_pt_bound_cond_coeff[b];

		// swapped action - unswapped action = beta * (D_A - D_B)
		aux_p[a] = param->d_beta * ((c_a - c_b) * (P1[gc_a] - P1[gc_b]) + (c_a * c_a - c_b * c_b) * (P2[gc_a] - P2[gc_b]));
		#ifdef MULTICANONICAL_MODE
		aux_p[a] += delta_topo_potential_swap(&(GC[gc_a]), &(GC[gc_b]), param);
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

	if(num_swaps > 0)
		{
		// sums of the swap action differences, one pair per configuration (slot), on the device of the slot
		double *P1, *P2;
		allocate_array_double(&P1, param->d_N_replica_pt, __FILE__, __LINE__);
		allocate_array_double(&P2, param->d_N_replica_pt, __FILE__, __LINE__);
		swap_plaquette_sums_replica(GC, geo, param, swap_rectangle, acc_counters, P1, P2);

		swap_group(GC, param, acc_counters, is_even_first, num_swaps_1, P1, P2);     // first group
		swap_group(GC, param, acc_counters, 1 - is_even_first, num_swaps_2, P1, P2); // second group: reverse parity

		free(P1);
		free(P2);
		}

	// the defect coefficients read on the devices depend on replica_index
	copy_replica_index_to_device(GC, param);
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
// on the device holding GC if gc_on_device (kernels in its queue); the translation counters stay on the host
void conf_translation(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	// extract random direction
	int const dir = (int) (STDIM * casuale(rng_host()));
	GC->translation[dir] += 1;

	int const dev = gc_on_device(GC);
	(void) dev;

	// The fields are contiguous blocks [volume][row] (lattice: STDIM links, Z: n_planes + 1 factors per site) and the
	// translation only moves whole rows: both loops below copy one double per iteration, so that consecutive threads
	// access consecutive addresses (a thread per link copying a whole matrix makes strided accesses and reaches a small
	// fraction of the memory bandwidth). Segments of the iteration index: lattice, [cold lattice,] Z.
	long const V = param->d_volume;
	long const w_lat = (long) (STDIM * sizeof(GAUGE_GROUP) / sizeof(double)); // doubles per site of a lattice
	long const w_Z = 2 * (long) (param->d_n_planes + 1);                       // doubles per site of Z
	long const n_lat = V * w_lat;
	#ifdef MULTICANONICAL_MODE
	long const n_lat_all = 2 * n_lat; // lattice and lattice_cold
	#else
	long const n_lat_all = n_lat;
	#endif

	// translation in direction +dir, including the Z factors and cold lattice (read from the auxiliary copies, equal
	// to the accepted configuration): row r = row nnm(r, dir) of the copy
	#ifdef OPENACC_MODE
	#pragma acc parallel loop gang vector if(dev) present(GC[0:1], geo[0:1]) firstprivate(dir, n_lat, n_lat_all, w_lat, w_Z) async(QUEUE_DEVICE)
	#endif
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long t = 0; t < n_lat_all + V * w_Z; t++)
		{
		double *dst;
		double const *src;
		long u, w;
		if(t < n_lat)
			{
			u = t;
			w = w_lat;
			dst = (double *) &(GC->lattice[0][0]);
			src = (double const *) &(GC->lattice_copy[0][0]);
			}
		#ifdef MULTICANONICAL_MODE
		else if(t < n_lat_all)
			{
			u = t - n_lat;
			w = w_lat;
			dst = (double *) &(GC->lattice_cold[0][0]);
			src = (double const *) &(GC->lattice_copy_cold[0][0]);
			}
		#endif
		else
			{
			u = t - n_lat_all;
			w = w_Z;
			dst = (double *) &(GC->Z[0][0]);
			src = (double const *) &(GC->Z_copy[0][0]);
			}
		long const r = u / w;
		dst[u] = src[nnm(geo, r, dir) * w + (u - r * w)];
		}

	// swap translated cold lattice with previous cold lattice
	#ifdef MULTICANONICAL_MODE
	swap_lattice_cold(GC);
	#endif

	// update the auxiliary confs: lattice_copy = lattice, Z_copy = Z
	#ifdef OPENACC_MODE
	#pragma acc parallel loop gang vector if(dev) present(GC[0:1]) firstprivate(n_lat) async(QUEUE_DEVICE)
	#endif
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long t = 0; t < n_lat + V * w_Z; t++)
		{
		if(t < n_lat)
			((double *) &(GC->lattice_copy[0][0]))[t] = ((double const *) &(GC->lattice[0][0]))[t];
		else
			((double *) &(GC->Z_copy[0][0]))[t - n_lat] = ((double const *) &(GC->Z[0][0]))[t - n_lat];
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

	// allocated at the first swap on the devices (its size depends on the swap rectangle)
	acc_counters->swap_partial = NULL;
	acc_counters->swap_partial_size = 0;

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
	free_double_on_devices(acc_counters->swap_partial, acc_counters->swap_partial_size);

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
