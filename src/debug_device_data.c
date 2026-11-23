#ifndef DEBUG_DEVICE_DATA_C
#define DEBUG_DEVICE_DATA_C

#include "../include/macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OPENACC_MODE
#include <openacc.h>
#endif

#include "../include/device.h"
#include "../include/device_data.h"
#include "../include/gauge_conf.h"
#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"

// Test of the data on the devices (init_device_data and copy functions): compute regions read the shared data and
// the configurations through the device pointer tables and the results are compared with the host.
// Without OpenACC the compute regions run on the host and the test checks only its own logic.

static int n_fail = 0;
static int n_test = 0;


static void check(int const ok, char const *const what, double const device_value, double const host_value)
	{
	n_test++;
	if(!ok) n_fail++;
	printf("%-48s device % .15e  host % .15e  %s\n", what, device_value, host_value, ok ? "OK" : "FAILED");
	}


static int close_enough(double const a, double const b)
	{
	return fabs(a - b) <= 1.0e-10 * (1.0 + fabs(b));
	}


// weight of the link (r, i) in the checksums
#define LINK_WEIGHT(r, i) ((double) (1 + ((r) * STDIM + (i)) % 7))


static double lattice_sum_host(GAUGE_GROUP **lat, long const V)
	{
	double s = 0.0;
	for(long r = 0; r < V; r++)
		for(int i = 0; i < STDIM; i++)
			s += retr(&(lat[r][i])) * LINK_WEIGHT(r, i);
	return s;
	}


// which: 0 lattice, 1 lattice_copy, 2 lattice_cold, 3 lattice_copy_cold (current device = device of gc)
static double lattice_sum_device(Gauge_Conf *gc, int const which, long const V)
	{
	double s = 0.0;
	#ifdef OPENACC_MODE
	#pragma acc parallel loop reduction(+:s) present(gc[0:1])
	#endif
	for(long r = 0; r < V; r++)
		{
		GAUGE_GROUP **lat = (which == 0) ? gc->lattice : (which == 1) ? gc->lattice_copy : (which == 2) ? gc->lattice_cold : gc->lattice_copy_cold;
		for(int i = 0; i < STDIM; i++)
			s += retr(&(lat[r][i])) * LINK_WEIGHT(r, i);
		}
	return s;
	}


// 1 if the two lattices have the same matrix elements (not memcmp: the padding of GAUGE_GROUP is not initialized)
static int same_lattice(GAUGE_GROUP **a, GAUGE_GROUP **b, long const V)
	{
	for(long r = 0; r < V; r++)
		for(int i = 0; i < STDIM; i++)
			{
			GAUGE_GROUP diff;
			equal(&diff, &(a[r][i]));
			minus_equal(&diff, &(b[r][i]));
			if(norm(&diff) != 0.0) return 0;
			}
	return 1;
	}


static GAUGE_GROUP **lattice_of(Gauge_Conf *gc, int const which)
	{
	return (which == 0) ? gc->lattice : (which == 1) ? gc->lattice_copy : (which == 2) ? gc->lattice_cold : gc->lattice_copy_cold;
	}


static double Z_sum_host(double complex **Z, long const V, int const n)
	{
	double s = 0.0;
	for(long r = 0; r < V; r++)
		for(int j = 0; j < n; j++)
			s += (creal(Z[r][j]) + 2.0 * cimag(Z[r][j])) * LINK_WEIGHT(r, j);
	return s;
	}


// which: 0 Z, 1 Z_copy
static double Z_sum_device(Gauge_Conf *gc, int const which, long const V, int const n)
	{
	double s = 0.0;
	#ifdef OPENACC_MODE
	#pragma acc parallel loop reduction(+:s) present(gc[0:1])
	#endif
	for(long r = 0; r < V; r++)
		{
		for(int j = 0; j < n; j++)
			{
			double complex const z = (which == 0) ? gc->Z[r][j] : gc->Z_copy[r][j];
			s += (creal(z) + 2.0 * cimag(z)) * LINK_WEIGHT(r, j);
			}
		}
	return s;
	}


static void test_geometry_data(Geometry *geo, GParam *param)
	{
	long const V = param->d_volume;
	char what[STD_STRING_LENGTH];

	FOR_EACH_DEVICE(d)
		{
		select_device(d);

		// nearest neighbours and mu-components through the 2d/3d tables, exact integer counts
		long ok_dev = 0, sum_dev = 0, mask_dev = 0;
		#ifdef OPENACC_MODE
		#pragma acc parallel loop reduction(+:ok_dev, sum_dev, mask_dev) present(geo[0:1])
		#endif
		for(long r = 0; r < V; r++)
			{
			for(int mu = 0; mu < STDIM; mu++)
				{
				ok_dev += (geo->d_nnm[geo->d_nnp[r][mu]][mu] == r) ? 1 : 0;
				ok_dev += (geo->d_mutimespace[mu][geo->d_mucomp[mu][r]][geo->d_muorth[mu][r]] == r) ? 1 : 0;
				}
			sum_dev += geo->d_sip_to_si[r];
			mask_dev += geo->d_defect_mask[r];
			}

		long sum_host = 0, mask_host = 0;
		for(long r = 0; r < V; r++)
			{
			sum_host += geo->d_sip_to_si[r];
			mask_host += geo->d_defect_mask[r];
			}

		sprintf(what, "device %d: nnp/nnm/mutimespace", d);
		check(ok_dev == 2 * STDIM * V, what, (double) ok_dev, (double) (2 * STDIM * V));
		sprintf(what, "device %d: sip_to_si", d);
		check(sum_dev == sum_host, what, (double) sum_dev, (double) sum_host);
		sprintf(what, "device %d: defect_mask", d);
		check(mask_dev == mask_host, what, (double) mask_dev, (double) mask_host);

		// global lookup tables
		#if STDIM == 4
		long tab_dev = 0, tab_host = 0;
		#ifdef OPENACC_MODE
		#pragma acc serial copy(tab_dev)
		#endif
			{
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					tab_dev += (i * 4 + j + 1) * (g_signed_ord_orth_dir[i][j][0] + 7 * g_signed_ord_orth_dir[i][j][1]);
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 3; j++)
					tab_dev += (i * 3 + j + 11) * g_indep_perm_dir[i][j];
			}
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++)
				tab_host += (i * 4 + j + 1) * (g_signed_ord_orth_dir[i][j][0] + 7 * g_signed_ord_orth_dir[i][j][1]);
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 3; j++)
				tab_host += (i * 3 + j + 11) * g_indep_perm_dir[i][j];
		sprintf(what, "device %d: global lookup tables", d);
		check(tab_dev == tab_host, what, (double) tab_dev, (double) tab_host);
		#endif
		}
	}


static long rect_sum_host(Rectangle const *rect)
	{
	long s = 0;
	for(long n = 0; n < rect->d_vol_rect; n++)
		s += (n % 3 + 1) * rect->rect_sites[n];
	return s;
	}


// kind: 0 swap_rect, 1 update_rect[l], 2 clover_rect[l], 3 topcharge_rect[l], 4 cooling_rect[l][j]
static long rect_sum_device(Rect_Utils *rect_aux, int const kind, int const l, int const j)
	{
	long s = 0;
	#ifdef OPENACC_MODE
	#pragma acc serial copy(s) present(rect_aux[0:1])
	#endif
		{
		Rectangle const *rect = (kind == 0) ? &(rect_aux->swap_rect) : (kind == 1) ? &(rect_aux->update_rect[l]) :
		                        (kind == 2) ? &(rect_aux->clover_rect[l]) : (kind == 3) ? &(rect_aux->topcharge_rect[l]) :
		                        &(rect_aux->cooling_rect[l][j]);
		for(long n = 0; n < rect->d_vol_rect; n++)
			s += (n % 3 + 1) * rect->rect_sites[n];
		}
	return s;
	}


static void test_rect_data(Rect_Utils *rect_aux, GParam *param)
	{
	char what[STD_STRING_LENGTH];

	FOR_EACH_DEVICE(d)
		{
		select_device(d);

		long s_dev = rect_sum_device(rect_aux, 0, 0, 0), s_host = rect_sum_host(&(rect_aux->swap_rect));
		sprintf(what, "device %d: swap_rect", d);
		check(s_dev == s_host, what, (double) s_dev, (double) s_host);

		for(int l = 0; l < param->d_N_hierarc_levels; l++)
			{
			s_dev = rect_sum_device(rect_aux, 1, l, 0);
			s_host = rect_sum_host(&(rect_aux->update_rect[l]));
			sprintf(what, "device %d: update_rect[%d]", d, l);
			check(s_dev == s_host, what, (double) s_dev, (double) s_host);
			#ifdef THETA_MODE
			s_dev = rect_sum_device(rect_aux, 2, l, 0);
			s_host = rect_sum_host(&(rect_aux->clover_rect[l]));
			sprintf(what, "device %d: clover_rect[%d]", d, l);
			check(s_dev == s_host, what, (double) s_dev, (double) s_host);
			#endif
			#ifdef MULTICANONICAL_MODE
			s_dev = rect_sum_device(rect_aux, 3, l, 0);
			s_host = rect_sum_host(&(rect_aux->topcharge_rect[l]));
			sprintf(what, "device %d: topcharge_rect[%d]", d, l);
			check(s_dev == s_host, what, (double) s_dev, (double) s_host);
			if(rect_aux->cooling_rect != NULL)
				{
				for(int j = 0; j < param->d_topo_coolsteps; j++)
					{
					s_dev = rect_sum_device(rect_aux, 4, l, j);
					s_host = rect_sum_host(&(rect_aux->cooling_rect[l][j]));
					sprintf(what, "device %d: cooling_rect[%d][%d]", d, l, j);
					check(s_dev == s_host, what, (double) s_dev, (double) s_host);
					}
				}
			#endif
			}
		}
	}


// replica_defect_coeff of gc read on its device
static double defect_coeff_device(Gauge_Conf *gc, GParam *param)
	{
	double c = 0.0;
	#ifdef OPENACC_MODE
	#pragma acc serial copy(c) present(gc[0:1], param[0:1])
	#endif
	c = replica_defect_coeff(gc, param);
	return c;
	}


static void test_param_data(Gauge_Conf *GC, GParam *param)
	{
	char what[STD_STRING_LENGTH];
	int const Nr = param->d_N_replica_pt;

	// defect coefficients before and after a permutation of the replica indices
	for(int pass = 0; pass < 2; pass++)
		{
		if(pass == 1)
			{
			for(int g = 0; g < Nr; g++)
				GC[g].replica_index = (GC[g].replica_index + 1) % Nr;
			copy_replica_index_to_device(GC, param);
			}
		for(int g = 0; g < Nr; g++)
			{
			select_device(GC[g].device);
			double const c_dev = defect_coeff_device(&(GC[g]), param);
			double const c_host = replica_defect_coeff(&(GC[g]), param);
			sprintf(what, "slot %d: defect coeff%s", g, (pass == 1) ? " (after permutation)" : "");
			check(c_dev == c_host, what, c_dev, c_host);
			}
		}
	// restore the replica indices
	for(int g = 0; g < Nr; g++)
		GC[g].replica_index = (GC[g].replica_index - 1 + Nr) % Nr;
	copy_replica_index_to_device(GC, param);

	#ifdef MULTICANONICAL_MODE
	// topological potential, before and after a change on the host
	for(int pass = 0; pass < 2; pass++)
		{
		if(pass == 1)
			{
			for(int a = 0; a < Nr; a++)
				for(int i = 0; i < param->d_n_grid; i++)
					param->d_grid[a][i] += 1.0 + 0.001 * (double) a;
			copy_topo_grid_to_device(param);
			}
		double s_host = 0.0;
		for(int a = 0; a < Nr; a++)
			for(int i = 0; i < param->d_n_grid; i++)
				s_host += param->d_grid[a][i] * (double) (1 + (a + i) % 5);
		FOR_EACH_DEVICE(d)
			{
			select_device(d);
			double s_dev = 0.0;
			int const n_grid = param->d_n_grid;
			#ifdef OPENACC_MODE
			#pragma acc serial copy(s_dev) present(param[0:1])
			#endif
				{
				for(int a = 0; a < Nr; a++)
					for(int i = 0; i < n_grid; i++)
						s_dev += param->d_grid[a][i] * (double) (1 + (a + i) % 5);
				}
			// not bitwise: on the GPU the products and sums are contracted into fused multiply-adds
			sprintf(what, "device %d: topo grid%s", d, (pass == 1) ? " (after change)" : "");
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			}
		}
	for(int a = 0; a < Nr; a++)
		for(int i = 0; i < param->d_n_grid; i++)
			param->d_grid[a][i] -= 1.0 + 0.001 * (double) a;
	copy_topo_grid_to_device(param);
	#endif
	}


static void test_gauge_conf_data(Gauge_Conf *GC, GParam *param, Meas_Utils *meas_aux)
	{
	char what[STD_STRING_LENGTH];
	long const V = param->d_volume;
	int const nZ = param->d_n_planes + 1;

	for(int g = 0; g < param->d_N_replica_pt; g++)
		{
		Gauge_Conf *gc = &(GC[g]);
		select_device(gc->device);

		// fields copied by init_device_data
		#ifdef MULTICANONICAL_MODE
		int const n_lat_all = 4; // lattice, lattice_copy, lattice_cold, lattice_copy_cold
		#else
		int const n_lat_all = 2; // lattice, lattice_copy
		#endif
		char const *const lat_name[4] = {"lattice", "lattice_copy", "lattice_cold", "lattice_copy_cold"};
		for(int w = 0; w < n_lat_all; w++)
			{
			double const s_dev = lattice_sum_device(gc, w, V);
			double const s_host = lattice_sum_host(lattice_of(gc, w), V);
			sprintf(what, "slot %d (device %d): %s", g, gc->device, lat_name[w]);
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			}
		for(int w = 0; w < 2; w++)
			{
			double const s_dev = Z_sum_device(gc, w, V, nZ);
			double const s_host = Z_sum_host((w == 0) ? gc->Z : gc->Z_copy, V, nZ);
			sprintf(what, "slot %d (device %d): %s", g, gc->device, (w == 0) ? "Z" : "Z_copy");
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			}
		#ifdef RNG_PHILOX_MODE
		long r_dev = 0, r_host = 0;
		#ifdef OPENACC_MODE
		#pragma acc parallel loop reduction(+:r_dev) present(gc[0:1])
		#endif
		for(long r = 0; r < V; r++)
			r_dev += (long) gc->rng_site[r].ctr[0] + (long) gc->rng_site[r].ctr[2] + 3L * (long) gc->rng_site[r].ctr[3] + (long) gc->rng_site[r].key[0];
		for(long r = 0; r < V; r++)
			r_host += (long) gc->rng_site[r].ctr[0] + (long) gc->rng_site[r].ctr[2] + 3L * (long) gc->rng_site[r].ctr[3] + (long) gc->rng_site[r].key[0];
		sprintf(what, "slot %d (device %d): rng_site", g, gc->device);
		check(r_dev == r_host, what, (double) r_dev, (double) r_host);
		#endif

		// write on the device and copy back: lattice_copy = lattice (host lattice_copy zeroed before)
		for(long r = 0; r < V; r++)
			for(int i = 0; i < STDIM; i++)
				zero(&(gc->lattice_copy[r][i]));
		copy_gauge_conf_to_device(gc, param, GC_COPY);
		#ifdef OPENACC_MODE
		#pragma acc parallel loop collapse(2) present(gc[0:1])
		#endif
		for(long r = 0; r < V; r++)
			for(int i = 0; i < STDIM; i++)
				equal(&(gc->lattice_copy[r][i]), &(gc->lattice[r][i]));
		copy_gauge_conf_to_host(gc, param, GC_COPY);
		int const same_copy = same_lattice(gc->lattice_copy, gc->lattice, V);
		sprintf(what, "slot %d (device %d): device write + copy back", g, gc->device);
		check(same_copy, what, (double) same_copy, 1.0);

		// auxiliary lattices of the measurements on the device of the slot
		if(meas_aux[0].lattice_aux[0] != NULL)
			{
			Meas_Utils *maux = &(meas_aux[0]);
			for(long r = 0; r < V; r++)
				for(int i = 0; i < STDIM; i++)
					zero(&(maux->lattice_aux[3][r][i]));
			#ifdef OPENACC_MODE
			#pragma acc parallel loop collapse(2) present(gc[0:1], maux[0:1])
			#endif
			for(long r = 0; r < V; r++)
				for(int i = 0; i < STDIM; i++)
					equal(&(maux->lattice_aux[3][r][i]), &(gc->lattice[r][i]));
			#ifdef OPENACC_MODE
			acc_update_self(maux->lattice_aux[3][0], (size_t) V * STDIM * sizeof(GAUGE_GROUP));
			#endif
			int const same_aux = same_lattice(maux->lattice_aux[3], gc->lattice, V);
			sprintf(what, "slot %d (device %d): lattice_aux on device", g, gc->device);
			check(same_aux, what, (double) same_aux, 1.0);
			}

		#ifdef MULTICANONICAL_MODE
		// swap of the cold lattices: lattice_copy_cold = identity on host and device, then swap
		for(long r = 0; r < V; r++)
			for(int i = 0; i < STDIM; i++)
				one(&(gc->lattice_copy_cold[r][i]));
		copy_gauge_conf_to_device(gc, param, GC_COLD);
		swap_lattice_cold(gc);
		// the swap is in the queue of the device of the slot (wait_devices would leave the last device selected)
		#ifdef OPENACC_MODE
		#pragma acc wait(QUEUE_DEVICE)
		#endif
		for(int w = 2; w < 4; w++)
			{
			double const s_dev = lattice_sum_device(gc, w, V);
			double const s_host = lattice_sum_host(lattice_of(gc, w), V);
			sprintf(what, "slot %d (device %d): %s after swap", g, gc->device, lat_name[w]);
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			}

		// views on the device (after the swap): lattice and Z of the view = the fields of the slot
		// (VIEW_FLOW: tested by the gradient flow in test_measures)
		char const *const view_name[3] = {"VIEW_COLD", "VIEW_COPY", "VIEW_COLD_COPY"};
		int const view_lat[3] = {2, 1, 3}; // lattice_cold, lattice_copy, lattice_copy_cold
		for(int k = 0; k <= VIEW_COLD_COPY; k++)
			{
			Gauge_Conf *v = gauge_conf_view(gc, (View_Kind) k);
			#ifdef OPENACC_MODE
			#pragma acc wait(QUEUE_DEVICE)
			#endif
			double s_dev = lattice_sum_device(v, 0, V);
			double s_host = lattice_sum_host(lattice_of(gc, view_lat[k]), V);
			sprintf(what, "slot %d (device %d): %s lattice", g, gc->device, view_name[k]);
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			s_dev = Z_sum_device(v, 0, V, nZ);
			s_host = Z_sum_host((k == VIEW_COLD) ? gc->Z : gc->Z_copy, V, nZ);
			sprintf(what, "slot %d (device %d): %s Z", g, gc->device, view_name[k]);
			check(close_enough(s_dev, s_host), what, s_dev, s_host);
			}
		swap_lattice_cold(gc);
		#ifdef OPENACC_MODE
		#pragma acc wait(QUEUE_DEVICE)
		#endif
		#endif
		}
	}


// sums of the swaps on the device (swap_plaquette_sums, random links): (c_a - c_b) P1 + (c_a^2 - c_b^2) P2 against the
// host sum of swap_action_site over the swap rectangle, for the coefficients of replicas 0, 1 and for an arbitrary pair
static void test_swap_sums(Gauge_Conf *GC, Geometry *geo, GParam *param, Rect_Utils *rect_aux)
	{
	char what[STD_STRING_LENGTH];
	int const mu = param->d_defect_dir;
	Rectangle const *const rect = &(rect_aux->swap_rect);
	int const n_c = (param->d_N_replica_pt > 1) ? 2 : 1;
	double const c_pairs[2][2] = {{0.3, 1.7}, {param->d_pt_bound_cond_coeff[0], param->d_pt_bound_cond_coeff[n_c - 1]}};

	for(int g = 0; g < param->d_N_replica_pt; g++)
		{
		Gauge_Conf *gc = &(GC[g]);
		for(long r = 0; r < param->d_volume; r++)
			for(int i = 0; i < STDIM; i++)
				rand_matrix(&(gc->lattice[r][i]), rng_host());
		copy_gauge_conf_to_device(gc, param, GC_LATTICE);

		double p1, p2;
		swap_plaquette_sums(gc, geo, param, rect, &p1, &p2);

		for(int k = 0; k < n_c; k++)
			{
			double const c_a = c_pairs[k][0];
			double const c_b = c_pairs[k][1];
			double d_host = 0.0;
			for(long n = 0; n < rect->d_vol_rect; n++)
				for(int i = 0; i < STDIM - 1; i++)
					d_host += swap_action_site(gc, geo, param, rect->rect_sites[n], mu, geo->d_orth_dir[mu][i], c_a, c_b);
			double const d_dev = (c_a - c_b) * p1 + (c_a * c_a - c_b * c_b) * p2;
			sprintf(what, "slot %d (device %d): swap D(%.3f, %.3f)", g, gc->device, c_a, c_b);
			check(close_enough(d_dev, d_host), what, d_dev, d_host);
			}
		}
	}


// observables of the configuration: plaquettes, clover energy, charge, Polyakov loops (N_TEST_OBS values)
#define N_TEST_OBS (4 + 2 * STDIM)
static void test_observables(Gauge_Conf *gc, Geometry *geo, GParam *param, double *obs)
	{
	plaquette(gc, geo, param, &(obs[0]), &(obs[1]));
	clover_disc_energy(gc, geo, param, &(obs[2]));
	obs[3] = topcharge(gc, geo, param);
	for(int mu = 0; mu < STDIM; mu++)
		polyakov(gc, geo, param, mu, 1, &(obs[4 + 2 * mu]), &(obs[5 + 2 * mu]));
	}


// operation k of test_measures on gc; out[0], out[1]: scalars to compare besides the links
#define N_TEST_OPS 8
static void test_operation(int const k, Gauge_Conf *gc, Geometry *geo, GParam *param, Meas_Utils *maux, double *out)
	{
	Cooling_Type const types[5] = {LEX_DIR_LEXEO_SITE, LEXEO_SITE_LEX_DIR, RND_DIR_RNDEO_SITE, LEX_SITE_LEX_DIR, LEX_DIR_LEX_SITE};
	out[0] = out[1] = 0.0;
	if(k < 5)
		{
		cooling(gc, geo, param, 2, types[k]);
		}
	else if(k == 5)
		{
		for(int n = 0; n < 3; n++)
			gradflow_RKstep(gc, geo, param, 0.02, maux);
		}
	else if(k == 6)
		{
		double t = 0.0, dt = 0.01;
		int accepted, n_acc = 0;
		for(int n = 0; n < 6; n++)
			{
			gradflow_RKstep_adaptive(gc, geo, param, &t, &dt, &accepted, maux);
			n_acc += accepted;
			}
		out[0] = t;
		out[1] = (double) n_acc;
		}
	else
		{
		cooling(gc, geo, param, 1, LEX_DIR_LEXEO_SITE);
		restore_gauge_conf(gc, param);
		}
	}


// measures on the device against the same functions on the host (force_host_compute, branch if() false of the
// compute regions), for every slot on the same random configuration: observables, the five cooling types, gradient
// flow with fixed and adaptive step (auxiliary lattices of the device, view VIEW_FLOW), restore_gauge_conf
static void test_measures(Gauge_Conf *GC, Geometry *geo, GParam *param, Meas_Utils *meas_aux)
	{
	char what[STD_STRING_LENGTH];
	long const V = param->d_volume;
	int const have_aux = (meas_aux[0].lattice_aux[0] != NULL);
	char const *const obs_name[4] = {"plaquette (spatial)", "plaquette (temporal)", "clover energy", "charge"};
	char const *const op_name[N_TEST_OPS] = {"cooling LEX_DIR_LEXEO_SITE", "cooling LEXEO_SITE_LEX_DIR", "cooling RND_DIR_RNDEO_SITE",
	                                         "cooling LEX_SITE_LEX_DIR", "cooling LEX_DIR_LEX_SITE", "gradient flow",
	                                         "adaptive gradient flow", "restore after cooling"};

	for(int g = 0; g < param->d_N_replica_pt; g++)
		{
		Gauge_Conf *gc = &(GC[g]);

		// reference configuration: random links, lattice_copy = lattice, on the host and on the device
		for(long r = 0; r < V; r++)
			for(int i = 0; i < STDIM; i++)
				{
				rand_matrix(&(gc->lattice[r][i]), rng_host());
				equal(&(gc->lattice_copy[r][i]), &(gc->lattice[r][i]));
				}
		copy_gauge_conf_to_device(gc, param, GC_LATTICE | GC_COPY);

		double obs[2][N_TEST_OBS];
		for(int side = 0; side < 2; side++) // 0: device, 1: host
			{
			int const previous = force_host_compute(side);
			test_observables(gc, geo, param, obs[side]);
			force_host_compute(previous);
			}
		for(int n = 0; n < N_TEST_OBS; n++)
			{
			if(n < 4) sprintf(what, "slot %d (device %d): %s", g, gc->device, obs_name[n]);
			else sprintf(what, "slot %d (device %d): Polyakov %d %s", g, gc->device, (n - 4) / 2, ((n - 4) % 2 == 0) ? "re" : "im");
			check(close_enough(obs[0][n], obs[1][n]), what, obs[0][n], obs[1][n]);
			}

		for(int k = 0; k < N_TEST_OPS; k++)
			{
			if((k == 5 || k == 6) && !have_aux) continue;

			// same random choices (RND_DIR_RNDEO_SITE) on both sides
			rng_t const rng_saved = *rng_host();
			int stdim_saved[STDIM];
			long parity_saved[2][2];
			memcpy(stdim_saved, gc->stdim_shuffle, sizeof(stdim_saved));
			memcpy(parity_saved, gc->parity_shuffle, sizeof(parity_saved));

			double links[2], out[2][2];
			for(int side = 0; side < 2; side++)
				{
				*rng_host() = rng_saved;
				memcpy(gc->stdim_shuffle, stdim_saved, sizeof(stdim_saved));
				memcpy(gc->parity_shuffle, parity_saved, sizeof(parity_saved));
				for(long r = 0; r < V; r++)
					for(int i = 0; i < STDIM; i++)
						equal(&(gc->lattice[r][i]), &(gc->lattice_copy[r][i]));
				copy_gauge_conf_to_device(gc, param, GC_LATTICE);

				int const previous = force_host_compute(side);
				test_operation(k, gc, geo, param, &(meas_aux[0]), out[side]);
				force_host_compute(previous);
				if(side == 0) copy_gauge_conf_to_host(gc, param, GC_LATTICE);
				links[side] = lattice_sum_host(gc->lattice, V);
				}
			sprintf(what, "slot %d (device %d): %s", g, gc->device, op_name[k]);
			check(close_enough(links[0], links[1]), what, links[0], links[1]);
			if(k == 6)
				{
				sprintf(what, "slot %d (device %d): adaptive flow time", g, gc->device);
				check(close_enough(out[0][0], out[1][0]), what, out[0][0], out[1][0]);
				sprintf(what, "slot %d (device %d): adaptive accepted steps", g, gc->device);
				check(out[0][1] == out[1][1], what, out[0][1], out[1][1]);
				}
			}
		}
	}


void real_main(char *in_file)
	{
	Gauge_Conf *GC;
	Geometry geo;
	GParam param;
	Rect_Utils rect_aux;
	Meas_Utils *meas_aux;

	readinput(in_file, &param);
	initrand(param.d_randseed);
	init_geometry(&geo, &param);
	init_rect_utils(&rect_aux, &param);
	init_meas_utils_replica(&meas_aux, &param);
	init_device_utils(&param);
	init_gauge_conf_replica(&GC, &geo, &param);

	int const n_meas = num_meas_utils_replica(&param);
	init_device_data(GC, &geo, &param, &rect_aux, meas_aux, n_meas);
	print_device_utils(stdout, &g_dev);

	test_geometry_data(&geo, &param);
	test_rect_data(&rect_aux, &param);
	test_param_data(GC, &param);
	test_gauge_conf_data(GC, &param, meas_aux);
	test_swap_sums(GC, &geo, &param, &rect_aux); // changes the links
	test_measures(GC, &geo, &param, meas_aux);   // changes the links

	free_device_data(GC, &geo, &param, &rect_aux, meas_aux, n_meas);

	free_gauge_conf_replica(GC, &param);
	free_meas_utils_replica(meas_aux, &param);
	free_rect_utils(&rect_aux, &param);
	free_geometry(&geo, &param);
	free_gparam(&param);

	if(n_fail == 0)
		printf("\nALL %d TESTS PASSED\n", n_test);
	else
		printf("\n%d OF %d TESTS FAILED\n", n_fail, n_test);
	}


int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		printf("Usage: %s input_file\n\n", argv[0]);
		printf("Test of the data on the devices; input file of yang_mills_local_pt_agf.\n");
		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");
	real_main(argv[1]);

	return (n_fail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
	}

#endif
