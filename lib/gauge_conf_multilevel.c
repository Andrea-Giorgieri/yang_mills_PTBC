#ifndef GAUGE_CONF_MULTI_C
#define GAUGE_CONF_MULTI_C

#include "../include/macro.h"
#include "../include/gauge_conf.h"

#include <complex.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/gauge_group.h"
#include "../include/gparam.h"
#include "../include/tens_prod.h"


static inline int get_ml_level(GParam const *const param, int const dt)
	{
	if(dt > param->d_ml_step[0])
		return -1;
	for(int level = 0; level < NLEVELS; level++)
		if(param->d_ml_step[level] == dt)
			return level;
	REQUIRE(0, "failed to determine multilevel level");
	}


void multihit(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param, long const r, int const dir, int const num_hit, GAUGE_GROUP *G)
	{
	if(num_hit > 0)
		{
		GAUGE_GROUP staple, partial;

		zero(G);
		equal(&partial, &GC->lattice[r][dir]);
		#ifndef THETA_MODE
		calcstaples_wilson(GC, geo, param, r, dir, &staple);
		#else
		// compute_clovers in direction "dir" HAS TO BE CALLED BEFORE!
		calcstaples_with_topo(GC, geo, param, r, dir, &staple);
		#endif

		for(int i = 0; i < num_hit; i++)
			{
			single_heatbath(&partial, &staple, param);
			plus_equal(G, &partial);

			unitarize(&partial);
			}
		times_equal_real(G, 1.0 / (double) num_hit);
		}
	else
		{
		equal(G, &GC->lattice[r][dir]);
		}
	}


// compute polyakov loop on a single slice
void compute_local_poly(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param)
	{
	int const num_hit = param->d_ml_num_hit;

	#ifdef THETA_MODE
	// clovers are eventually needed by the multihit
	compute_clovers(GC, geo, param, 0);
	#endif

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		GAUGE_GROUP matrix;

		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);

		one(&GC->loc_poly[slice][r]);
		for(int i = 0; i < param->d_ml_step[NLEVELS - 1]; i++)
			{
			int const t = slice * param->d_ml_step[NLEVELS - 1] + i;
			multihit(GC, geo, param, sisp_and_t_to_si(geo, r, t), 0, num_hit, &matrix);
			times_equal(&GC->loc_poly[slice][r], &matrix);
			}
		}
	}


// perform a complete update on the given level
void update_for_multilevel(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const level)
	{
	#ifdef DEBUG
	ASSERT(param->d_min_size > 1, "this function cannot be used in the completely reduced case");
	#endif

	long const num_even = param->d_n_even;         // number of even sites
	long const even_volume = param->d_even_volume; // volume of largest even sublattice
	long const volume = param->d_volume;           // full lattice volume

	// heatbath
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = 0; s < num_even; s++)
			{
			long const r = geo->d_sip_to_si[s];
			int t;
			long rsp;
			si_to_sisp_and_t(&rsp, &t, geo, r);
			if(t % param->d_ml_step[level] != 0 || dir == 0)
				heatbath(GC, geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = num_even; s < even_volume; s++)
			{
			long const r = geo->d_sip_to_si[s];
			int t;
			long rsp;
			si_to_sisp_and_t(&rsp, &t, geo, r);
			if(t % param->d_ml_step[level] != 0 || dir == 0)
				heatbath(GC, geo, param, r, dir);
			}

		for(long s = even_volume; s < volume; s++)
			{
			long const r = geo->d_sip_to_si[s];
			int t;
			long rsp;
			si_to_sisp_and_t(&rsp, &t, geo, r);
			if(t % param->d_ml_step[level] != 0 || dir == 0)
				heatbath(GC, geo, param, r, dir);
			}
		}

	// overrelax
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		for(int j = 0; j < param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = 0; s < num_even; s++)
				{
				long const r = geo->d_sip_to_si[s];
				int t;
				long rsp;
				si_to_sisp_and_t(&rsp, &t, geo, r);
				if(t % param->d_ml_step[level] != 0 || dir == 0)
					overrelaxation(GC, geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = num_even; s < even_volume; s++)
				{
				long const r = geo->d_sip_to_si[s];
				int t;
				long rsp;
				si_to_sisp_and_t(&rsp, &t, geo, r);
				if(t % param->d_ml_step[level] != 0 || dir == 0)
					overrelaxation(GC, geo, param, r, dir);
				}

			for(long s = even_volume; s < volume; s++)
				{
				long const r = geo->d_sip_to_si[s];
				int t;
				long rsp;
				si_to_sisp_and_t(&rsp, &t, geo, r);
				if(t % param->d_ml_step[level] != 0 || dir == 0)
					overrelaxation(GC, geo, param, r, dir);
				}
			}
		}

	accept_gauge_conf(GC, param);
	}


// initialize ml_polycorr[level] to 0
static inline void zero_polycorr_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		zero_TensProd(&GC->ml_polycorr[level][slice][r]);
		}
	}


// initialize ml_polycorr[level] and ml_polyplaq[level] to 0
static inline void zero_polycorr_polyplaq_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		zero_TensProd(&GC->ml_polycorr[level][slice][r]);
		if(slice == 0) zero_TensProd(&GC->ml_polyplaq[level][r]);
		}
	}


// initialize ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level] to 0
static inline void zero_polycorr_polyplaq_polyplaqconn_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		zero_TensProd(&GC->ml_polycorr[level][slice][r]);
		if(slice == 0)
			{
			zero_TensProd(&GC->ml_polyplaq[level][r]);
			zero_TensProd(&GC->ml_polyplaqconn[level][r]);
			}
		}
	}


// normalize ml_polycorr[level]
static inline void normalize_polycorr_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		times_equal_real_TensProd(&GC->ml_polycorr[level][slice][r], 1.0 / (double) param->d_ml_upd[level]);
		}
	}


// normalize ml_polycorr[level] for long simulations
static inline void normalize_polycorr_long_level0(Gauge_Conf *const GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[0]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[0];
		int const slice = (int) (raux % param->d_ml_num_slices[0]);
		times_equal_real_TensProd(&GC->ml_polycorr[0][slice][r], 1.0 / ((double) param->d_ml_upd[0] * param->d_ml_level0_repeat));
		}
	}


// normalize ml_polycorr[level] and ml_polyplaq[level]
static inline void normalize_polycorr_polyplaq_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		times_equal_real_TensProd(&GC->ml_polycorr[level][slice][r], 1.0 / (double) param->d_ml_upd[level]);
		if(slice == 0) times_equal_real_TensProd(&GC->ml_polyplaq[level][r], 1.0 / (double) param->d_ml_upd[level]);
		}
	}


// normalize ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
static inline void normalize_polycorr_polyplaq_polyplaqconn_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		times_equal_real_TensProd(&GC->ml_polycorr[level][slice][r], 1.0 / (double) param->d_ml_upd[level]);
		if(slice == 0)
			{
			times_equal_real_TensProd(&GC->ml_polyplaq[level][r], 1.0 / (double) param->d_ml_upd[level]);
			times_equal_real_TensProd(&GC->ml_polyplaqconn[level][r], 1.0 / (double) param->d_ml_upd[level]);
			}
		}
	}


// normalize polycorr[level], polyplaq[level] and polyplaqconn[level] for long simulations
static inline void normalize_polycorr_polyplaq_polyplaqconn_long_level0(Gauge_Conf *const GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[0]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[0];
		int const slice = (int) (raux % param->d_ml_num_slices[0]);

		times_equal_real_TensProd(&GC->ml_polycorr[0][slice][r], 1.0 / ((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));

		if(slice == 0)
			{
			times_equal_real_TensProd(&GC->ml_polyplaq[0][r], 1.0 / ((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));
			times_equal_real_TensProd(&GC->ml_polyplaqconn[0][r], 1.0 / ((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));
			}
		}
	}


static inline void accumulate_polycorr_intermediate_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		one_TensProd(&TP);
		for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
			{
			times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
			}
		plus_equal_TensProd(&GC->ml_polycorr[level][slice][r], &TP);
		}
	}


static inline void accumulate_polycorr_polyplaq_intermediate_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		one_TensProd(&TP);
		for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
			{
			times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
			}
		plus_equal_TensProd(&GC->ml_polycorr[level][slice][r], &TP);
		if(slice == 0)
			{
			equal_TensProd(&TP, &GC->ml_polyplaq[level + 1][r]);
			for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
				{
				times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
				}
			plus_equal_TensProd(&GC->ml_polyplaq[level][r], &TP);
			}
		}
	}


static inline void accumulate_polycorr_polyplaq_polyplaqconn_intermediate_level(Gauge_Conf *const GC, GParam const *const param, int const level)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[level]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[level];
		int const slice = (int) (raux % param->d_ml_num_slices[level]);
		one_TensProd(&TP);
		for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
			{
			times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
			}
		plus_equal_TensProd(&GC->ml_polycorr[level][slice][r], &TP);
		if(slice == 0)
			{
			equal_TensProd(&TP, &GC->ml_polyplaq[level + 1][r]);
			for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
				{
				times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
				}
			plus_equal_TensProd(&GC->ml_polyplaq[level][r], &TP);
			equal_TensProd(&TP, &GC->ml_polyplaqconn[level + 1][r]);
			for(int j = 0; j < param->d_ml_step[level] / param->d_ml_step[level + 1]; j++)
				{
				times_equal_TensProd(&TP, &GC->ml_polycorr[level + 1][slice * param->d_ml_step[level] / param->d_ml_step[level + 1] + j][r]);
				}
			plus_equal_TensProd(&GC->ml_polyplaqconn[level][r], &TP);
			}
		}
	}


static inline void accumulate_polycorr_innermost_level(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);
		long r1 = sisp_and_t_to_si(geo, r, 0); // r is a 3d index, r1 is the 4d index value of (r,t=0)
		for(int j = 0; j < param->d_dist_poly; j++) r1 = nnp(geo, r1, 1);
		long r2;
		int t_tmp;
		si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
		TensProd_init(&TP, &GC->loc_poly[slice][r], &GC->loc_poly[slice][r2]);
		plus_equal_TensProd(&GC->ml_polycorr[NLEVELS - 1][slice][r], &TP);
		}
	}


static inline void accumulate_polycorr_polyplaq_innermost_level(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);
		long r1 = sisp_and_t_to_si(geo, r, 0);
		for(int j = 0; j < param->d_dist_poly; j++) r1 = nnp(geo, r1, 1);
		long r2;
		int t_tmp;
		si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
		TensProd_init(&TP, &GC->loc_poly[slice][r], &GC->loc_poly[slice][r2]);
		plus_equal_TensProd(&GC->ml_polycorr[NLEVELS - 1][slice][r], &TP);
		if(slice == 0)
			{
			times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
			plus_equal_TensProd(&GC->ml_polyplaq[NLEVELS - 1][r], &TP);
			}
		}
	}


static inline void accumulate_polycorr_polyplaq_polyplaqconn_innermost_level(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		TensProd TP;
		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);
		long r1 = sisp_and_t_to_si(geo, r, 0);
		for(int j = 0; j < param->d_dist_poly; j++) r1 = nnp(geo, r1, 1);
		long r2;
		int t_tmp;
		si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
		TensProd_init(&TP, &GC->loc_poly[slice][r], &GC->loc_poly[slice][r2]);
		plus_equal_TensProd(&GC->ml_polycorr[NLEVELS - 1][slice][r], &TP);
		if(slice == 0)
			{
			times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
			plus_equal_TensProd(&GC->ml_polyplaq[NLEVELS - 1][r], &TP);
			TensProd_init(&TP, &GC->loc_plaqconn[r], &GC->loc_poly[slice][r2]);
			plus_equal_TensProd(&GC->ml_polyplaqconn[NLEVELS - 1][r], &TP);
			}
		}
	}


// multilevel for polyakov correlator
void multilevel_polycorr(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const dt)
	{
	// d_size[0] >= d_ml_step[0] > d_ml_step[1] > ...

	// determine the level to be used
	int const level = get_ml_level(param, dt);

	if(level == -1)
		{
		multilevel_polycorr(GC, geo, param, param->d_ml_step[0]);
		return;
		}

	REQUIRE(level >= 0 && level < NLEVELS, "error in the multilevel");

	// initialize ml_polycorr[level] to zero
	zero_polycorr_level(GC, param, level);

	// perform the update
	for(int upd = 0; upd < param->d_ml_upd[level]; upd++)
		{
		update_for_multilevel(GC, geo, param, level);

		if(level == NLEVELS - 1)
			{
			// compute Polyakov loop restricted to the slice
			compute_local_poly(GC, geo, param);

			// compute the tensor products and update ml_polycorr[level]
			accumulate_polycorr_innermost_level(GC, geo, param);
			}
		else
			{
			// recursive call to next level
			multilevel_polycorr(GC, geo, param, param->d_ml_step[level + 1]);

			// update polycorr[level] with polycorr[level+1]
			accumulate_polycorr_intermediate_level(GC, param, level);
			}
		}

	// normalize ml_polycorr[level]
	normalize_polycorr_level(GC, param, level);
	}


// multilevel for polyakov correlator to be used in long simulations
void multilevel_polycorr_long_zero(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const iteration)
	{
	// initialize ml_polycorr[0] to 0 if needed
	if(iteration == 0)
		{
		zero_polycorr_level(GC, param, 0);
		}

	// perform the update
	for(int upd = 0; upd < param->d_ml_upd[0]; upd++)
		{
		// update on level zero
		update_for_multilevel(GC, geo, param, 0);

		#if NLEVELS == 1
		// compute Polyakov loop restricted to the slice
		compute_local_poly(GC, geo, param);

		// compute the tensor products and update ml_polycorr[0]
		accumulate_polycorr_innermost_level(GC, geo, param);

		#else
		// important: call inner levels with the "non-long" version
		multilevel_polycorr(GC, geo, param, param->d_ml_step[1]);

		// update polycorr[0] with polycorr[1]
		accumulate_polycorr_intermediate_level(GC, param, 0);
		#endif
		}

	if(iteration == param->d_ml_level0_repeat - 1) // iteration starts from zero
		{
		// normalize polycorr[0]
		normalize_polycorr_long_level0(GC, param);
		}
	}


// compute polyakov loop and plaquette on a single slice
void compute_local_poly_and_plaq(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param)
	{
	int const num_hit = param->d_ml_num_hit;

	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0); // clovers are eventually needed by the multihit
	#endif

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);

		GAUGE_GROUP matrix;
		one(&GC->loc_poly[slice][r]);
		for(int i = 0; i < param->d_ml_step[NLEVELS - 1]; i++)
			{
			int const t = slice * param->d_ml_step[NLEVELS - 1] + i;
			multihit(GC, geo, param, sisp_and_t_to_si(geo, r, t), 0, num_hit, &matrix);
			times_equal(&GC->loc_poly[slice][r], &matrix);
			}

		if(slice == 0)
			{
			// moves to the correct position of the plaquette:
			// polyakov loop are separated along direction 1 and the transverse direction is 2
			long r4 = sisp_and_t_to_si(geo, r, 1); // t = 1
			for(int j = 0; j < param->d_dist_poly / 2; j++) r4 = nnp(geo, r4, 1);
			for(int j = 0; j < param->d_trasv_dist; j++) r4 = nnp(geo, r4, 2);

			GC->loc_plaq[r] = plaquettep_complex(GC, geo, param, r4, param->d_plaq_dir[0], param->d_plaq_dir[1]);
			}
		}
	}


// multilevel for flux width computation using the disconnected correlator
void multilevel_tube_disc(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const dt)
	{
	// d_size[0] >= d_ml_step[0] > d_ml_step[1] > ...

	// determine the level to be used
	const int level = get_ml_level(param, dt);

	if(level == -1)
		{
		multilevel_tube_disc(GC, geo, param, param->d_ml_step[0]);
		return;
		}

	REQUIRE(level >= 0 && level < NLEVELS, "invalid multilevel level");

	// initialize ml_polycorr[level] and ml_polyplaq[level]
	zero_polycorr_polyplaq_level(GC, param, level);

	// perform the update
	for(int upd = 0; upd < param->d_ml_upd[level]; ++upd)
		{
		update_for_multilevel(GC, geo, param, level);

		if(level == NLEVELS - 1)
			{
			// compute the Polyakov loop and plaquette restricted to the slice
			compute_local_poly_and_plaq(GC, geo, param);

			// compute the tensor products and update ml_polycorr[level] and ml_polyplaq[level]
			accumulate_polycorr_polyplaq_innermost_level(GC, geo, param);
			}
		else
			{
			// recursive call to next level
			multilevel_tube_disc(GC, geo, param, param->d_ml_step[level + 1]);

			// update ml_polycorr[level] and ml_polyplaq[level] with ml_polycorr[level + 1] and ml_polyplaq[level + 1]
			accumulate_polycorr_polyplaq_intermediate_level(GC, param, level);
			}
		}

	// normalize ml_polycorr[level] and ml_polyplaq[level]
	normalize_polycorr_polyplaq_level(GC, param, level);
	}


// compute polyakov loop, plaquette and connected plaquette on a single slice
void compute_local_poly_plaq_and_plaqconn(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param)
	{
	int const num_hit = param->d_ml_num_hit;

	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0); // clovers are eventually needed by the multihit
	#endif

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long raux = 0; raux < param->d_space_vol[0] * param->d_ml_num_slices[NLEVELS - 1]; raux++)
		{
		long const r = raux / param->d_ml_num_slices[NLEVELS - 1];
		int const slice = (int) (raux % param->d_ml_num_slices[NLEVELS - 1]);

		GAUGE_GROUP matrix;
		one(&GC->loc_poly[slice][r]);
		for(int i = 0; i < param->d_ml_step[NLEVELS - 1]; i++)
			{
			int const t = slice * param->d_ml_step[NLEVELS - 1] + i;
			multihit(GC, geo, param, sisp_and_t_to_si(geo, r, t), 0, num_hit, &matrix);
			times_equal(&GC->loc_poly[slice][r], &matrix);
			}

		if(slice == 0)
			{
			long r4 = sisp_and_t_to_si(geo, r, 0); // t=0 starting point

			multihit(GC, geo, param, r4, 0, num_hit, &matrix);
			equal(&GC->loc_plaqconn[r], &matrix);

			r4 = nnp(geo, r4, 0); // now we are in the t=1 plane

			#ifdef DEBUG
			long const r4start = r4;
			#endif

			// moves to the correct position of the plaquette:
			// polyakov loop are separated along direction 1 and the transverse direction is 2
			for(int j = 0; j < param->d_dist_poly / 2; j++)
				{
				times_equal(&GC->loc_plaqconn[r], &GC->lattice[r4][1]);
				r4 = nnp(geo, r4, 1);
				}
			for(int j = 0; j < param->d_trasv_dist; j++)
				{
				times_equal(&GC->loc_plaqconn[r], &GC->lattice[r4][2]);
				r4 = nnp(geo, r4, 2);
				}
			plaquettep_matrix(GC, geo, param, r4, param->d_plaq_dir[0], param->d_plaq_dir[1], &matrix);

			GC->loc_plaq[r] = retr(&matrix) + I * imtr(&matrix);
			times_equal(&GC->loc_plaqconn[r], &matrix);

			for(int j = 0; j < param->d_trasv_dist; j++)
				{
				r4 = nnm(geo, r4, 2);
				times_equal_dag(&GC->loc_plaqconn[r], &GC->lattice[r4][2]);
				}
			for(int j = 0; j < param->d_dist_poly / 2; j++)
				{
				r4 = nnm(geo, r4, 1);
				times_equal_dag(&GC->loc_plaqconn[r], &GC->lattice[r4][1]);
				}

			#ifdef DEBUG
			ASSERT(r4start == r4, "Polyakov loop is not closed: start = %ld, end = %ld", r4start, r4);
			#endif

			for(int j = 1; j < param->d_ml_step[NLEVELS - 1]; j++)
				{
				multihit(GC, geo, param, r4, 0, num_hit, &matrix);
				times_equal(&GC->loc_plaqconn[r], &matrix);
				r4 = nnp(geo, r4, 0);
				}
			}
		}
	}


// multilevel for flux width computation using the connected correlator
void multilevel_tube_conn(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const dt)
	{
	// d_size[0] >= d_ml_step[0] > d_ml_step[1] > ...

	// determine the level to be used
	const int level = get_ml_level(param, dt);

	if(level == -1)
		{
		multilevel_tube_conn(GC, geo, param, param->d_ml_step[0]);
		return;
		}

	REQUIRE(level >= 0 && level < NLEVELS, "invalid multilevel level");

	// initialize ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
	zero_polycorr_polyplaq_polyplaqconn_level(GC, param, level);

	// perform the update
	for(int upd = 0; upd < param->d_ml_upd[level]; ++upd)
		{
		update_for_multilevel(GC, geo, param, level);

		if(level == NLEVELS - 1)
			{
			// compute the Polyakov loop, plaquette and connected plaquette restricted to the slice
			compute_local_poly_plaq_and_plaqconn(GC, geo, param);

			// compute the tensor products and update ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
			accumulate_polycorr_polyplaq_polyplaqconn_innermost_level(GC, geo, param);
			}
		else
			{
			// recursive call to next level
			multilevel_tube_conn(GC, geo, param, param->d_ml_step[level + 1]);

			// update ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
			// with ml_polycorr[level + 1], ml_polyplaq[level + 1] and ml_polyplaqconn[level + 1]
			accumulate_polycorr_polyplaq_polyplaqconn_intermediate_level(GC, param, level);
			}
		}

	// normalize ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
	normalize_polycorr_polyplaq_polyplaqconn_level(GC, param, level);
	}


// multilevel for flux tube with connected correlator to be used in long simulations
void multilevel_tube_conn_long_zero(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const iteration)
	{
	// initialize ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] to 0 if needed
	if(iteration == 0)
		{
		zero_polycorr_polyplaq_polyplaqconn_level(GC, param, 0);
		}

	// perform the update
	for(int upd = 0; upd < param->d_ml_upd[0]; upd++)
		{
		// update on level zero
		update_for_multilevel(GC, geo, param, 0);

		#if NLEVELS == 1
		// compute Polyakov loop, plaquette and connected plaquette restricted to the slice
		compute_local_poly_plaq_and_plaqconn(GC, geo, param);

		// compute the tensor products and update ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0]
		accumulate_polycorr_polyplaq_polyplaqconn_innermost_level(GC, geo, param);
		#else
		// important: call inner levels with the "non-long" version
		multilevel_tube_conn(GC, geo, param, param->d_ml_step[1]);

		// update polycorr[0] with polycorr[1]
		accumulate_polycorr_polyplaq_polyplaqconn_intermediate_level(GC, param, 0);
		#endif
		}

	if(iteration == param->d_ml_level0_repeat - 1) // iteration starts from zero
		{
		// normalize polycorr[0], polyplaq[0] and polyplaqconn[0]
		normalize_polycorr_polyplaq_polyplaqconn_long_level0(GC, param);
		}
	}


void perform_multilevel_long_update_zero(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const iteration,
                                         Multilevel_Obs const ml_obs)
	{
	switch(ml_obs)
		{
		case NONE:
			break;
		case POLYCORR_LONG:
			multilevel_polycorr_long_zero(GC, geo, param, iteration);
			break;
		case TUBE_CONN_LONG:
			multilevel_tube_conn_long_zero(GC, geo, param, iteration);
			break;
		default:
			REQUIRE(0, "unknown multilevel observable (%d)\n", (int)ml_obs);
		}
	}

#endif
