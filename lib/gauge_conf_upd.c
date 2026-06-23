#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
#include<omp.h>
#endif
#include<stdlib.h>

#include"../include/memalign.h"
#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// compute the staple in position r, direction i and save it in M
void calcstaples_wilson(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        long r,
                        int i,
                        GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#else
	(void) param;
	#endif

	double complex factor;
	GAUGE_GROUP link1, link2, link3, link12, stap;

	zero(M);

	long const rpi = nnp(geo, r, i);

	for(int l = i + 1; l < i + STDIM; l++)
		{
		int j = l % STDIM;

		long const rpj = nnp(geo, r, j);
		long const rmj = nnm(geo, r, j);
		long const rpimj = nnm(geo, rpi, j);

		//
		//    i ^
		//      |    (1)
		//      +----->-----+
		//      |           |
		//                  |
		//      |           V (2)
		//                  |
		//      |           |
		//      +-----<-----+--> j
		//      r    (3)
		//

		// staple
		equal(&link1, &(GC->lattice[rpi][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rpj][i])); // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));   // link3 = (3)

		times_dag2(&link12, &link1, &link2); // link12 = link1 * link2^{dag}
		times_dag2(&stap, &link12, &link3);  // stap = link12 * link3^{dag}

		// twist (clockwise plaquette) modification
		factor = GC->Z[r][dirs_to_si(i, j)]; // Z_\mu\nu(x)
		times_equal_complex(&stap, factor);  // Z_\mu\nu(x) * staple

		// accumulate M
		plus_equal(M, &stap);

		//
		//    i ^
		//      |    (1)
		//      |-----<-----+
		//      |           |
		//      |
		//  (2) V           |
		//      |
		//      |           |
		//      +----->-----+--> j
		//           (3)    r
		//

		// staple
		equal(&link1, &(GC->lattice[rpimj][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rmj][i]));   // link2 = (2)
		equal(&link3, &(GC->lattice[rmj][j]));   // link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12 = link1^{dag} * link2^{dag}
		times(&stap, &link12, &link3);        // stap = link12 * link3

		// twist (anticlockwise plaquette) modification
		factor = GC->Z[rmj][dirs_to_si(j, i)]; // Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		times_equal_complex(&stap, factor);    // Z_\mu\nu(x-\nu) * staple

		// accumulate M
		plus_equal(M, &stap);
		}
	}

// compute the staple for the trace deformed theory:
// in practice a Polyakov loop without a link
void calcstaples_tracedef(Gauge_Conf const *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          long r,
                          int i,
                          GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(i == 0, "trace-deformed staple of a non-temporal link (mu = %d)", i);
	#endif

	if(i != 0)
		{
		zero(M);
		}
	else
		{
		GAUGE_GROUP aux;

		one(&aux);

		long int r_next = r;
		for(int j = 1; j < param->d_size[0]; j++)
			{
			r_next = nnp(geo, r_next, 0);
			times_equal(&aux, &(GC->lattice[r_next][0]));
			}
		equal(M, &aux);
		}
	}


// compute all the clovers in directions orthogonal to "dir"
void compute_clovers(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     int dir)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		GAUGE_GROUP aux;

		for(int i = 0; i < STDIM; i++)
			{
			for(int j = i + 1; j < STDIM; j++)
				{
				if(i != dir && j != dir)
					{
					clover(GC, geo, param, r, i, j, &aux);

					equal(&(GC->clover_array[r][i][j]), &aux);
					minus_equal_dag(&(GC->clover_array[r][i][j]), &aux); // clover_array[r][i][j] = aux - aux^{dag}

					equal(&(GC->clover_array[r][j][i]), &(GC->clover_array[r][i][j]));
					times_equal_real(&(GC->clover_array[r][j][i]), -1.0); // clover_array[r][j][i] = -clover_array[r][i][j]
					}
				}
			}
		}
	}


// compute the staple in position r, direction i and save it in M
// when an imaginary theta term is present
void calcstaples_with_topo(Gauge_Conf const *const GC,
                           Geometry const *const geo,
                           GParam const *const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, M);
	#else
	// the theta term is written as
	// theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood[\mu][\nu][0] sood[\mu][\nu][1]} )
	// the independent permutations are 0123 0231 0312
	// Q_{\mu\nu} is the clover on the plane (\mu,\nu)
	// GC->clover_array[r][mu][nu] = (Q - Q^{dag})_{\mu\nu}(r)

	GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, topo_M, aux;
	double complex factor;

	zero(M);
	zero(&topo_M);

	long const rpi = nnp(geo, r, i);
	double const Zr = creal(GC->Z[r][param->d_n_planes]);
	double const Zrpi = creal(GC->Z[rpi][param->d_n_planes]);

	for(int l = i + 1; l < i + STDIM; l++)
		{
		int const j = l % STDIM;
		int const sood0 = g_signed_ord_orth_dir[i][j][0];
		int const sood1 = g_signed_ord_orth_dir[i][j][1];

		long const rpj = nnp(geo, r, j);
		long const rmj = nnm(geo, r, j);
		long const rpipj = nnp(geo, rpi, j);
		long const rpimj = nnm(geo, rpi, j);

		double const Zrpj = creal(GC->Z[rpj][param->d_n_planes]);
		double const Zrmj = creal(GC->Z[rmj][param->d_n_planes]);
		double const Zrpipj = creal(GC->Z[rpipj][param->d_n_planes]);
		double const Zrpimj = creal(GC->Z[rpimj][param->d_n_planes]);

		GAUGE_GROUP const *CAr = &(GC->clover_array[r][sood0][sood1]);
		GAUGE_GROUP const *CArpi = &(GC->clover_array[rpi][sood0][sood1]);
		GAUGE_GROUP const *CArpj = &(GC->clover_array[rpj][sood0][sood1]);
		GAUGE_GROUP const *CArmj = &(GC->clover_array[rmj][sood0][sood1]);
		GAUGE_GROUP const *CArpipj = &(GC->clover_array[rpipj][sood0][sood1]);
		GAUGE_GROUP const *CArpimj = &(GC->clover_array[rpimj][sood0][sood1]);

		//
		//    i ^
		//      |    (1)
		//  (b) +----->-----+ (c)
		//      |           |
		//                  |
		//      |           V (2)
		//                  |
		//      |           |
		//  (a) +-----<-----+--> j
		//      r    (3)   (d)
		//

		// non-topo staple
		equal(&link1, &(GC->lattice[rpi][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rpj][i])); // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));   // link3 = (3)

		times_dag2(&link12, &link1, &link2); // link12 = link1 * link2^{dag}
		times_dag2(&stap, &link12, &link3);  // stap = link12 * link3^{dag}

		// clover insertion in (a)
		times(&aux, &stap, CAr);    // stap * clover
		times_equal_real(&aux, Zr); // bulk factor
		equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, CArpi, &stap);    // clover * stap
		times_equal_real(&aux, Zrpi); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link1, CArpipj);   // link1 * clover
		times_equal_dag(&aux, &link2);  // *= link2^{dag}
		times_equal_dag(&aux, &link3);  // *= link3^{dag}
		times_equal_real(&aux, Zrpipj); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times(&aux, &link12, CArpj);   // link1 * link2 * clover
		times_equal_dag(&aux, &link3); // *= link3^{dag}
		times_equal_real(&aux, Zrpj);  // bulk factor
		plus_equal(&topo_stap, &aux);

		// twist modification (clockwise plaquette)
		factor = GC->Z[r][dirs_to_si(i, j)];     // Z_\mu\nu(x)
		times_equal_complex(&topo_stap, factor); // Z_\mu\nu(x) * topo_staple
		times_equal_complex(&stap, factor);      // Z_\mu\nu(x) * staple

		// accumulate M and topo_M
		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);

		//
		//   i ^
		//     |   (1)
		// (d) +----<------+ (a)
		//     |           |
		//     |
		// (2) V           |
		//     |
		//     |           | (b)
		// (c) +------>----+---> j
		//           (3)   r
		//

		// non-topo staple
		equal(&link1, &(GC->lattice[rpimj][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rmj][i]));   // link2 = (2)
		equal(&link3, &(GC->lattice[rmj][j]));   // link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12 = link1^{dag} * link2^{dag}
		times(&stap, &link12, &link3);        // stap = link12 * link3

		// clover insertion in (a)
		times(&aux, CArpi, &stap);    // clover * stap
		times_equal_real(&aux, Zrpi); // bulk factor
		equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &stap, CAr);    // stap * clover
		times_equal_real(&aux, Zr); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link12, CArmj);  // link1^{dag} * link2^{dag} * clover
		times_equal(&aux, &link3);    // *= link3
		times_equal_real(&aux, Zrmj); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times_dag1(&aux, &link1, CArpimj); // link1^{dag} * clover
		times_equal_dag(&aux, &link2);     // *= link2^{dag}
		times_equal(&aux, &link3);         // *= link3
		times_equal_real(&aux, Zrpimj);    // bulk factor
		plus_equal(&topo_stap, &aux);

		// twist modification (anticlockwise plaquette)
		factor = GC->Z[rmj][dirs_to_si(j, i)];   // Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		times_equal_complex(&topo_stap, factor); // Z_\nu\mu(x-\nu) * topo_staple
		times_equal_complex(&stap, factor);      // Z_\nu\mu(x-\nu) * staple

		// accumulate M and topo_M (with minus sign from definition of topological charge)
		plus_equal(M, &stap);
		minus_equal(&topo_M, &topo_stap);
		}

	// multiply topo_staple by theta coefficient and sum to staple
	double const topo_staple_coeff = (param->d_theta) * ((double) NCOLOR) / (param->d_beta * 128.0 * PI * PI);
	times_equal_real(&topo_M, topo_staple_coeff);
	plus_equal(M, &topo_M);

	#endif
	}

// update functions for parallel tempering on defect

// compute all the clovers in directions orthogonal to "dir" for all replica
void compute_clovers_replica(Gauge_Conf const *const GC,
                             Geometry const *const geo,
                             GParam const *const param,
                             int dir)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < ((param->d_N_replica_pt) * param->d_volume); s++)
		{
		GAUGE_GROUP aux;
		// s = i_r * volume + r
		long const r = s % (param->d_volume);                // space-time index
		int const i_r = (int) ((s - r) / (param->d_volume)); // replica index

		for(int i = 0; i < STDIM; i++)
			{
			for(int j = i + 1; j < STDIM; j++)
				{
				if(i != dir && j != dir)
					{
					clover(&(GC[i_r]), geo, param, r, i, j, &aux);

					equal(&(GC[i_r].clover_array[r][i][j]), &aux);
					minus_equal_dag(&(GC[i_r].clover_array[r][i][j]), &aux); // clover_array[r][i][j] = aux - aux^{dag}

					equal(&(GC[i_r].clover_array[r][j][i]), &(GC[i_r].clover_array[r][i][j]));
					times_equal_real(&(GC[i_r].clover_array[r][j][i]), -1.0); // clover_array[r][j][i] = -clover_array[r][i][j]
					}
				}
			}
		}
	}

// compute all the clovers in directions orthogonal to "dir" for all replica on a given rectangle
void compute_clovers_replica_rect(Gauge_Conf const *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  int dir,
                                  Rectangle const *const clover_rectangle)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < ((param->d_N_replica_pt) * (clover_rectangle->d_vol_rect)); s++)
		{
		GAUGE_GROUP aux;

		// s = i_r * volume_rect + n
		long const n = s % (clover_rectangle->d_vol_rect);                // space-time index on the rectangle
		long const r = clover_rectangle->rect_sites[n];                   // space-time index on the lattice
		int const i_r = (int) ((s - n) / (clover_rectangle->d_vol_rect)); // replica index

		for(int i = 0; i < STDIM; i++)
			{
			for(int j = i + 1; j < STDIM; j++)
				{
				if(i != dir && j != dir)
					{
					clover(&(GC[i_r]), geo, param, r, i, j, &aux);

					equal(&(GC[i_r].clover_array[r][i][j]), &aux);
					minus_equal_dag(&(GC[i_r].clover_array[r][i][j]), &aux); // clover_array[r][i][j] = aux - aux^{dag}

					equal(&(GC[i_r].clover_array[r][j][i]), &(GC[i_r].clover_array[r][i][j]));
					times_equal_real(&(GC[i_r].clover_array[r][j][i]), -1.0); // clover_array[r][j][i] = -clover_array[r][i][j]
					}
				}
			}
		}
	}

// evaluate non-topo staples with defect and twist factors in position r and direction i and save it in M
void calcstaples_wilson_with_defect(Gauge_Conf const *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    long r,
                                    int i,
                                    GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#else
	(void) param;
	#endif

	GAUGE_GROUP link1, link2, link3, link12, stap;
	double complex factor;

	zero(M);

	long const rpi = nnp(geo, r, i);

	for(int l = i + 1; l < i + STDIM; l++)
		{
		int const j = l % STDIM;

		long const rpj = nnp(geo, r, j);
		long const rmj = nnm(geo, r, j);
		long const rpimj = nnm(geo, rpi, j);

		//
		//    i ^
		//      |    (1)
		//      +----->-----+
		//      |           |
		//                  |
		//      |           V (2)
		//                  |
		//      |           |
		//      +-----<-----+--> j
		//      r    (3)
		//

		// staple
		equal(&link1, &(GC->lattice[rpi][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rpj][i])); // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));   // link3 = (3)

		times_dag2(&link12, &link1, &link2); // link12 = link1 * link2^{dag}
		times_dag2(&stap, &link12, &link3);  // stap = link12 * link3^{dag}

		// boundary condition and twist (clockwise plaquette) modification
		factor = GC->Z[r][dirs_to_si(i, j)];                                         // Z_\mu\nu(x)
		factor *= (GC->C[r][i]) * (GC->C[rpi][j]) * (GC->C[rpj][i]) * (GC->C[r][j]); // *= K_\mu\nu(x)
		times_equal_complex(&stap, factor);                                          // K_\mu\nu(x) * Z_\mu\nu(x) * staple

		// accumulate M
		plus_equal(M, &stap);

		//
		//    i ^
		//      |    (1)
		//      |-----<-----+
		//      |           |
		//      |
		//  (2) V           |
		//      |
		//      |           |
		//      +----->-----+--> j
		//           (3)    r
		//

		// staple
		equal(&link1, &(GC->lattice[rpimj][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rmj][i]));   // link2 = (2)
		equal(&link3, &(GC->lattice[rmj][j]));   // link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12 = link1^{dag} * link2^{dag}
		times(&stap, &link12, &link3);        // stap = link12 * link3

		// boundary condition and twist (anticlockwise plaquette) modification
		factor = GC->Z[rmj][dirs_to_si(j, i)];                                           // Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		factor *= (GC->C[rmj][i]) * (GC->C[rpimj][j]) * (GC->C[r][i]) * (GC->C[rmj][j]); // *= K_\mu\nu(x-\nu)
		times_equal_complex(&stap, factor);                                              // K_\mu\nu(x-\nu) * Z_\nu\mu(x-\nu) * staple

		// accumulate M
		plus_equal(M, &stap);
		}
	}


// evaluate topo staples with defect on the non-topo contributions
// in position r and direction i and save it in M
void calcstaples_with_topo_with_defect(Gauge_Conf const *const GC,
                                       Geometry const *const geo,
                                       GParam const *const param,
                                       long r,
                                       int i,
                                       GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, M);
	#else
	// the theta term is written as
	// theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood[\mu][\nu][0] sood[\mu][\nu][1]} )
	// the independent permutations are 0123 0231 0312
	// Q_{\mu\nu} is the clover on the plane (\mu,\nu)
	// GC->clover_array[r][mu][nu] = (Q - Q^{dag})_{\mu\nu}(r)

	GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, topo_M, aux;

	double const topo_staple_coeff = (param->d_theta) * ((double) NCOLOR) / (param->d_beta * 128.0 * PI * PI);
	double complex factor;

	zero(M);
	zero(&topo_M);

	long const rpi = nnp(geo, r, i);
	double const Zr = creal(GC->Z[r][param->d_n_planes]);
	double const Zrpi = creal(GC->Z[rpi][param->d_n_planes]);

	for(int l = i + 1; l < i + STDIM; l++)
		{
		int j = l % STDIM;
		int const sood0 = g_signed_ord_orth_dir[i][j][0];
		int const sood1 = g_signed_ord_orth_dir[i][j][1];

		long const rpj = nnp(geo, r, j);
		long const rmj = nnm(geo, r, j);
		long const rpipj = nnp(geo, rpi, j);
		long const rpimj = nnm(geo, rpi, j);

		double const Zrpj = creal(GC->Z[rpj][param->d_n_planes]);
		double const Zrmj = creal(GC->Z[rmj][param->d_n_planes]);
		double const Zrpipj = creal(GC->Z[rpipj][param->d_n_planes]);
		double const Zrpimj = creal(GC->Z[rpimj][param->d_n_planes]);

		GAUGE_GROUP const *CAr = &(GC->clover_array[r][sood0][sood1]);
		GAUGE_GROUP const *CArpi = &(GC->clover_array[rpi][sood0][sood1]);
		GAUGE_GROUP const *CArpj = &(GC->clover_array[rpj][sood0][sood1]);
		GAUGE_GROUP const *CArmj = &(GC->clover_array[rmj][sood0][sood1]);
		GAUGE_GROUP const *CArpipj = &(GC->clover_array[rpipj][sood0][sood1]);
		GAUGE_GROUP const *CArpimj = &(GC->clover_array[rpimj][sood0][sood1]);

		//
		//    i ^
		//      |    (1)
		//  (b) +----->-----+ (c)
		//      |           |
		//                  |
		//      |           V (2)
		//                  |
		//      |           |
		//  (a) +-----<-----+--> j
		//      r    (3)   (d)
		//

		// non-topo staple
		equal(&link1, &(GC->lattice[rpi][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rpj][i])); // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));   // link3 = (3)

		times_dag2(&link12, &link1, &link2); // link12 = link1 * link2^{dag}
		times_dag2(&stap, &link12, &link3);  // stap = link12 * link3^{dag}

		// clover insertion in (a)
		times(&aux, &stap, CAr);    // stap * clover
		times_equal_real(&aux, Zr); // bulk factor
		equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, CArpi, &stap);    // clover * stap
		times_equal_real(&aux, Zrpi); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link1, CArpipj);   // link1 * clover
		times_equal_dag(&aux, &link2);  // *= link2^{dag}
		times_equal_dag(&aux, &link3);  // *= link3^{dag}
		times_equal_real(&aux, Zrpipj); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times(&aux, &link12, CArpj);   // link1 * link2 * clover
		times_equal_dag(&aux, &link3); // *= link3^{dag}
		times_equal_real(&aux, Zrpj);  // bulk factor
		plus_equal(&topo_stap, &aux);

		// boundary condition (only affects non-topo staple) and twist (clockwise plaquette) modification
		factor = GC->Z[r][dirs_to_si(i, j)];                                         // Z_\mu\nu(x)
		times_equal_complex(&topo_stap, factor);                                     // Z_\mu\nu(x) * topo_staple
		factor *= (GC->C[r][i]) * (GC->C[rpi][j]) * (GC->C[rpj][i]) * (GC->C[r][j]); // *= K_\mu\nu(x)
		times_equal_complex(&stap, factor);                                          // K_\mu\nu(x) * Z_\mu\nu(x) * staple

		// accumulate M and topo_M
		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);

		//
		//    i ^
		//     |    (1)
		// (d) +----<------+ (a)
		//     |           |
		//     |
		// (2) V           |
		//     |
		//     |           | (b)
		// (c) +------>----+---> j
		//           (3)   r
		//

		// non-topo staple
		equal(&link1, &(GC->lattice[rpimj][j])); // link1 = (1)
		equal(&link2, &(GC->lattice[rmj][i]));   // link2 = (2)
		equal(&link3, &(GC->lattice[rmj][j]));   // link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12 = link1^{dag} * link2^{dag}
		times(&stap, &link12, &link3);        // stap = link12 * link3

		// clover insertion in (a)
		times(&aux, CArpi, &stap);    // clover * stap
		times_equal_real(&aux, Zrpi); // bulk factor
		equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &stap, CAr);    // stap * clover
		times_equal_real(&aux, Zr); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link12, CArmj);  // link1^{dag} * link2^{dag} * clover
		times_equal(&aux, &link3);    // *= link3
		times_equal_real(&aux, Zrmj); // bulk factor
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times_dag1(&aux, &link1, CArpimj); // link1^{dag} * clover
		times_equal_dag(&aux, &link2);     // *= link2^{dag}
		times_equal(&aux, &link3);         // *= link3
		times_equal_real(&aux, Zrpimj);    // bulk factor
		plus_equal(&topo_stap, &aux);

		// boundary condition (only affects non-topo staple) and twist (anticlockwise plaquette) modification
		factor = GC->Z[rmj][dirs_to_si(j, i)];                                           // Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		times_equal_complex(&topo_stap, factor);                                         // Z_\nu\mu(x-\nu) * topo_staple
		factor *= (GC->C[rmj][i]) * (GC->C[rpimj][j]) * (GC->C[r][i]) * (GC->C[rmj][j]); // *= K_\mu\nu(x-\nu)
		times_equal_complex(&stap, factor);                                              // K_\mu\nu(x-\nu) * Z_\nu\mu(x-\nu) * staple

		// accumulate M and topo_M (with minus sign from definition of topological charge)
		plus_equal(M, &stap);
		minus_equal(&topo_M, &topo_stap);
		}

	// multiply topo_staple by theta coefficient and sum to staple
	times_equal_real(&topo_M, topo_staple_coeff);
	plus_equal(M, &topo_M);

	#endif
	}


// perform an update with heatbath
void heatbath(Gauge_Conf *const GC,
              Geometry const *const geo,
              GParam const *const param,
              long r,
              int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap);
	#endif

	single_heatbath(&(GC->lattice[r][i]), &stap, param);
	}


// perform an update with overrelaxation
void overrelaxation(Gauge_Conf *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    long r,
                    int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap);
	#endif

	single_overrelaxation(&(GC->lattice[r][i]), &stap);
	}


// Compute the variation of the Wilson action for the proposed conf (debug only)
double delta_action(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param)
	{
	Gauge_Conf aux;
	aux.lattice = GC->lattice_copy;
	aux.Z = GC->Z_copy;

	double res1 = 0.0, res2 = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : res1, res2)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		GAUGE_GROUP stap1, stap2;
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(&aux, geo, param, r, i, &stap1);
		calcstaples_wilson(GC, geo, param, r, i, &stap2);

		// compute action
		times_equal(&stap1, &(aux.lattice[r][i]));
		times_equal(&stap2, &(GC->lattice[r][i]));
		res1 += retr(&stap1);
		res2 += retr(&stap2);
		}

	return (param->d_beta) * (res1 - res2) * 0.25;
	}


// perform an update with metropolis
// return 1 if the proposed update is accepted
int metropolis(Gauge_Conf *const GC,
               Geometry const *const geo,
               GParam const *const param,
               long r,
               int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap, new_link, tmp_matrix, rnd_matrix;
	double action_new, action_old;

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap);
	#endif

	// compute old action
	times(&tmp_matrix, &(GC->lattice[r][i]), &stap);
	action_old = param->d_beta * (1.0 - retr(&tmp_matrix));

	// compute the new link
	one(&tmp_matrix);
	rand_matrix(&rnd_matrix);
	times_equal_real(&rnd_matrix, param->d_epsilon_metro);
	plus_equal(&rnd_matrix, &tmp_matrix);
	unitarize(&rnd_matrix); // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
	if(casuale() < 0.5)
		{
		times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
		}
	else
		{
		times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
		}

	// new action
	times(&tmp_matrix, &new_link, &stap);
	action_new = param->d_beta * (1.0 - retr(&tmp_matrix));

	int acc = 0;
	if(casuale() < exp(action_old - action_new))
		{
		equal(&(GC->lattice[r][i]), &new_link);
		acc = 1;
		}

	return acc;
	}


// perform an update with metropolis with trace deformations
// return 1 if the proposed update is accepted
int metropolis_with_tracedef(Gauge_Conf *const GC,
                             Geometry const *const geo,
                             GParam const *const param,
                             long r,
                             int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap_w, stap_td, new_link, tmp_matrix, rnd_matrix, poly;
	double action_new, action_old;
	double rpart, ipart;

	// compute old action
	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap_w);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap_w);
	#endif
	times(&tmp_matrix, &(GC->lattice[r][i]), &stap_w);
	action_old = param->d_beta * (1.0 - retr(&tmp_matrix));
	if(i == 0) // just if we are updating a temporal link
		{
		// "staple" for trace deformation
		calcstaples_tracedef(GC, geo, param, r, i, &stap_td);

		// trace deformation contribution to action_old
		times(&poly, &(GC->lattice[r][i]), &stap_td);
		one(&tmp_matrix);
		for(int j = 0; j < (int) floor(NCOLOR / 2.0); j++)
			{
			times_equal(&tmp_matrix, &poly);
			rpart = NCOLOR * retr(&tmp_matrix);
			ipart = NCOLOR * imtr(&tmp_matrix);
			action_old += param->d_h[j] * (rpart * rpart + ipart * ipart);
			}
		}

	// compute the update to be proposed
	one(&tmp_matrix);
	rand_matrix(&rnd_matrix);
	times_equal_real(&rnd_matrix, param->d_epsilon_metro);
	plus_equal(&rnd_matrix, &tmp_matrix);
	unitarize(&rnd_matrix); // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
	if(casuale() < 0.5)
		{
		times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
		}
	else
		{
		times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
		}

	// compute the new action
	times(&tmp_matrix, &new_link, &stap_w);
	action_new = param->d_beta * (1.0 - retr(&tmp_matrix));
	if(i == 0) // just if we are updating a temporal link
		{
		// trace deformation contribution to action_new
		times(&poly, &new_link, &stap_td);
		one(&tmp_matrix);
		for(int j = 0; j < (int) floor(NCOLOR / 2.0); j++)
			{
			times_equal(&tmp_matrix, &poly);
			rpart = NCOLOR * retr(&tmp_matrix);
			ipart = NCOLOR * imtr(&tmp_matrix);
			action_new += param->d_h[j] * (rpart * rpart + ipart * ipart);
			}
		}

	int acc = 0;
	if(casuale() < exp(action_old - action_new))
		{
		equal(&(GC->lattice[r][i]), &new_link);
		acc = 1;
		}

	return acc;
	}


// perform an update with heatbath in the presence of a defect
void heatbath_with_defect(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          long r,
                          int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
	#endif

	single_heatbath(&(GC->lattice[r][i]), &stap, param);
	}


// perform an update with overrelaxation
void overrelaxation_with_defect(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                long r,
                                int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	#endif

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
	#endif

	single_overrelaxation(&(GC->lattice[r][i]), &stap);
	}


// perform a complete update
void update(Gauge_Conf *const GC,
            Geometry const *const geo,
            GParam const *const param,
            Acc_Utils *acc_counters)
	{
	#ifdef DEBUG
	ASSERT(param->d_min_size > 1, "this function cannot be used in the completely reduced case");
	#endif

	#ifndef MULTICANONICAL_MODE
	(void) acc_counters; // to avoid compiler warning of unused variable
	#endif

	/* Check if there's at least one even dimension of the lattice, i.e. check if d_volume is even.
	If there's at least one even dimension: d_volume/2 even sites and d_volume/2 odd sites.
	Otherwise: (d_volume+1)/2 even sites and (d_volume-1)/2 odd sites. */
	long const is_even = (param->d_volume) % 2;
	long const num_even = (param->d_volume + is_even) / 2; // number of even sites

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
			heatbath(GC, geo, param, s, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = num_even; s < (param->d_volume); s++)
			{
			heatbath(GC, geo, param, s, dir);
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
				overrelaxation(GC, geo, param, s, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = num_even; s < (param->d_volume); s++)
				{
				overrelaxation(GC, geo, param, s, dir);
				}
			}
		}

	// Metropolis test if using multicanonical
	int acc = 1;
	#ifdef MULTICANONICAL_MODE
	acc = multicanonic_metropolis_step_all_links(GC, geo, param);
	acc_counters->num_accepted_metro_multicanonic[0] += acc;
	acc_counters->num_metro_multicanonic[0] += 1;
	#endif
	// update or restore the lattice and auxiliary copies
	if(acc == 1) accept_gauge_conf(GC, param);
	else restore_gauge_conf(GC, param);

	GC->update_index++;
	}

// update all replica in the presence of a defect
void update_with_defect(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                        Acc_Utils *acc_counters)
	{
	#ifdef DEBUG
	ASSERT(param->d_min_size > 1, "this function cannot be used in the completely reduced case");
	#endif

	#ifndef MULTICANONICAL_MODE
	(void) acc_counters; // to avoid compiler warning of unused variable
	#endif

	/* Check if there's at least one even dimension of the lattice, i.e. check if d_volume is even.
	If there's at least one even dimension: d_volume/2 even sites and d_volume/2 odd sites.
	Otherwise: (d_volume+1)/2 even sites and (d_volume-1)/2 odd sites. */
	long const is_even = (param->d_volume) % 2;
	long const num_even = (param->d_volume + is_even) / 2; // number of even sites
	long const num_odd = (param->d_volume - is_even) / 2;  // number of odd sites

	// heatbath
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = 0; s < ((param->d_N_replica_pt) * num_even); s++)
			{
			// s = i * num_even + r
			long const r = s % num_even;              // site index
			int const i = (int) ((s - r) / num_even); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = 0; s < ((param->d_N_replica_pt) * num_odd); s++)
			{
			// s = i * num_odd + aux ; aux = r - num_even
			long const aux = s % num_odd;
			long const r = num_even + aux;             // site index
			int const i = (int) ((s - aux) / num_odd); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		}

	// overrelax
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif

		for(int j = 0; j < param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = 0; s < (param->d_N_replica_pt) * num_even; s++)
				{
				// s = i * num_even + r
				long const r = s % num_even;              // site index
				int const i = (int) ((s - r) / num_even); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = 0; s < (param->d_N_replica_pt) * num_odd; s++)
				{
				// s = i * num_odd + aux ; aux = r - num_even
				long const aux = s % num_odd;
				long const r = num_even + aux;             // site index
				int const i = (int) ((s - aux) / num_odd); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}

	// Metropolis test if using multicanonical
	for(int j = 0; j < (param->d_N_replica_pt); j++)
		{
		int acc = 1;
		// multicanonic Metropolis tests and acceptance counters update
		#ifdef MULTICANONICAL_MODE
		acc = multicanonic_metropolis_step_all_links(&(GC[j]), geo, param);
		acc_counters->num_accepted_metro_multicanonic[j] += acc;
		acc_counters->num_metro_multicanonic[j] += 1;
		#endif
		// update or restore the lattice and auxiliary copies
		if(acc == 1) accept_gauge_conf(&(GC[j]), param);
		else restore_gauge_conf(&(GC[j]), param);
		}
	}

// hierarchical update functions

// update all replica only on a given rectangle in the presence of a defect
void update_rectangle_with_defect(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                                  int const hierarc_level,
                                  Rect_Utils const *const rect_aux,
                                  Acc_Utils *acc_counters)
	{
	#ifndef MULTICANONICAL_MODE
	(void) acc_counters; // to avoid compiler warning of unused variable
	#endif

	/* Check if there's at least one even dimension of the rectangle, i.e. check if d_vol_rect is even.
		If there's at least one even dimension: d_vol_rect/2 even sites and d_vol_rect/2 odd sites.
		Otherwise: (d_vol_rect+1)/2 even sites and (d_vol_rect-1)/2 odd sites. */
	long const rect_volume = (rect_aux->update_rect[hierarc_level]).d_vol_rect;
	long const is_even = rect_volume % 2;
	long const num_even = (rect_volume + is_even) / 2; // number of even sites
	long const num_odd = (rect_volume - is_even) / 2;  // number of odd sites

	// heatbath
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, &(rect_aux->clover_rect[hierarc_level]));
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = 0; s < (num_even * (param->d_N_replica_pt)); s++)
			{
			// s = i * num_even + n
			long const n = s % num_even;                                         // site index on rectangle
			long const r = (rect_aux->update_rect[hierarc_level]).rect_sites[n]; // site index on lattice
			int const i = (int) ((s - n) / num_even);                            // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long s = 0; s < (num_odd * (param->d_N_replica_pt)); s++)
			{
			// s = i * num_odd + aux; aux = n - num_even
			long const aux = s % num_odd;
			long const n = aux + num_even;                                       // site index on rectangle
			long const r = (rect_aux->update_rect[hierarc_level]).rect_sites[n]; // site index on lattice
			int const i = (int) ((s - aux) / num_odd);                           // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		}

	// overrelax
	for(int dir = 0; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, &(rect_aux->clover_rect[hierarc_level]));
		#endif

		for(int j = 0; j < param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = 0; s < (num_even * (param->d_N_replica_pt)); s++)
				{
				// s = i * num_even + n
				long const n = s % num_even;                                         // site index on rectangle
				long const r = (rect_aux->update_rect[hierarc_level]).rect_sites[n]; // site index on lattice
				int const i = (int) ((s - n) / num_even);                            // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long s = 0; s < (num_odd * (param->d_N_replica_pt)); s++)
				{
				// s = i * num_odd + aux; aux = n - num_even
				long const aux = s % num_odd;
				long const n = aux + num_even;                                       // site index on rectangle
				long const r = (rect_aux->update_rect[hierarc_level]).rect_sites[n]; // site index on lattice
				int const i = (int) ((s - aux) / num_odd);                           // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}

	// Metropolis test if using multicanonical
	// TODO: check rectangle version (seems ok)
	for(int j = 0; j < (param->d_N_replica_pt); j++)
		{
		int acc = 1;
		// multicanonic Metropolis tests and acceptance counters update
		#ifdef MULTICANONICAL_MODE
		//acc = multicanonic_metropolis_step_all_links(&(GC[j]), geo, param);
		acc = multicanonic_metropolis_step_rectangle(&(GC[j]), geo, param, hierarc_level, rect_aux);
		acc_counters->num_accepted_metro_multicanonic[j] += acc;
		acc_counters->num_metro_multicanonic[j] += 1;
		#endif
		// update or restore the lattice and auxiliary copies
		if(acc == 1) accept_gauge_conf_rectangle(&(GC[j]), hierarc_level, rect_aux);
		else restore_gauge_conf_rectangle(&(GC[j]), hierarc_level, rect_aux);
		//if (acc == 1) accept_gauge_conf(&(GC[j]), param);
		//else restore_gauge_conf(&(GC[j]), param);
		}
	}

// perform a hierarchical update on all rectangles
void hierarchical_update_rectangle_with_defect(Gauge_Conf *const GC, Geometry const *const geo,
                                               GParam const *const param,
                                               int const hierarc_level,
                                               Rect_Utils const *const rect_aux,
                                               Acc_Utils *acc_counters)
	{
	if(hierarc_level < param->d_N_hierarc_levels)
		{
		for(int j = 0; j < param->d_N_sweep_rect[hierarc_level]; j++)
			{
			update_rectangle_with_defect(GC, geo, param, hierarc_level, rect_aux, acc_counters);
			if(param->d_N_replica_pt > 1)
				{
				swap(GC, geo, param, &(rect_aux->swap_rect), acc_counters);
				if(fabs(param->d_pt_bound_cond_coeff[0] - 1.0) < MIN_VALUE)
					conf_translation(&(GC[0]), geo, param);
				}
			hierarchical_update_rectangle_with_defect(GC, geo, param, hierarc_level + 1, rect_aux, acc_counters);
			}
		}
	}

// perform a single step of parallel tempering with hierarchical update
void parallel_tempering_with_hierarchical_update(Gauge_Conf *const GC, Geometry const *const geo,
                                                 GParam const *const param,
                                                 Rect_Utils const *const rect_aux,
                                                 Acc_Utils *acc_counters)
	{
	// set multicanonic Metropolis acceptance counters to zero to compute mean acc over single updating step
	#ifdef MULTICANONICAL_MODE
	for(int i = 0; i < param->d_N_replica_pt; i++)
		{
		acc_counters->num_accepted_metro_multicanonic[i] = 0;
		acc_counters->num_metro_multicanonic[i] = 0;
		}
	#endif

	// full update + hierarchical update + swaps and translations after every sweep
	update_with_defect(GC, geo, param, acc_counters);
	if(param->d_N_replica_pt > 1)
		{
		swap(GC, geo, param, &(rect_aux->swap_rect), acc_counters);
		if(fabs(param->d_pt_bound_cond_coeff[0] - 1.0) < MIN_VALUE)
			conf_translation(&(GC[0]), geo, param);
		}
	hierarchical_update_rectangle_with_defect(GC, geo, param, 0, rect_aux, acc_counters);

	// increase update index of all replica
	for(int i = 0; i < param->d_N_replica_pt; i++)
		GC[i].update_index++;

	// print mean multicanonic acceptance over a single updating step
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_acceptance(GC, param, acc_counters);
	#endif
	}


// perform a complete update with trace deformation
// TODO: check if ok with multicanonical
void update_with_trace_def(Gauge_Conf *const GC,
                           Geometry const *const geo,
                           GParam const *const param,
                           double *acc_td)
	{
	long const num_even = (param->d_volume + (param->d_volume % 2)) / 2;
	long const num_sp_even = (param->d_space_vol[0] + (param->d_space_vol[0] % 2)) / 2;

	int *a;
	allocate_array_int(&a, param->d_space_vol[0], __FILE__, __LINE__);
	for(long r = 0; r < param->d_space_vol[0]; r++) a[r] = 0;

	// heatbath on spatial links
	for(int dir = 1; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < num_even; r++)
			{
			heatbath(GC, geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = num_even; r < (param->d_volume); r++)
			{
			heatbath(GC, geo, param, r, dir);
			}
		}

	// metropolis on temporal links
	for(int t = 0; t < param->d_size[0]; t++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, 0);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < num_sp_even; r++)
			{
			long const r4 = sisp_and_t_to_si(geo, r, t);
			a[r] += metropolis_with_tracedef(GC, geo, param, r4, 0);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = num_sp_even; r < (param->d_space_vol[0]); r++)
			{
			long const r4 = sisp_and_t_to_si(geo, r, t);
			a[r] += metropolis_with_tracedef(GC, geo, param, r4, 0);
			}
		}

	long asum = 0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for reduction(+:asum)
	#endif
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		asum += (long) a[r];
		}
	*acc_td = ((double) asum) * param->d_inv_vol;

	// overrelax spatial links
	for(int dir = 1; dir < STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		for(int j = 0; j < param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < num_even; r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = num_even; r < (param->d_volume); r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}
			}
		}

	// Metropolis test if using multicanonical and final unitarization
	int acc = 1;

	#ifdef MULTICANONICAL_MODE
	acc = multicanonic_metropolis_step_all_links(GC, geo, param);
	#endif

	// update or restore the lattice auxiliary copies
	if(acc == 1) accept_gauge_conf(GC, param);
	else restore_gauge_conf(GC, param);

	free(a);

	GC->update_index++;
	}


// perform n cooling steps minimizing the action at theta=0
void cooling(Gauge_Conf *const GC,
             Geometry const *const geo,
             GParam const *const param,
             int const n)
	{
	long const num_even = (param->d_volume + (param->d_volume % 2)) / 2;
	for(int k = 0; k < n; k++)
		{
		// cooling
		for(int dir = 0; dir < STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < num_even; r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = num_even; r < (param->d_volume); r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));
		unitarize(&(GC->lattice[r][dir]));
		}
	}


// perform n cooling steps minimizing the action at theta=0 in the presence of the defect
void cooling_with_defect(Gauge_Conf *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         int const n)
	{
	long const num_even = (param->d_volume + (param->d_volume % 2)) / 2;

	for(int k = 0; k < n; k++)
		{
		// cooling
		for(int dir = 0; dir < STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < num_even; r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson_with_defect(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = num_even; r < (param->d_volume); r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson_with_defect(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));
		unitarize(&(GC->lattice[r][dir]));
		}
	}


// perform topo_coolsteps cooling steps around the defect minimizing the action at theta=0
void hierarchical_cooling(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          Rectangle const *const cooling_rect)
	{
	for(int k = 0; k < param->d_topo_coolsteps; k++)
		{
		// cooling
		long const is_even = ((cooling_rect[k]).d_vol_rect) % 2;
		long const num_even = ((cooling_rect[k]).d_vol_rect + is_even) / 2; // number of even sites

		for(int dir = 0; dir < STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long n = 0; n < num_even; n++)
				{
				GAUGE_GROUP staple;
				long const r = (cooling_rect[k]).rect_sites[n];
				calcstaples_wilson(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long n = num_even; n < ((cooling_rect[k]).d_vol_rect); n++)
				{
				GAUGE_GROUP staple;
				long const r = (cooling_rect[k]).rect_sites[n];
				calcstaples_wilson(GC, geo, param, r, dir, &staple);
				cool(&(GC->lattice[r][dir]), &staple);
				}
			}
		}

	// final unitarization
	int k = param->d_topo_coolsteps - 1; // largest rectangle
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long n = 0; n < (STDIM * ((cooling_rect[k]).d_vol_rect)); n++)
		{
		long const s = n % ((cooling_rect[k]).d_vol_rect);
		long const r = (cooling_rect[k]).rect_sites[s];
		int const dir = (int) ((n - s) / ((cooling_rect[k]).d_vol_rect));
		unitarize(&(GC->lattice[r][dir]));
		}
	}


// perform a single step of the Runge Kutta integrator for the Wilson flow
// as described in Luscher arXiv:1006.4518 app. C
void gradflow_RKstep(Gauge_Conf *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     double dt,
                     Meas_Utils *meas_aux)
	{
	GAUGE_GROUP staple, aux, link;

	// initialize
	equal_lattice(meas_aux->lattice_aux[0], (GAUGE_GROUP const * const *) GC->lattice, param);

	// just to call calcstaples_wilson on aux lattice 0
	Gauge_Conf helper;
	helper.lattice = meas_aux->lattice_aux[0];
	helper.Z = GC->Z;

	// now GC = lattice0 = W_0, lattice1 = uninitialized
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(&helper, geo, param, r, dir, &staple); // staple   = staple(W_0)
		equal(&link, &(meas_aux->lattice_aux[0][r][dir]));        // link     = link(W_0)
		times(&aux, &link, &staple);                              // aux      = force(W_0)
		times_equal_real(&aux, -dt / 4.0);                        // aux      = -1/4*dt*force(W_0) = 1/4*Z_0
		equal(&(meas_aux->lattice_aux[1][r][dir]), &aux);         // lattice1 = 1/4*Z_0
		taexp(&aux);                                              // aux      = exp(1/4*Z_0)
		times(&(GC->lattice[r][dir]), &aux, &link);               // GC       = exp(1/4*Z_0)*W_0 = W_1
		}

	// now GC = W_1, lattice0 = W_0, lattice1 = 1/4*Z_0
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(GC, geo, param, r, dir, &staple);                           // staple   = staple(W_1)
		equal(&link, &(GC->lattice[r][dir]));                                          // link     = link(W_1)
		times(&aux, &link, &staple);                                                   // aux      = force(W_1)
		times_equal_real(&aux, -dt * 8.0 / 9.0);                                       // aux      = -8/9*dt*force(W_1) = 8/9*Z_1
		minus_equal_times_real(&aux, &(meas_aux->lattice_aux[1][r][dir]), 17.0 / 9.0); // aux      = 8/9*Z_1 - 17/36*Z_0
		equal(&(meas_aux->lattice_aux[1][r][dir]), &aux);                              // lattice1 = 8/9*Z_1 - 17/36*Z_0
		taexp(&aux);                                                                   // aux      = exp(8/9*Z_1 - 17/36*Z_0)
		times(&(meas_aux->lattice_aux[0][r][dir]), &aux, &link);                       // lattice0 = exp(8/9*Z_1 - 17/36*Z_0)*W_1 = W_2
		}

	// now GC = W_1, lattice0 = W_2, lattice1 = 8/9*Z_1-17/36*Z_0
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(&helper, geo, param, r, dir, &staple); // staple = staple(W_2)
		equal(&link, &(meas_aux->lattice_aux[0][r][dir]));        // link   = link(W_2)
		times(&aux, &link, &staple);                              // aux    = force(W_2)
		times_equal_real(&aux, -dt * 3.0 / 4.0);                  // aux    = -3/4*dt*force(W_2) = 3/4*Z_2
		minus_equal(&aux, &(meas_aux->lattice_aux[1][r][dir]));   // aux    = 3/4*Z_2 - 8/9*Z_1 + 17/36*Z_0
		taexp(&aux);                                              // aux    = exp(3/4*Z_2 - 8/9*Z_1 + 17/36*Z_0)
		times(&(GC->lattice[r][dir]), &aux, &link);               // GC     = exp(3/4*Z_2 - 8/9*Z_1 + 17/36*Z_0)*W_2 = W_3
		unitarize(&(GC->lattice[r][dir]));
		}
	// now GC = W_3, lattice0 = W_2, lattice1 = 8/9*Z_1-17/36*Z_0
	}


double gradflow_RKstep_adaptive_aux(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double dt,
                                    Meas_Utils *meas_aux)
	{
	GAUGE_GROUP staple, aux, aux2, link;

	// initialize
	equal_equal_lattice(meas_aux->lattice_aux[0], meas_aux->lattice_aux[3], (GAUGE_GROUP const * const *) GC->lattice, param);
	for(int j = 0; j < NTHREADS; j++) meas_aux->local_max_dist[j] = 0.0;

	// just to call calcstaples_wilson on aux lattice 0
	Gauge_Conf helper;
	helper.lattice = meas_aux->lattice_aux[0];
	helper.Z = GC->Z;

	// now GC = lattice0 = lattice3 = W_0, lattice1 = lattice2 = uninitialized
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(&helper, geo, param, r, dir, &staple); // staple   = staple(W_0)
		equal(&link, &(meas_aux->lattice_aux[0][r][dir]));        // link     = link(W_0)
		times(&aux, &link, &staple);                              // aux      = force(W_0)
		times_equal_real(&aux, -dt);                              // aux      = -dt*force(W_0) = Z_0
		equal(&(meas_aux->lattice_aux[1][r][dir]), &aux);         // lattice1 = Z_0
		times_equal_real(&aux, 1.0 / 4.0);                        // aux      = 1/4*Z_0
		taexp(&aux);                                              // aux      = exp(1/4*Z_0)
		times(&(GC->lattice[r][dir]), &aux, &link);               // GC       = exp(1/4*Z_0)*W_0 = W_1
		}

	// now GC = W_1, lattice0 = lattice3 = W_0, lattice1 = Z_0, lattice2 = uninitialized
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, aux2, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(GC, geo, param, r, dir, &staple); // staple = staple(W_1)
		equal(&link, &(GC->lattice[r][dir]));                // link   = link(W_1)
		times(&aux, &link, &staple);                         // aux    = force(W_1)
		times_equal_real(&aux, -dt);                         // aux    = -dt*force(W_1) = Z_1
		equal(&aux2, &aux);                                  // aux2   = Z_1

		times_equal_real(&aux2, 2.0);                                                           // aux2   = 2*Z_1
		minus_equal(&aux2, &(meas_aux->lattice_aux[1][r][dir]));                                // aux2   = 2*Z_1-Z_0
		taexp(&aux2);                                                                           // aux2   = exp(2*Z_1-Z_0)
		times(&(meas_aux->lattice_aux[2][r][dir]), &aux2, &(meas_aux->lattice_aux[0][r][dir])); // lattice2 = exp(2*Z_1-Z_0)*W_0

		times_equal_real(&aux, 8.0 / 9.0);                                              // aux      = 8/9*Z_1
		minus_equal_times_real(&aux, &(meas_aux->lattice_aux[1][r][dir]), 17.0 / 36.0); // aux      = 8/9*Z_1 - 17/36*Z_0
		equal(&(meas_aux->lattice_aux[1][r][dir]), &aux);                               // lattice1 = 8/9*Z_1 - 17/36*Z_0
		taexp(&aux);                                                                    // aux      = exp(8/9*Z_1 - 17/36*Z_0)
		times(&(meas_aux->lattice_aux[0][r][dir]), &aux, &link);                        // lattice0 = exp(8/9*Z_1 - 17/36*Z_0)*W_1 = W_2
		}

	// now GC = W_1, lattice0 = W_2, lattice1 = (8/9)Z_1-(17/36)Z_0, lattice2 = W'_2, lattice3 = W_0
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(staple, aux, link)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson(&helper, geo, param, r, dir, &staple); // staple = staple(W_2)
		equal(&link, &(meas_aux->lattice_aux[0][r][dir]));        // link   = link(W_2)
		times(&aux, &link, &staple);                              // aux    = force(W_2)
		times_equal_real(&aux, -dt * 3.0 / 4.0);                  // aux    = -3/4*dt*force(W_2) = 3/4*Z_2
		minus_equal(&aux, &(meas_aux->lattice_aux[1][r][dir]));   // aux    = (3/4)Z_2-(8/9)Z_1+(17/36)Z_0
		taexp(&aux);                                              // aux    = exp((3/4)Z_2-(8/9)Z_1+(17/36)Z_0)
		times(&(GC->lattice[r][dir]), &aux, &link);               // GC     = exp((3/4)Z_2-(8/9)Z_1+(17/36)Z_0)*W_2 = W_3
		}

	// now GC = W_3, lattice0 = W_2, lattice1 = (8/9)Z_1-(17/36)Z_0, lattice2 = W'_2, lattice3 = W_0
	// error calculation, dist(W_3, W'_2), and final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const dir = (int) ((s - r) / (param->d_volume));
		#ifdef OPENMP_MODE
		int const thread_num = omp_get_thread_num();
		#else
		int const thread_num = 0;
		#endif

		minus_equal(&(meas_aux->lattice_aux[2][r][dir]), &(GC->lattice[r][dir]));
		double const dist = norm(&(meas_aux->lattice_aux[2][r][dir]));
		if(dist > meas_aux->local_max_dist[thread_num]) meas_aux->local_max_dist[thread_num] = dist;
		unitarize(&(GC->lattice[r][dir]));
		}

	// reduction to find maximum dist over threads
	double max_dist = MIN_VALUE;
	for(int j = 0; j < NTHREADS; j++)
		{
		if(meas_aux->local_max_dist[j] > max_dist)
			{
			max_dist = meas_aux->local_max_dist[j];
			}
		}

	return max_dist / ((double) NCOLOR * (double) NCOLOR);
	}


// perform a single step of the Runge Kutta integrator for the Wilson flow
// with adaptive integration step as described in Fritzsch-Ramos arXiv:1301.4388 app. D
void gradflow_RKstep_adaptive(Gauge_Conf *const GC,
                              Geometry const *const geo,
                              GParam const *const param,
                              double *t,
                              double *dt,
                              int *accepted,
                              Meas_Utils *meas_aux)
	{
	// integration, distance calculation and unitarization
	double const max_dist = gradflow_RKstep_adaptive_aux(GC, geo, param, *dt, meas_aux);

	// accept-reject step: if the integration step is accepted advance t, else reset gauge conf
	if(max_dist < param->d_agf_delta)
		{
		*accepted = 1;
		*t = *t + *dt;
		}
	else
		{
		*accepted = 0;
		equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) meas_aux->lattice_aux[3], param);
		}

	//TODO: remove, debug only
	//fprintf(stdout, "%ld %d %.12g %.12g %.12g\n", GC->update_index, *accepted, *t, *dt, max_dist);
	//fflush(stdout);

	// calculation of new integration step
	*dt *= 0.95 * pow(param->d_agf_delta / max_dist, 1.0 / 3.0);
	if(*dt > 0.1) *dt = 0.1;
	}

double gradflow_RKstep_adaptive_no_advance(Gauge_Conf *const GC,
                                           Geometry const *const geo,
                                           GParam const *const param,
                                           double const *t,
                                           double const *dt,
                                           double *dt_new,
                                           int *accepted,
                                           Meas_Utils *meas_aux)
	{
	(void) t;

	// integration, distance calculation and unitarization
	double const max_dist = gradflow_RKstep_adaptive_aux(GC, geo, param, *dt, meas_aux);

	// no accept-reject step
	if(max_dist < param->d_agf_delta)
		{
		*accepted = 1;
		}
	else
		{
		*accepted = 0;
		}

	// calculation of new integration step
	*dt_new = *dt * 0.95 * pow(param->d_agf_delta / max_dist, 1.0 / 3.0);

	return max_dist;
	}


void gradflow_RKstep_adaptive_check(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double *t,
                                    double *dt,
                                    int *accepted,
                                    Meas_Utils *meas_aux,
                                    Gauge_Conf *const GC_reset,
                                    Gauge_Conf *const conf_rk3dth)
	{
	double dt_new, dth;

	dth = *(dt) / 2;
	equal_lattice(GC_reset->lattice, (GAUGE_GROUP const * const *) GC->lattice, param);
	double const dist_rk2dth_rk3dth_1 = gradflow_RKstep_adaptive_no_advance(GC, geo, param, t, &dth, &dt_new, accepted, meas_aux);
	double const dist_rk2dth_rk3dth_2 = gradflow_RKstep_adaptive_no_advance(GC, geo, param, t, &dth, &dt_new, accepted, meas_aux);
	equal_lattice(conf_rk3dth->lattice, (GAUGE_GROUP const * const *) GC->lattice, param);
	equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) GC_reset->lattice, param);
	double const dist_rk2dt_rk3dt = gradflow_RKstep_adaptive_no_advance(GC, geo, param, t, dt, &dt_new, accepted, meas_aux);
	double const dist_rk3dth_rk3dt = lattice_max_dist((GAUGE_GROUP const * const *) GC->lattice, (GAUGE_GROUP const * const *) conf_rk3dth->lattice, param);

	fprintf(stdout, "%ld %d %18.12e %18.12e %18.12e\n", GC->update_index, *accepted, *t, *dt, dist_rk2dt_rk3dt);
	fprintf(stdout, "%ld %d %18.12g %18.12g %18.12g %18.12g %18.12g %18.12g\n", GC->update_index, *accepted, *t, *dt, dist_rk2dt_rk3dt, dist_rk2dth_rk3dth_1, dist_rk2dth_rk3dth_2, dist_rk3dth_rk3dt);
	fflush(stdout);

	//if (*accepted == 1)
	//	{
	//	*t += *dt;
	//	}
	//else
	//	{
	//	equal_lattice(GC->lattice, GC_reset->lattice, param);
	//	}
	//*dt = dt_new;
	}


void gradflow_RKstep_adaptive_debug(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double *t,
                                    double *dt,
                                    int *accepted,
                                    double *total_error,
                                    Meas_Utils *meas_aux)
	{
	// integration, distance calculation and unitarization
	double const max_dist = gradflow_RKstep_adaptive_aux(GC, geo, param, *dt, meas_aux);
	*total_error = max_dist;

	//if the integration step is accepted, advance t and reset dt. Else, reset gauge conf and decrease dt
	if(max_dist < param->d_agf_delta && *dt < 1.01 * param->d_agf_meas_each)
		{
		*accepted = 1;
		*t = *t + *dt;
		*dt = param->d_agf_step;
		//*total_error += max_dist;
		}
	else
		{
		*accepted = 0;
		equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) meas_aux->lattice_aux[3], param);
		*dt -= param->d_agf_meas_each;
		}
	}


void gradflow_RKstep_adaptive_debug2(Gauge_Conf *const GC,
                                     Geometry const *const geo,
                                     GParam const *const param,
                                     double *t,
                                     double *dt,
                                     int *accepted,
                                     double *total_error,
                                     Meas_Utils *meas_aux)
	{
	// integration, distance calculation and unitarization
	double const max_dist = gradflow_RKstep_adaptive_aux(GC, geo, param, *dt, meas_aux);
	*total_error = max_dist;

	//if the integration step is accepted, nothing. Else, reset gauge conf
	if(max_dist < param->d_agf_delta) // && *dt < 1.01*param->d_agf_meas_each)
		{
		*accepted = 1;
		(void) t;
		}
	else
		{
		*accepted = 0;
		equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) meas_aux->lattice_aux[3], param);
		//*dt = *dt-param->d_agf_meas_each;
		}
	}


// n step of ape smearing with parameter alpha
// new=Proj[old + alpha *staple ]
void ape_smearing(Gauge_Conf *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  double alpha,
                  int n)
	{
	Gauge_Conf helper1;
	init_gauge_conf_from_gauge_conf(&helper1, GC, param); //helper1=GC

	for(int count = 0; count < n; count++)
		{
		if(count % 2 == 0) // smear(helper1)->GC
			{
			for(int dir = 0; dir < STDIM; dir++)
				{
				#ifdef OPENMP_MODE
				#pragma omp parallel for num_threads(NTHREADS)
				#endif
				for(long r = 0; r < param->d_volume; r++)
					{
					GAUGE_GROUP staple, link;

					calcstaples_wilson(&helper1, geo, param, r, dir, &staple);
					equal(&link, &(helper1.lattice[r][dir]));
					times_equal_real(&link, 1 - alpha);
					times_equal_real(&staple, alpha / 6.0);
					plus_equal_dag(&link, &staple);
					unitarize(&link);
					equal(&(GC->lattice[r][dir]), &link);
					}
				}
			}
		else // smear(GC)->helper1
			{
			for(int dir = 0; dir < STDIM; dir++)
				{
				#ifdef OPENMP_MODE
				#pragma omp parallel for num_threads(NTHREADS)
				#endif
				for(long r = 0; r < param->d_volume; r++)
					{
					GAUGE_GROUP staple, link;

					calcstaples_wilson(GC, geo, param, r, dir, &staple);
					equal(&link, &(GC->lattice[r][dir]));
					times_equal_real(&link, 1 - alpha);
					times_equal_real(&staple, alpha / 6.0);
					plus_equal_dag(&link, &staple);
					unitarize(&link);
					equal(&(helper1.lattice[r][dir]), &link);
					}
				}
			}
		}

	if(n > 0 && n % 2 == 0) equal_gauge_conf(GC, &helper1, param); // GC=helper1
	free_gauge_conf(&helper1, param);
	}


#endif
