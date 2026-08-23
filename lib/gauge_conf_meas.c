#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<complex.h>

#include"../include/gparam.h"
#include"../include/memalign.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"

#include<time.h> // DEBUG


// return 1/Nc*ReTr[P_ji(r)] with P_ji(r) the plaquette in position r with positive directions (j,i)
double plaquettep(Gauge_Conf const *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  long r,
                  int i,
                  int j)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	ASSERT(j < STDIM, "j too large: %d >= %d", j, STDIM);
	#else
	(void) param;
	#endif

	GAUGE_GROUP matrix;

	//
	//       ^ i
	//       |  (2)
	//       +---<---+
	//       |       |
	//   (3) V       ^ (1)
	//       |       |
	//       +--->---+---> j
	//       r  (4)
	//

	// anticlockwise plaquette P_ji
	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));

	//twist factor Z_ji = conj(Z_ij)
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j, i)]);

	return retr(&matrix);
	}


// return 1/Nc*Tr[P_ji(r)] with P_ji(r) the plaquette in position r with positive directions (j,i)
double complex plaquettep_complex(Gauge_Conf const *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  long r,
                                  int j,
                                  int i)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	ASSERT(j < STDIM, "j too large: %d >= %d", j, STDIM);
	#else
	(void) param;
	#endif

	GAUGE_GROUP matrix;

	//
	//       ^ i
	//       |  (2)
	//       +---<---+
	//       |       |
	//   (3) V       ^ (1)
	//       |       |
	//       +--->---+---> j
	//       r  (4)
	//

	// anticlockwise plaquette P_ji
	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));

	//twist factor Z_ji = conj(Z_ji)
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j, i)]);

	return retr(&matrix) + I * imtr(&matrix);
	}


// compute P_ij(r), the plaquette in position r with positive directions (i,j), and save it in matrix
void plaquettep_matrix(Gauge_Conf const *const GC,
                       Geometry const *const geo,
                       GParam const *const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	ASSERT(j < STDIM, "j too large: %d >= %d", j, STDIM);
	#else
	(void) param;
	#endif

	//
	//       ^ j
	//       |  (3)
	//       +---<---+
	//       |       |
	//   (4) V       ^ (2)
	//       |       |
	//       +--->---+---> i
	//       r  (1)
	//

	// anticlockwise plaquette P_ij
	equal(matrix, &(GC->lattice[r][i]));
	times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(matrix, &(GC->lattice[r][j]));

	//twist factor Z_ij
	times_equal_complex(matrix, GC->Z[r][dirs_to_si(i, j)]);
	}


// compute C_ij(r), the four-leaf clover in position r with positive directions (i,j), and save it in M
void clover(Gauge_Conf const *const GC,
            Geometry const *const geo,
            GParam const *const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	ASSERT(r < param->d_volume, "r too large: %ld >= %ld", r, param->d_volume);
	ASSERT(i < STDIM, "i too large: %d >= %d", i, STDIM);
	ASSERT(j < STDIM, "j too large: %d >= %d", j, STDIM);
	#else
	(void) param;
	#endif

	GAUGE_GROUP aux;

	int const si_ij = dirs_to_si(i, j);

	long const rpi = nnp(geo, r, i);
	long const rmi = nnm(geo, r, i);
	long const rpj = nnp(geo, r, j);
	long const rmj = nnm(geo, r, j);
	long const rpimj = nnm(geo, rpi, j);
	long const rmipj = nnp(geo, rmi, j);
	long const rmimj = nnm(geo, rmi, j);

	//
	//                    j ^
	//                      |
	//               (14)   |     (3)
	//          +-----<-----++-----<-----+
	//          |           ||           |
	//          |           ||           |
	//     (15) V      (13) ^V (4)       ^ (2)
	//          |           ||           |
	//          |    (16)   || r  (1)    |
	//      rmi +----->-----++----->-----+------> i
	//          +-----<-----++-----<-----+
	//          |    (9)    ||    (8)    |
	//          |           ||           |
	//     (10) V      (12) ^V (5)       ^ (7)
	//          |           ||           |
	//          |           ||           |
	//          +------>----++----->-----+
	//                (11)   rmj  (6)
	//

	// P_ij(r)
	equal(&aux, &(GC->lattice[r][i]));             // 1
	times_equal(&aux, &(GC->lattice[rpi][j]));     // 2
	times_equal_dag(&aux, &(GC->lattice[rpj][i])); // 3
	times_equal_dag(&aux, &(GC->lattice[r][j]));   // 4
	times_equal_complex(&aux, GC->Z[r][si_ij]);    // Z_ij
	equal(M, &aux);

	// P_-ji(r)
	equal_dag(&aux, &(GC->lattice[rmj][j]));      // 5
	times_equal(&aux, &(GC->lattice[rmj][i]));    // 6
	times_equal(&aux, &(GC->lattice[rpimj][j]));  // 7
	times_equal_dag(&aux, &(GC->lattice[r][i]));  // 8
	times_equal_complex(&aux, GC->Z[rmj][si_ij]); // Z_-ji
	plus_equal(M, &aux);

	// P_-i-j(r)
	equal_dag(&aux, &(GC->lattice[rmi][i]));         // 9
	times_equal_dag(&aux, &(GC->lattice[rmimj][j])); // 10
	times_equal(&aux, &(GC->lattice[rmimj][i]));     // 11
	times_equal(&aux, &(GC->lattice[rmj][j]));       // 12
	times_equal_complex(&aux, GC->Z[rmimj][si_ij]);  // Z_-i-j
	plus_equal(M, &aux);

	// P_j-i(r)
	equal(&aux, &(GC->lattice[r][j]));               // 13
	times_equal_dag(&aux, &(GC->lattice[rmipj][i])); // 14
	times_equal_dag(&aux, &(GC->lattice[rmi][j]));   // 15
	times_equal(&aux, &(GC->lattice[rmi][i]));       // 16
	times_equal_complex(&aux, GC->Z[rmi][si_ij]);    // Z_j-i
	plus_equal(M, &aux);
	}


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const *const GC,
               Geometry const *const geo,
               GParam const *const param,
               double *plaqs,
               double *plaqt)
	{
	double ps = 0.0, pt = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : pt) reduction(+ : ps)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		for(int j = 1; j < STDIM; j++)
			{
			pt += plaquettep(GC, geo, param, r, 0, j);
			}

		for(int i = 1; i < STDIM; i++)
			{
			for(int j = i + 1; j < STDIM; j++)
				{
				ps += plaquettep(GC, geo, param, r, i, j);
				}
			}
		}

	#if STDIM > 2
	double const inv_ns = 2.0 / (((double) STDIM - 1) * (STDIM - 2));
	*plaqs = ps * param->d_inv_vol * inv_ns;
	#else
	*plaqs = 0.0;
	#endif

	#if STDIM > 1
	double const inv_nt = 1.0 / ((double) STDIM - 1);
	*plaqt = pt * param->d_inv_vol * inv_nt;
	#else
	*plaqt = 0.0;
	#endif
	}


// compute the clover discretization of
// sum_{\mu\nu} Tr(F_{\mu\nu}F_{\mu\nu})/2
double loc_clover_energy(Gauge_Conf const *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         long const r)
	{
	double res = 0.0;

	for(int i = 0; i < STDIM; i++)
		{
		for(int j = i + 1; j < STDIM; j++)
			{
			GAUGE_GROUP aux1, aux2;
			clover(GC, geo, param, r, i, j, &aux1);
			ta(&aux1);
			equal(&aux2, &aux1);
			times_equal(&aux1, &aux2);
			res += retr(&aux1);
			}
		}

	return -(double) NCOLOR * res / 16;
	}


// compute and write to binary file the clover energy density
void energy_density(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    Meas_Utils *const meas_aux,
                    int const meas_count)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		meas_aux->scalar_density[r] = loc_clover_energy(GC, geo, param, r);;
		//TODO: try parallel writing with pwrite() if available
		}

	FILE *fp = meas_aux->energydensityfilep;
	fwrite(&(GC->update_index), sizeof(long), 1, fp);
	fwrite(&meas_count, sizeof(int), 1, fp);
	fwrite(meas_aux->scalar_density, sizeof(double), (size_t) param->d_volume, fp);
	}


// compute the average clover energy density
void clover_disc_energy(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        double *const energy)
	{
	double res = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : res)
	#endif
	for(long r = 0; r < param->d_volume; r++) res += loc_clover_energy(GC, geo, param, r);
	*energy = res * param->d_inv_vol;
	}


// compute the clover energy density
// for all the slices along direction mu
void clover_energy_slices(Gauge_Conf const *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          int const mu,
                          double *slices,
                          int meas_count,
                          FILE *filep)
	{
	int const L_mu = param->d_size[mu];
	int const offset = (int) (GC->translation[mu] % L_mu);
	for(int i = 0; i < L_mu; i++) slices[i] = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+: slices[:L_mu])
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		int const x_mu = periodic_condition(geo->d_mucomp[mu][r] - offset, L_mu);
		slices[x_mu] += loc_clover_energy(GC, geo, param, r);
		}
	for(int i = 0; i < L_mu; i++) slices[i] *= param->d_inv_space_vol[mu];

	fprintf(filep, "%ld %d ", GC->update_index, meas_count);
	for(int i = 0; i < L_mu; i++) fprintf(filep, " % 18.12e", slices[i]);
	fprintf(filep, "\n");
	}


// compute the total action (debug only)
void action(Gauge_Conf const *const GC,
            Geometry const *const geo,
            GParam const *const param,
            double *S_wilson, double *S_theta, double *S_total, double *V_mc)
	{
	// compute Wilson and theta terms of the action as sum of forces on the links
	double Sw = 0.0, St = 0.0;
	#ifdef THETA_MODE
	for(int dir = 0; dir < STDIM; dir++)
		{
		compute_clovers(GC, geo, param, dir);
		}
	#endif
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : Sw, St, S)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		GAUGE_GROUP staple_wilson, staple_topo;
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));

		calcstaples_wilson_with_defect(GC, geo, param, r, i, &staple_wilson);
		calcstaples_with_topo_with_defect(GC, geo, param, r, i, &staple_topo);
		minus_equal(&staple_topo, &staple_wilson);

		times_equal(&staple_wilson, &(GC->lattice[r][i]));
		times_equal(&staple_topo, &(GC->lattice[r][i]));
		Sw += retr(&staple_wilson);
		St += retr(&staple_topo);
		}
	double const S0 = (double) (param->d_n_planes * param->d_volume) / 2;
	*S_wilson = param->d_beta * (S0 - Sw / 4.0); // Wilson term: degree 4 in the links, sum of forces divided by 4
	*S_theta = -param->d_beta * St / 8.0;        // Theta term: degree 8 in the links, sum of forces divided by 8

	// compute the total action from average plaquette and topological charge
	double plaqs, plaqt;
	plaquette(GC, geo, param, &plaqs, &plaqt);
	double const sum_plaq = (double) param->d_volume * ((STDIM - 1) * (STDIM - 2) * plaqs / 2 + (STDIM - 1) * plaqt);
	#ifdef THETA_MODE
	*S_total = param->d_beta * (S0 - sum_plaq) - param->d_theta * topcharge(GC, geo, param);
	#else
	*S_total = param->d_beta * (S0 - sum_plaq);
	#endif

	// compute multicanonical potential
	#ifdef MULTICANONICAL_MODE
	*V_mc = compute_topo_potential(GC->replica_index, topcharge(GC, geo, param), param);
	#else
	*V_mc = 0.0;
	#endif
	}


// compute the polyakov loop density and write to binary file
void polyakov_density(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param,
                      int mu,
                      Meas_Utils *const meas_aux,
                      int const meas_count)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long rsp = 0; rsp < param->d_space_vol[mu]; rsp++)
		{
		long r = sisp_and_mu_to_si(geo, rsp, 0, mu);
		GAUGE_GROUP matrix;
		one(&matrix);
		for(int i = 0; i < param->d_size[mu]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][mu]));
			r = nnp(geo, r, mu);
			}

		meas_aux->polyre_density[rsp] = retr(&matrix);
		meas_aux->polyim_density[rsp] = imtr(&matrix);
		//TODO: try parallel writing with pwrite() if available
		}

	FILE *fp = meas_aux->polyakovdensityfilep[mu];
	fwrite(&(GC->update_index), sizeof(long), 1, fp);
	fwrite(&meas_count, sizeof(int), 1, fp);
	fwrite(meas_aux->polyre_density, sizeof(double), (size_t) param->d_space_vol[mu], fp);
	fwrite(meas_aux->polyim_density, sizeof(double), (size_t) param->d_space_vol[mu], fp);
	}


// compute the mean trace of the Polyakov loop winding pwr times along direction mu
void polyakov(Gauge_Conf const *const GC,
              Geometry const *const geo,
              GParam const *const param,
              int mu,
              int pwr,
              double *repoly,
              double *impoly)
	{
	#ifdef DEBUG
	ASSERT(pwr > 0, "power of Polyakov loop must be at least 1, got %d", pwr);
	#endif

	double rep = 0.0, imp = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : rep) reduction(+ : imp)
	#endif
	for(long rsp = 0; rsp < param->d_space_vol[mu]; rsp++)
		{
		long r = sisp_and_mu_to_si(geo, rsp, 0, mu);
		GAUGE_GROUP matrix;
		one(&matrix);
		for(int i = 0; i < param->d_size[mu]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][mu]));
			r = nnp(geo, r, mu);
			}
		if(pwr > 1)
			{
			GAUGE_GROUP matrix2;

			equal(&matrix2, &matrix);
			for(int j = 1; j < pwr; j++) times_equal(&matrix, &matrix2);
			}

		rep += retr(&matrix);
		imp += imtr(&matrix);
		}

	*repoly = rep * (param->d_inv_space_vol[mu]);
	*impoly = imp * (param->d_inv_space_vol[mu]);
	}

// compute the mean trace of the Polyakov loop winding along multiple directions
void multipolyakov(Gauge_Conf const *const GC,
                   Geometry const *const geo,
                   GParam const *const param,
                   double *repoly,
                   double *impoly)
	{
	double rep = 0.0, imp = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : rep) reduction(+ : imp)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		GAUGE_GROUP matrix;
		one(&matrix);
		for(int j = 0; j < param->d_multipolyakov_order; j++)
			{
			int const mu = param->d_multipolyakov_dirs[j];
			for(int i = 0; i < param->d_size[mu]; i++)
				{
				times_equal(&matrix, &(GC->lattice[r][mu]));
				r = nnp(geo, r, mu);
				}
			}

		rep += retr(&matrix);
		imp += imtr(&matrix);
		}

	*repoly = rep * (param->d_inv_vol);
	*impoly = imp * (param->d_inv_vol);
	}


// compute the mean trace of the Polyakov loop in the adjoint representation:
// Tr_adj(P)/(N^2-1) = (|Tr_F(P)|^2 - 1)/(N^2 - 1) = (|tr_F(P)/N|^2 - 1/N^2)/(1 - 1/N^2)
void polyakov_adj(Gauge_Conf const *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  double *repoly,
                  double *impoly)
	{
	*impoly = 0.0;

	#if NCOLOR != 1

	double const inv_n2 = 1.0 / ((double) NCOLOR * NCOLOR);
	double const adj_tr_norm = 1.0 / (1.0 - inv_n2);
	double rep = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : rep)
	#endif
	for(long rsp = 0; rsp < param->d_space_vol[0]; rsp++)
		{
		long r = sisp_and_t_to_si(geo, rsp, 0);

		GAUGE_GROUP matrix;
		one(&matrix);
		for(int i = 0; i < param->d_size[0]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][0]));
			r = nnp(geo, r, 0);
			}
		double const tr_N_re = retr(&matrix);
		double const tr_N_im = imtr(&matrix);
		rep += (tr_N_re * tr_N_re + tr_N_im * tr_N_im - inv_n2);
		}

	*repoly = rep * adj_tr_norm * param->d_inv_space_vol[0];

	#else

	*repoly = 0.0;

	#endif
	}


// compute the mean trace of all the powers of the Polyakov loop winding along direction mu
void polyakov_powers(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     int mu,
                     double *repoly_pwrs,
                     double *impoly_pwrs)
	{
	double rep[MAX_POLY_PWR], imp[MAX_POLY_PWR];

	for(int i = 0; i < MAX_POLY_PWR; i++)
		{
		rep[i] = 0.0;
		imp[i] = 0.0;
		}

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : rep[:MAX_POLY_PWR]) reduction(+ : imp[:MAX_POLY_PWR])
	#endif
	for(long rsp = 0; rsp < param->d_space_vol[mu]; rsp++)
		{
		long r = sisp_and_mu_to_si(geo, rsp, 0, mu);
		GAUGE_GROUP matrix, matrix2;
		one(&matrix);
		for(int i = 0; i < param->d_size[mu]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][mu]));
			r = nnp(geo, r, mu);
			}
		equal(&matrix2, &matrix);
		for(int i = 0; i < MAX_POLY_PWR; i++)
			{
			rep[i] += retr(&matrix);
			imp[i] += imtr(&matrix);
			times_equal(&matrix, &matrix2);
			}
		}
	for(int i = 0; i < MAX_POLY_PWR; i++)
		{
		repoly_pwrs[i] = rep[i] * param->d_inv_space_vol[mu];
		impoly_pwrs[i] = imp[i] * param->d_inv_space_vol[mu];
		}
	}


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     long r)
	{
	#if (STDIM==4 && NCOLOR>1)

	// topcharge normalization: Nc from retr() = 1/N ReTr[], bulk factor Z if using OBC, 1 / (128 * pi ** 2) from definition
	double charge_norm = NCOLOR / (128.0 * PI * PI);
	#ifdef THETA_MODE
	if(param->d_meas_effective_charge)
		charge_norm *= creal(GC->Z[r][param->d_n_planes]);
	#endif

	int sign = -1;
	double loc_charge = 0.0;

	for(int i = 0; i < 3; i++)
		{
		int const mu = g_indep_perm_dir[0][i];
		int const nu = g_indep_perm_dir[1][i];
		int const rho = g_indep_perm_dir[2][i];
		int const sigma = g_indep_perm_dir[3][i];

		GAUGE_GROUP C_mu_nu, C_rho_sigma, C_tmp;

		clover(GC, geo, param, r, mu, nu, &C_mu_nu);
		clover(GC, geo, param, r, rho, sigma, &C_rho_sigma);

		times_dag2(&C_tmp, &C_rho_sigma, &C_mu_nu); // C_tmp = C_rho_sigma * C_mu_nu^{dag}
		double const real1 = retr(&C_tmp);

		times(&C_tmp, &C_rho_sigma, &C_mu_nu); // C_tmp = C_rho_sigma * C_mu_nu
		double const real2 = retr(&C_tmp);

		loc_charge += (double) sign * (real1 - real2);
		sign = -sign;
		}
	return loc_charge * charge_norm;

	#elif (STDIM==2 && NCOLOR==1)

	GAUGE_GROUP u1matrix;

	plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
	return atan2(cimag(u1matrix.comp), creal(u1matrix.comp)) / PI2;

	#else

	(void) GC;
	(void) geo;
	(void) param;
	(void) r;
	REQUIRE(0, "unsupported configuration for topological charge: STDIM=%d, NCOLOR=%d", STDIM, NCOLOR);
	return 0.0;

	#endif
	}


// compute and write to binary file the charge density
void charge_density(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    Meas_Utils *const meas_aux,
                    int const meas_count)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		meas_aux->scalar_density[r] = loc_topcharge(GC, geo, param, r);;
		//TODO: try parallel writing with pwrite() if available
		}

	FILE *fp = meas_aux->chargedensityfilep;
	fwrite(&(GC->update_index), sizeof(long), 1, fp);
	fwrite(&meas_count, sizeof(int), 1, fp);
	fwrite(meas_aux->scalar_density, sizeof(double), (size_t) param->d_volume, fp);
	}


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const *const GC,
                 Geometry const *const geo,
                 GParam const *const param)
	{
	double res = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : res)
	#endif
	for(long r = 0; r < (param->d_volume); r++) res += loc_topcharge(GC, geo, param, r);

	return res;
	}


// sum loc_topcharge over STDIM-1 dirs and then the abs value of the result over the remaining dir
double topcharge_prime(Gauge_Conf const *const GC,
                       Geometry const *const geo,
                       GParam const *const param,
                       int const dir)
	{
	double res = 0.0;

	double *tmp;
	allocate_array_double(&tmp, param->d_size[dir], __FILE__, __LINE__);
	for(int i = 0; i < param->d_size[dir]; i++) tmp[i] = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel num_threads(NTHREADS)
	#endif
	{
		double *tmp_private;
		allocate_array_double(&tmp_private, param->d_size[dir], __FILE__, __LINE__);
		for(int j = 0; j < param->d_size[dir]; j++) tmp_private[j] = 0;

		#ifdef OPENMP_MODE
		#pragma omp for
		#endif
		for(long r = 0; r < (param->d_volume); r++)
			{
			int cartcoord[STDIM];
			si_to_cart(cartcoord, r, param);
			tmp_private[cartcoord[dir]] += loc_topcharge(GC, geo, param, r);
			}
		#ifdef OPENMP_MODE
		#pragma omp critical
		#endif
		{
			for(int j = 0; j < param->d_size[dir]; j++) tmp[j] += tmp_private[j];
		}

		free(tmp_private);
	}

	for(int i = 0; i < param->d_size[dir]; i++) res += fabs(tmp[i]);
	free(tmp);

	return res;
	}


// chi^\prime = (1/8) int d^4x |x|^2 <q(x)q(0)> = < (1/8) int d^4x |x|^2 q(x) q(0) > = < G2 >
// This function computes the quantity (q(0)/8) sum_{x} d(x,0)^2 q(x) = a^2 G2, whose mean over the ensemble is a^2 chi^\prime
// d(x,y) = lattice distance between sites x and y keeping periodic boundary conditions into account (i.e., the shortest distance between x and y)
double topo_chi_prime(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param)
	{
	double res = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+: res)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		res += square_distance(r, 0, param) * loc_topcharge(GC, geo, param, r);
		}
	res *= loc_topcharge(GC, geo, param, 0) / (2 * STDIM); // res *= q(0) / 8

	return res;
	}


void topcharge_slices(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param,
                      int const mu,
                      double *slices,
                      int meas_count,
                      FILE *filep)
	{
	int const L_mu = param->d_size[mu];
	int const offset = (int) (GC->translation[mu] % L_mu);
	for(int i = 0; i < L_mu; i++) slices[i] = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:slices[:L_mu])
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		int const x_mu = periodic_condition(geo->d_mucomp[mu][r] - offset, L_mu);
		slices[x_mu] += loc_topcharge(GC, geo, param, r);
		}

	fprintf(filep, "%ld %d ", GC->update_index, meas_count);
	for(int i = 0; i < L_mu; i++) fprintf(filep, " % 18.12e", slices[i]);
	fprintf(filep, "\n");
	}


void topcharge_p_slices(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        int const mu,
                        int const nu,
                        Meas_Utils *const meas_aux,
                        int meas_count)
	{
	int const L_mu = param->d_size[mu];
	int const L_nu = param->d_size[nu];
	int const offset_mu = (int) (GC->translation[mu] % L_mu);
	int const offset_nu = (int) (GC->translation[nu] % L_nu);
	double *slicesre = meas_aux->real_slices;
	double *slicesim = meas_aux->imag_slices;
	FILE *const filep = meas_aux->q_slices_filep;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		meas_aux->scalar_density[r] = loc_topcharge(GC, geo, param, r);
	double const *const charge = meas_aux->scalar_density;

	for(int k = 0; k < L_nu; k++)
		{
		double const p = k * PI2 / L_nu;
		for(int i = 0; i < L_mu; i++)
			{
			slicesre[i] = 0.0;
			slicesim[i] = 0.0;
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) reduction(+:slicesre[:L_mu]) reduction(+:slicesim[:L_mu])
		#endif
		for(long r = 0; r < (param->d_volume); r++)
			{
			int const x_mu = periodic_condition(geo->d_mucomp[mu][r] - offset_mu, L_mu);
			int const x_nu = periodic_condition(geo->d_mucomp[nu][r] - offset_nu, L_nu);
			double const ph = p * x_nu;
			slicesre[x_mu] += cos(ph) * charge[r];
			slicesim[x_mu] += sin(ph) * charge[r];
			}

		fprintf(filep, "%ld %d %d ", GC->update_index, meas_count, k);
		for(int i = 0; i < L_mu; i++) fprintf(filep, " % 18.12e % 18.12e", slicesre[i], slicesim[i]);
		fprintf(filep, "\n");
		}
	}


// TODO: just to check how cooling destroys topological correlations, remove
void check_correlation_decay_cooling(Gauge_Conf const *const GC, Geometry const *const geo, GParam const *const param, double *ratio)
	{
	Gauge_Conf helperconf;
	init_gauge_conf_from_gauge_conf(&helperconf, GC, param);

	for(int i = 0; i < (param->d_coolrepeat); i++)
		{
		cooling(&helperconf, geo, param, param->d_coolsteps, LEX_DIR_LEXEO_SITE);
		double const Q = fabs(topcharge(&helperconf, geo, param));
		double satd = 0.0;

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) reduction(+ : satd)
		#endif
		for(long r = 0; r < (param->d_volume); r++)
			{
			satd += fabs(loc_topcharge(&helperconf, geo, param, r));
			}
		ratio[i] = 1.0 - Q / satd;
		}
	free_gauge_conf(&helperconf, param);
	}


// compute the correlator of the local topological charge
// after "ncool" cooling steps up to spatial distance "dist"
void loc_topcharge_corr(Gauge_Conf *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        int ncool,
                        int dist,
                        double *res)
	{
	#if STDIM < 2
	for(int i = 0; i < dist; i++)
		res[i] = 0.0;
	#else
	// TODO: refactor with meas_aux
	double *charge_array;
	allocate_array_double(&charge_array, param->d_volume, __FILE__, __LINE__);

	// compute the local topological charge
	if(ncool > 0) cooling(GC, geo, param, ncool, LEX_DIR_LEXEO_SITE);

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++) charge_array[r] = loc_topcharge(GC, geo, param, r);

	if(ncool > 0) restore_gauge_conf(GC, param);

	// compute correlators
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int i = 0; i < dist; i++)
		{
		double corr = 0.0;
		for(long r1 = 0; r1 < param->d_volume; r1++)
			{
			double aux = 0.0;
			for(int dir = 1; dir < STDIM; dir++)
				{
				long r2 = r1;
				for(int j = 0; j < i; j++) r2 = nnp(geo, r2, dir);
				aux += charge_array[r2];
				}
			corr += aux * charge_array[r1];
			}
		res[i] = corr * param->d_inv_vol / (double) (STDIM - 1);
		}

	// free memory
	free(charge_array);
	#endif
	}


void perform_measures_aux(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                          int const meas_count, Meas_Utils *meas_aux)
	{
	if(param->d_plaquette_meas == 1)
		{
		double plaqs, plaqt;
		plaquette(GC, geo, param, &plaqs, &plaqt);
		meas_aux->meanplaq[meas_count] = ((STDIM - 2.0) * plaqs + 2.0 * plaqt) / STDIM;
		}
	if(param->d_clover_energy_meas == 1) clover_disc_energy(GC, geo, param, &(meas_aux->clover_energy[meas_count]));
	if(param->d_energy_density_meas == 1) energy_density(GC, geo, param, meas_aux, meas_count + 1);
	if(param->d_charge_meas == 1) meas_aux->charge[meas_count] = topcharge(GC, geo, param);
	if(param->d_charge_density_meas == 1) charge_density(GC, geo, param, meas_aux, meas_count + 1);
	if(param->d_polyakov_meas == 1) for(int i = 0; i < STDIM; i++) polyakov(GC, geo, param, i, 1, &(meas_aux->polyre[meas_count][i]), &(meas_aux->polyim[meas_count][i]));
	if(param->d_multipolyakov_order >= 1) multipolyakov(GC, geo, param, &(meas_aux->multipolyre[meas_count]), &(meas_aux->multipolyim[meas_count]));
	if(param->d_polyakov_density_meas == 1) for(int i = 0; i < STDIM; i++) polyakov_density(GC, geo, param, i, meas_aux, meas_count + 1);
	if(param->d_energy_slices_meas == 1) clover_energy_slices(GC, geo, param, 0, meas_aux->real_slices, meas_count + 1, meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) topcharge_slices(GC, geo, param, 0, meas_aux->real_slices, meas_count + 1, meas_aux->q_slices_filep);
	if(param->d_charge_p_slices_meas == 1) topcharge_p_slices(GC, geo, param, 0, param->d_test_flag, meas_aux, meas_count + 1);
	if(param->d_chi_prime_meas == 1) meas_aux->chi_prime[meas_count] = topo_chi_prime(GC, geo, param);
	if(param->d_charge_prime_meas == 1) for(int i = 0; i < STDIM; i++) meas_aux->charge_prime[meas_count][i] = topcharge_prime(GC, geo, param, i);
	if(param->d_action_meas == 1) action(GC, geo, param, &(meas_aux->action1[meas_count]), &(meas_aux->action2[meas_count]), &(meas_aux->action3[meas_count]), &(meas_aux->potential[meas_count]));
	}


void perform_measures_localobs_hot(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                                   Meas_Utils *meas_aux)
	{
	int i;
	double plaqs = 0.0, plaqt = 0.0, clover_energy = 0.0, charge = 0.0, chi_prime = 0.0, multipolyre = 0.0, multipolyim = 0.0; // =0.0 to suppress gcc warning
	double polyre[STDIM], polyim[STDIM], polyre_pwrs[MAX_POLY_PWR], polyim_pwrs[MAX_POLY_PWR], charge_prime[STDIM];
	double action1 = 0.0, action2 = 0.0, action3 = 0.0, potential = 0.0;

	// perform meas
	if(param->d_plaquette_meas == 1) plaquette(GC, geo, param, &plaqs, &plaqt);
	if(param->d_clover_energy_meas == 1) clover_disc_energy(GC, geo, param, &clover_energy);
	if(param->d_energy_density_meas == 1) energy_density(GC, geo, param, meas_aux, 0);
	if(param->d_charge_meas == 1) charge = topcharge(GC, geo, param);
	if(param->d_charge_density_meas == 1) charge_density(GC, geo, param, meas_aux, 0);
	if(param->d_polyakov_meas == 1) for(i = 0; i < STDIM; i++) polyakov(GC, geo, param, i, 1, &(polyre[i]), &(polyim[i]));
	if(param->d_polyakov_powers_meas == 1) polyakov_powers(GC, geo, param, 0, polyre_pwrs, polyim_pwrs);
	if(param->d_multipolyakov_order >= 1) multipolyakov(GC, geo, param, &(multipolyre), &(multipolyim));
	if(param->d_polyakov_density_meas == 1) for(i = 0; i < STDIM; i++) polyakov_density(GC, geo, param, i, meas_aux, 0);
	if(param->d_chi_prime_meas == 1) chi_prime = topo_chi_prime(GC, geo, param);
	if(param->d_charge_prime_meas == 1) for(i = 0; i < STDIM; i++) charge_prime[i] = topcharge_prime(GC, geo, param, i);
	if(param->d_energy_slices_meas == 1) clover_energy_slices(GC, geo, param, 0, meas_aux->real_slices, 0, meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) topcharge_slices(GC, geo, param, 0, meas_aux->real_slices, 0, meas_aux->q_slices_filep);
	if(param->d_charge_p_slices_meas == 1) topcharge_p_slices(GC, geo, param, 0, param->d_test_flag, meas_aux, 0);
	if(param->d_action_meas == 1) action(GC, geo, param, &action1, &action2, &action3, &potential);

	// print meas (density profiles already printed)
	fprintf(meas_aux->datafilep, "%ld ", GC->update_index);
	if(param->d_plaquette_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", plaqs, plaqt);
	if(param->d_clover_energy_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", clover_energy);
	if(param->d_charge_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", charge);
	if(param->d_polyakov_meas == 1) for(i = 0; i < STDIM; i++) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", polyre[i], polyim[i]);
	if(param->d_polyakov_powers_meas == 1) for(i = 0; i < MAX_POLY_PWR; i++) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", polyre_pwrs[i], polyim_pwrs[i]);
	if(param->d_multipolyakov_order >= 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", multipolyre, multipolyim);
	if(param->d_chi_prime_meas == 1) fprintf(meas_aux->chiprimefilep, "%ld 0 % 18.12e\n", GC->update_index, chi_prime);
	if(param->d_charge_prime_meas == 1) for(i = 0; i < STDIM; i++) fprintf(meas_aux->datafilep, "% 18.12e ", charge_prime[i]);
	if(param->d_action_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e % 18.12e % 18.12e ", action1, action2, action3, potential);

	// flush data files
	fflush(meas_aux->datafilep);
	if(param->d_energy_density_meas == 1) fflush(meas_aux->energydensityfilep);
	if(param->d_charge_density_meas == 1) fflush(meas_aux->chargedensityfilep);
	if(param->d_polyakov_density_meas == 1) for(i = 0; i < STDIM; i++) fflush(meas_aux->polyakovdensityfilep[i]);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


void perform_measures_localobs(Gauge_Conf *const GC,
                               Geometry const *const geo,
                               GParam const *const param,
                               Meas_Utils *meas_aux)
	{
	// measures without smoothing
	perform_measures_localobs_hot(GC, geo, param, meas_aux);

	// measures with adaptive gradient flow
	if(param->d_agf_num_meas > 0) perform_measures_localobs_adaptive_gradflow(GC, geo, param, meas_aux);

	// measures with fixed-step gradient flow
	if(param->d_gf_num_meas > 0) perform_measures_localobs_gradflow(GC, geo, param, meas_aux);

	// measures with cooling
	if(param->d_coolrepeat > 0) perform_measures_localobs_cooling(GC, geo, param, meas_aux, param->d_cooling_type);

	// multicanonical topcharge and weight
	#ifdef MULTICANONICAL_MODE
	double x = GC->stored_topcharge;
	double V = compute_topo_potential(GC->replica_index, x, param);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", x, exp(V));
	#endif

	// newline and flush data files
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


void perform_measures_localobs_cooling(Gauge_Conf *const GC,
                                       Geometry const *const geo,
                                       GParam const *const param,
                                       Meas_Utils *meas_aux, Cooling_Type const type)
	{
	for(int meas_count = 0; meas_count < param->d_coolrepeat; meas_count++)
		{
		cooling(GC, geo, param, param->d_coolsteps, type);
		perform_measures_aux(GC, geo, param, meas_count, meas_aux);
		}

	// restore gauge conf before cooling from GC->lattice_copy
	restore_gauge_conf(GC, param);

	// print meas cooling
	print_measures_aux(param->d_coolrepeat, GC->update_index, param, meas_aux);

	//TODO: remove, debug only
	//for(int meas_count=0; meas_count < param->d_coolrepeat; meas_count++)
	//	{
	//	cooling_with_defect(GC, geo, param, param->d_coolsteps);
	//	perform_measures_aux(GC, geo, param, meas_count, meas_aux);
	//	}
	//restore_gauge_conf(GC, param);
	//print_measures_aux(param->d_coolrepeat, GC->update_index, param, meas_aux);
	}


void perform_measures_localobs_gradflow(Gauge_Conf *const GC,
                                        Geometry const *const geo,
                                        GParam const *const param,
                                        Meas_Utils *meas_aux)
	{
	// count starts from 1 to avoid problems with %
	for(int count = 1; count < (param->d_ngfsteps + 1); count++)
		{
		gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
		if(count % param->d_gf_meas_each == 0)
			{
			int const meas_count = count / param->d_gf_meas_each - 1;
			perform_measures_aux(GC, geo, param, meas_count, meas_aux);
			}
		}

	// restore gauge conf before gradflow from GC->lattice_copy
	restore_gauge_conf(GC, param);

	// print meas gradflow
	print_measures_aux(param->d_gf_num_meas, GC->update_index, param, meas_aux);
	}


void perform_measures_localobs_adaptive_gradflow(Gauge_Conf *const GC,
                                                 Geometry const *const geo,
                                                 GParam const *const param,
                                                 Meas_Utils *meas_aux)
	{
	int accepted;
	double gftime = 0.0;
	double gftime_step = param->d_agf_step;
	int meas_count = 0;

	while(meas_count < param->d_agf_num_meas)
		{
		gradflow_RKstep_adaptive(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux);

		// if step accepted, perform measures
		if(accepted == 1 && fabs(gftime - param->d_agf_meas_each * (meas_count + 1)) - param->d_agf_time_bin < MIN_VALUE)
			{
			perform_measures_aux(GC, geo, param, meas_count, meas_aux);
			meas_count = meas_count + 1;
			}

		// adapt step to the time of next measure if this would be skipped
		if((gftime + gftime_step - param->d_agf_meas_each * (meas_count + 1)) > param->d_agf_time_bin)
			{
			gftime_step = param->d_agf_meas_each * (meas_count + 1) - gftime;
			}
		}

	// restore gauge conf before gradflow from GC->lattice_copy
	restore_gauge_conf(GC, param);

	// print meas gradflow
	print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);
	}


// Compare different definitions of the local integration error
void perform_measures_localobs_with_adaptive_gradflow_debug1(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux)
	{
	// meas no gradflow
	perform_measures_localobs_hot(GC, geo, param, meas_aux);

	// meas gradflow
	if(param->d_agf_num_meas > 0)
		{
		// allocate memory
		Gauge_Conf helper1, helper2;
		init_gauge_conf_from_gauge_conf(&helper1, GC, param);
		init_gauge_conf_from_gauge_conf(&helper2, GC, param);

		// gradflow starts
		int accepted;
		double gftime = 0.0;
		double gftime_step = param->d_agf_step;
		int meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			gradflow_RKstep_adaptive_check(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux, &helper1, &helper2);

			// if step accepted, perform measures
			if(accepted == 1 && fabs(gftime - param->d_agf_meas_each * (meas_count + 1)) - param->d_agf_time_bin < MIN_VALUE)
				{
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}

			// adapt step to the time of next if this would be skipped
			if((gftime + gftime_step - param->d_agf_meas_each * (meas_count + 1)) > param->d_agf_time_bin)
				{
				gftime_step = param->d_agf_meas_each * (meas_count + 1) - gftime;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);

		// free memory
		free_gauge_conf(&helper1, param);
		free_gauge_conf(&helper2, param);
		}

	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	double x = GC->stored_topcharge;
	double V = compute_topo_potential(GC->replica_index, x, param);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", x, exp(V));
	#endif

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


// Test the local error of the adaptive integration along a fixed-step integration
void perform_measures_localobs_with_adaptive_gradflow_debug2(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux)
	{
	// meas no gradflow
	perform_measures_localobs_hot(GC, geo, param, meas_aux);

	// meas gradflow
	if(param->d_agf_num_meas > 0)
		{
		// allocate memory
		Gauge_Conf helper1, helper2;
		init_gauge_conf_from_gauge_conf(&helper1, GC, param);
		init_gauge_conf_from_gauge_conf(&helper2, GC, param);

		// gradflow starts
		double gftime = 0.0;
		double gftime_step_max = 0.40;
		double gftime_step_min = 0.01;
		int accepted;
		int meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			double gftime_step = gftime_step_max;
			while(gftime_step > gftime_step_min)
				{
				gradflow_RKstep_adaptive_check(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux, &helper1, &helper2);
				equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) helper1.lattice, param);
				gftime_step *= 0.8;
				}
			gradflow_RKstep_adaptive_check(GC, geo, param, &gftime, &gftime_step_min, &accepted, meas_aux, &helper1, &helper2);
			gftime += gftime_step_min;

			// perform measures
			if(gftime - param->d_agf_meas_each * (meas_count + 1) > -gftime_step_min / 2)
				{
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}

			if(gftime_step_min < gftime / 50.0)
				{
				gftime_step_min = gftime / 50.0;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);

		// free memory
		free_gauge_conf(&helper1, param);
		free_gauge_conf(&helper2, param);
		}

	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	double x = GC->stored_topcharge;
	double V = compute_topo_potential(GC->replica_index, x, param);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", x, exp(V));
	#endif

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


// Test the local error of the adaptive integration along an adaptive-step integration
void perform_measures_localobs_with_adaptive_gradflow_debug3(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux)
	{
	// meas no gradflow
	perform_measures_localobs_hot(GC, geo, param, meas_aux);

	// meas gradflow
	if(param->d_agf_num_meas > 0)
		{
		// allocate memory
		Gauge_Conf helper1, helper2, accepted_conf;
		init_gauge_conf_from_gauge_conf(&helper1, GC, param);
		init_gauge_conf_from_gauge_conf(&helper2, GC, param);
		init_gauge_conf_from_gauge_conf(&accepted_conf, GC, param);

		// gradflow starts
		double gftime = 0.0;
		double gftime_step_max = 0.4;
		double gftime_step_min = 0.01;
		double gftime_step_accepted = 0.0;
		int accepted;
		int meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			int accepted_flag = 0;
			double gftime_step = gftime_step_max;
			while(gftime_step > gftime_step_min)
				{
				gradflow_RKstep_adaptive_check(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux, &helper1, &helper2);
				if(accepted == 1 && accepted_flag == 0)
					{
					equal_lattice(accepted_conf.lattice, (GAUGE_GROUP const * const *) GC->lattice, param);
					gftime_step_accepted = gftime_step;
					accepted_flag = 1;
					}
				equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) helper1.lattice, param);
				gftime_step *= 0.8;
				}
			gradflow_RKstep_adaptive_check(GC, geo, param, &gftime, &gftime_step_min, &accepted, meas_aux, &helper1, &helper2);
			if(accepted_flag == 1)
				{
				equal_lattice(GC->lattice, (GAUGE_GROUP const * const *) accepted_conf.lattice, param);
				}
			else
				{
				gftime_step_accepted = gftime_step_min;
				}
			gftime += gftime_step_accepted;

			// perform measures
			if(gftime - param->d_agf_meas_each * (meas_count + 1) > param->d_agf_meas_each / 10000.)
				{
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}

			if(gftime_step_min < gftime / 50.0)
				{
				gftime_step_min = gftime / 50.0;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);

		// free memory
		free_gauge_conf(&helper1, param);
		free_gauge_conf(&helper2, param);
		}

	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	double x = GC->stored_topcharge;
	double V = compute_topo_potential(GC->replica_index, x, param);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", x, exp(V));
	#endif

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


// Compare the total error of the adaptive integration with a fixed-step integration
void perform_measures_localobs_with_adaptive_gradflow_debug4(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux)
	{
	// meas no gradflow
	perform_measures_localobs_hot(GC, geo, param, meas_aux);

	// meas gradflow
	if(param->d_agf_num_meas > 0)
		{
		// allocate memory
		GAUGE_GROUP ***helper_confs;
		allocate_array_GAUGE_GROUP_pointer_pointer(&helper_confs, param->d_agf_num_meas, __FILE__, __LINE__);
		for(int i = 0; i < param->d_agf_num_meas; i++)
			{
			allocate_array_GAUGE_GROUP_pointer(&(helper_confs[i]), param->d_volume, __FILE__, __LINE__);
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < (param->d_volume); r++)
				{
				allocate_array_GAUGE_GROUP(&(helper_confs[i][r]), STDIM, __FILE__, __LINE__);
				}
			}

		// adaptive step integration
		double gftime_step = param->d_agf_step;
		double gftime = 0.0;
		int accepted = 0;
		int meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			gradflow_RKstep_adaptive(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux);

			// if step accepted, save conf
			if(accepted == 1 && fabs(gftime - param->d_agf_meas_each * (meas_count + 1)) - param->d_agf_time_bin < MIN_VALUE)
				{
				equal_lattice(helper_confs[meas_count], (GAUGE_GROUP const * const *) GC->lattice, param);
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}

			// adapt step to the time of next measure if this would be skipped
			if((gftime + gftime_step - param->d_agf_meas_each * (meas_count + 1)) > param->d_agf_time_bin)
				{
				gftime_step = param->d_agf_meas_each * (meas_count + 1) - gftime;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);

		// fixed step integration
		gftime_step = 0.01;
		gftime = 0.0;
		meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			gradflow_RKstep(GC, geo, param, gftime_step, meas_aux);
			gftime += gftime_step;
			gftime_step = 0.01;

			// save conf
			if(fabs(gftime - param->d_agf_meas_each * (meas_count + 1)) - param->d_agf_time_bin < MIN_VALUE)
				{
				double max_dist = lattice_max_dist((GAUGE_GROUP const * const *) GC->lattice, (GAUGE_GROUP const * const *) helper_confs[meas_count], param);
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;

				fprintf(stderr, "%ld %.12g %.12g\n", GC->update_index, gftime, max_dist);
				fflush(stderr);
				}

			// adapt step to the time of next measure if this would be skipped
			if((gftime + gftime_step - param->d_agf_meas_each * (meas_count + 1)) > param->d_agf_time_bin)
				{
				gftime_step = param->d_agf_meas_each * (meas_count + 1) - gftime;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_aux(param->d_agf_num_meas, GC->update_index, param, meas_aux);

		// free memory
		for(int i = 0; i < param->d_agf_num_meas; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < (param->d_volume); r++)
				{
				free(helper_confs[i][r]);
				}
			free(helper_confs[i]);
			}
		free(helper_confs);
		}

	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	double x = GC->stored_topcharge;
	double V = compute_topo_potential(GC->replica_index, x, param);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", x, exp(V));
	#endif

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if(param->d_energy_slices_meas == 1) fflush(meas_aux->e_slices_filep);
	if(param->d_charge_slices_meas == 1) fflush(meas_aux->q_slices_filep);
	if(param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	}


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                FILE *datafilep)
	{
	const int max_hit = 50;
	const int dir = 1;

	double complex *poly_array;
	allocate_array_double_complex(&poly_array, param->d_space_vol[0], __FILE__, __LINE__);

	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0);
	#endif

	fprintf(datafilep, "Multihit optimization: \n");
	fprintf(datafilep, "the smaller the value the better the multihit\n");

	for(int mh = 1; mh < max_hit; mh++)
		{
		time_t time1, time2;
		time(&time1);

		// polyakov loop computation
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			GAUGE_GROUP matrix, tmp;
			one(&matrix);
			for(int i = 0; i < param->d_size[0]; i++)
				{
				multihit(GC, geo, param, sisp_and_t_to_si(geo, r, i), 0, mh, &tmp);
				times_equal(&matrix, &tmp);
				}
			poly_array[r] = retr(&matrix) + I * imtr(&matrix);
			}

		// average correlator computation
		double complex poly_corr = 0.0 + I * 0.0;
		double poly_corr_abs = 0.0;
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			long r1 = sisp_and_t_to_si(geo, r, 0);
			for(int i = 0; i < param->d_dist_poly; i++)
				{
				r1 = nnp(geo, r1, dir);
				}
			long r2;
			int t_tmp;
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1, t_tmp is unused

			poly_corr += poly_array[r] * conj(poly_array[r2]);
			poly_corr_abs += cabs(poly_array[r] * conj(poly_array[r2]));
			}
		poly_corr *= param->d_inv_space_vol[0];
		poly_corr_abs *= param->d_inv_space_vol[0] * sqrt(mh);

		// fluctuation of the average correlator computation
		double poly_corr_fluct = 0.0;
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			long r1 = sisp_and_t_to_si(geo, r, 0);
			for(int i = 0; i < param->d_dist_poly; i++)
				{
				r1 = nnp(geo, r1, dir);
				}
			long r2;
			int t_tmp;
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1, t_tmp is unused
			poly_corr_fluct += cabs(poly_array[r] * conj(poly_array[r2]) - poly_corr);
			}
		poly_corr_fluct *= param->d_inv_space_vol[0] * sqrt(mh);

		time(&time2);
		double const diff_sec = difftime(time2, time1);
		fprintf(datafilep, "%d	% 18.12e	% 18.12e (time:%g)\n", mh, poly_corr_abs, poly_corr_fluct, diff_sec);

		fflush(datafilep);
		}

	free(poly_array);
	}


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  FILE *datafilep)
	{
	double complex *poly_array;
	allocate_array_double_complex(&poly_array, param->d_space_vol[0], __FILE__, __LINE__);

	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	multilevel_polycorr(GC, geo, param, param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	// averages
	double complex poly_corr = 0.0 + I * 0.0;
	double poly_corr_abs = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		poly_array[r] = retr_TensProd(&(GC->ml_polycorr[0][0][r])) + I * imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

		poly_corr += poly_array[r];
		poly_corr_abs += cabs(poly_array[r]);
		}
	poly_corr *= param->d_inv_space_vol[0];
	poly_corr_abs *= param->d_inv_space_vol[0];

	// fluctuations
	// double poly_corr_fluct = 0.0;
	// for(long r = 0; r < param->d_space_vol[0]; r++)
	// 	{
	// 	poly_corr_fluct += cabs(poly_array[r] - poly_corr);
	// 	}
	// poly_corr_fluct *= param->d_inv_space_vol[0];

	// normalizations
	for(int i = 0; i < NLEVELS; i++)
		{
		poly_corr_abs *= sqrt(param->d_ml_upd[i]);
		// poly_corr_fluct *= sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs *= sqrt(param->d_multihit);
	// poly_corr_fluct *= sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e ", poly_corr_abs);
	for(int i = 0; i < NLEVELS; i++)
		{
		fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
		}
	fprintf(datafilep, "(1, %d) ", param->d_multihit);
	fprintf(datafilep, "\n");
	fflush(datafilep);

	free(poly_array);
	}


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf *const GC,
                               Geometry const *const geo,
                               GParam const *const param,
                               Meas_Utils *meas_aux)
	{
	#ifndef OPT_MULTIHIT
	#ifndef OPT_MULTILEVEL
	multilevel_polycorr(GC, geo, param, param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	double res = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res += retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	res *= param->d_inv_space_vol[0];

	fprintf(meas_aux->datafilep, "% 18.12e\n", res);
	fflush(meas_aux->datafilep);
	#endif
	#endif

	#ifdef OPT_MULTIHIT
	optimize_multihit_polycorr(GC, geo, param, meas_aux->datafilep);
	#endif

	#ifdef OPT_MULTILEVEL
	optimize_multilevel_polycorr(GC, geo, param, meas_aux->datafilep);
	#endif
	}


// to optimize the number of hits to be used in multilevel for the adjoint representation
void optimize_multihit_polycorradj(Gauge_Conf *const GC,
                                   Geometry const *const geo,
                                   GParam const *const param,
                                   FILE *datafilep)
	{
	int const max_hit = 50;
	int const dir = 1;
	double const inv_n2 = 1.0 / ((double) NCOLOR * NCOLOR);

	double *poly_array;
	allocate_array_double(&poly_array, param->d_space_vol[0], __FILE__, __LINE__);

	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0);
	#endif

	fprintf(datafilep, "Multihit optimization: \n");
	fprintf(datafilep, "the smaller the value the better the multihit\n");

	for(int mh = 1; mh < max_hit; mh++)
		{
		time_t time1, time2;
		time(&time1);

		// polyakov loop in the adjoint representation computation
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			#if NCOLOR == 1
			ASSERT(0, "adjoint representation is trivial for NCOLOR == 1 ");
			poly_array[r] = 0.0;
			#else
			GAUGE_GROUP matrix, tmp;
			one(&matrix);
			for(int i = 0; i < param->d_size[0]; i++)
				{
				multihit(GC, geo, param, sisp_and_t_to_si(geo, r, i), 0, mh, &tmp);
				times_equal(&matrix, &tmp);
				}

			//trace of the matrix in the fundamental representation
			double const tr_N_re = retr(&matrix);
			double const tr_N_im = imtr(&matrix);
			double const tr2_N2 = tr_N_re * tr_N_re + tr_N_im * tr_N_im;

			//trace of the matrix in adjoint representation
			poly_array[r] = (tr2_N2 - inv_n2) / (1.0 - inv_n2);
			#endif
			}

		// average correlator computation
		double poly_corr = 0.0;
		double poly_corr_abs = 0.0;
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			long r1 = sisp_and_t_to_si(geo, r, 0);
			for(int i = 0; i < param->d_dist_poly; i++)
				{
				r1 = nnp(geo, r1, dir);
				}
			int t_tmp;
			long r2;
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1, t_tmp is unused

			poly_corr += poly_array[r] * poly_array[r2];
			poly_corr_abs += fabs(poly_array[r] * poly_array[r2]);
			}
		poly_corr *= param->d_inv_space_vol[0];
		poly_corr_abs *= param->d_inv_space_vol[0];

		// fluctuation of the average correlator computation
		double poly_corr_fluct = 0.0;
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			long r1 = sisp_and_t_to_si(geo, r, 0);
			for(int i = 0; i < param->d_dist_poly; i++)
				{
				r1 = nnp(geo, r1, dir);
				}
			int t_tmp;
			long r2;
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1, t_tmp is unused

			poly_corr_fluct += fabs(poly_array[r] * poly_array[r2] - poly_corr);
			}
		poly_corr_fluct *= param->d_inv_space_vol[0];

		time(&time2);
		double const diff_sec = difftime(time2, time1);
		fprintf(datafilep, "%d	% 18.12e	% 18.12e (time:%g)\n", mh, poly_corr_abs * sqrt(mh), poly_corr_fluct * sqrt(mh), diff_sec);
		fflush(datafilep);
		}

	free(poly_array);
	}


// to optimize the multilevel (adjoint representation)
void optimize_multilevel_polycorradj(Gauge_Conf *const GC,
                                     Geometry const *const geo,
                                     GParam const *const param,
                                     FILE *datafilep)
	{
	double *poly_array;
	allocate_array_double(&poly_array, param->d_space_vol[0], __FILE__, __LINE__);

	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	multilevel_polycorradj(GC,
	                       geo,
	                       param,
	                       param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]));
			}
		}

	// averages
	double poly_corr = 0.0;
	double poly_corr_abs = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		poly_array[r] = retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));

		poly_corr += poly_array[r];
		poly_corr_abs += fabs(poly_array[r]);
		}
	poly_corr *= param->d_inv_space_vol[0];
	poly_corr_abs *= param->d_inv_space_vol[0];

	// fluctuations
	double poly_corr_fluct = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		poly_corr_fluct += fabs(poly_array[r] - poly_corr);
		}
	poly_corr_fluct *= param->d_inv_space_vol[0];

	// normalizations
	for(int i = 0; i < NLEVELS; i++)
		{
		poly_corr_abs *= sqrt(param->d_ml_upd[i]);
		poly_corr_fluct *= sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs *= sqrt(param->d_multihit);
	poly_corr_fluct *= sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e % 18.12e ", poly_corr_abs, poly_corr_fluct);
	for(int i = 0; i < NLEVELS; i++)
		{
		fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
		}
	fprintf(datafilep, "(1, %d) \n", param->d_multihit);
	fflush(datafilep);

	free(poly_array);
	}


// perform the computation of the polyakov loop correlator in the adjoint representation with the multilevel algorithm
void perform_measures_polycorradj(Gauge_Conf *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  Meas_Utils *meas_aux)
	{
	#ifndef OPT_MULTIHIT
	#ifndef OPT_MULTILEVEL

	multilevel_polycorradj(GC, geo, param, param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]));
			}
		}

	double res = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res += retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));
		}
	res *= param->d_inv_space_vol[0];

	fprintf(meas_aux->datafilep, "% 18.12e\n", res);
	fflush(meas_aux->datafilep);
	#endif
	#endif

	#ifdef OPT_MULTIHIT
	optimize_multihit_polycorradj(GC, geo, param, meas_aux->datafilep);
	#endif

	#ifdef OPT_MULTILEVEL
	optimize_multilevel_polycorradj(GC, geo, param, meas_aux->datafilep);
	#endif
	}


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *const GC,
                                       GParam const *const param,
                                       FILE *datafilep)
	{
	double complex *poly_array;
	allocate_array_double_complex(&poly_array, param->d_space_vol[0], __FILE__, __LINE__);

	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	// average
	double complex poly_corr = 0.0 + I * 0.0;
	double poly_corr_abs = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		poly_array[r] = retr_TensProd(&(GC->ml_polycorr[0][0][r])) + I * imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

		poly_corr += poly_array[r];
		poly_corr_abs += cabs(poly_array[r]);
		}
	poly_corr *= param->d_inv_space_vol[0];
	poly_corr_abs *= param->d_inv_space_vol[0];

	// fluctuation
	double poly_corr_fluct = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		poly_corr_fluct += cabs(poly_array[r] - poly_corr);
		}
	poly_corr_fluct *= param->d_inv_space_vol[0];

	// normalization
	for(int i = 0; i < NLEVELS; i++)
		{
		poly_corr_abs *= sqrt(param->d_ml_upd[i]);
		poly_corr_fluct *= sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs *= sqrt(param->d_ml_level0_repeat);
	poly_corr_fluct *= sqrt(param->d_ml_level0_repeat);

	poly_corr_abs *= sqrt(param->d_multihit);
	poly_corr_fluct *= sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e % 18.12e ", poly_corr_abs, poly_corr_fluct);
	for(int i = 0; i < NLEVELS; i++)
		{
		fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
		}
	fprintf(datafilep, "(1, %d) ", param->d_multihit);
	fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
	fprintf(datafilep, "\n");

	fflush(datafilep);

	free(poly_array);
	}


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *const GC,
                                    GParam const *const param,
                                    Meas_Utils *meas_aux)
	{
	#ifdef OPT_MULTILEVEL
	optimize_multilevel_polycorr_long(GC, param, meas_aux->datafilep);
	#else
	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	double res = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res += retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	res *= param->d_inv_space_vol[0];

	fprintf(meas_aux->datafilep, "% 18.12e\n", res);
	fflush(meas_aux->datafilep);
	#endif
	}


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                Meas_Utils *meas_aux)
	{
	multilevel_tube_disc(GC, geo, param, param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	double res_re = 0.0;
	double res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		res_im += imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	res_re = 0.0;
	res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polyplaq[0][r]));
		res_im += imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                Meas_Utils *meas_aux)
	{
	multilevel_tube_conn(GC, geo, param, param->d_size[0]);

	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]));
			times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	double res_re = 0.0;
	double res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		res_im += imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	res_re = 0.0;
	res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polyplaq[0][r]));
		res_im += imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	res_re = 0.0;
	res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		res_im += imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


// print the value of the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *const GC,
                                     GParam const *const param,
                                     Meas_Utils *meas_aux)
	{
	for(int i = 1; i < param->d_size[0] / param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS)
		#endif
		for(long r = 0; r < param->d_space_vol[0]; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]));
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]));
			times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]));
			}
		}

	double res_re = 0.0;
	double res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		res_im += imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	res_re = 0.0;
	res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polyplaq[0][r]));
		res_im += imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	res_re = 0.0;
	res_im = 0.0;
	for(long r = 0; r < param->d_space_vol[0]; r++)
		{
		res_re += retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		res_im += imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		}
	res_re *= param->d_inv_space_vol[0];
	res_im *= param->d_inv_space_vol[0];
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", res_re, res_im);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


int sprintf_header_datafile_aux(char *const header, char *const smoothing_method, GParam const *const param)
	{
	int j = sprintf(header, "( ");
	if(param->d_plaquette_meas == 1) j += sprintf(header + j, "plaq ");
	if(param->d_clover_energy_meas == 1) j += sprintf(header + j, "clover_energy ");
	if(param->d_charge_meas == 1) j += sprintf(header + j, "charge ");
	if(param->d_polyakov_meas == 1) for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "polyre_%d polyim_%d ", i, i);
	if(param->d_multipolyakov_order >= 1) j += sprintf(header + j, "multipolyre multipolyim ");
	if(param->d_charge_prime_meas == 1) for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "charge_prime_%d ", i);
	if(param->d_action_meas == 1)
		{
		for(int i = 1; i < 4; i++) j += sprintf(header + j, "action_%d ", i);
		j += sprintf(header + j, "potential ");
		}
	j += sprintf(header + j, ") x %s ", smoothing_method);
	return j;
	}


void header_datafile(char *const header, GParam const *const param)
	{
	char smoothing_method[STD_STRING_LENGTH];

	int j = sprintf(header, "# %d ", STDIM);
	for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "%d ", param->d_size[i]);
	j += sprintf(header + j, "\n# upd_index ");
	if(param->d_plaquette_meas == 1) j += sprintf(header + j, "plaqs plaqt ");
	if(param->d_clover_energy_meas == 1) j += sprintf(header + j, "clover_energy ");
	if(param->d_charge_meas == 1) j += sprintf(header + j, "charge ");
	if(param->d_polyakov_meas == 1) for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "polyre_%d polyim_%d ", i, i);
	if(param->d_polyakov_powers_meas == 1) for(int i = 0; i < MAX_POLY_PWR; i++) j += sprintf(header + j, "polyre_%d^%d polyim_%d^%d ", 0, i + 1, 0, i + 1);
	if(param->d_multipolyakov_order >= 1) j += sprintf(header + j, "multipolyre multipolyim ");
	if(param->d_charge_prime_meas == 1) for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "charge_prime_%d ", i);
	if(param->d_action_meas == 1)
		{
		for(int i = 1; i < 4; i++) j += sprintf(header + j, "action_%d ", i);
		j += sprintf(header + j, "potential ");
		}

	if(param->d_agf_num_meas > 0)
		{
		sprintf(smoothing_method, "%d agfrepeat each dt = %.10lf", param->d_agf_num_meas, param->d_agf_meas_each);
		j += sprintf_header_datafile_aux(header + j, smoothing_method, param);
		}
	if(param->d_gf_num_meas > 0)
		{
		sprintf(smoothing_method, "%d gfrepeat each dt = %.10lf", param->d_gf_num_meas, param->d_gf_meas_each * param->d_gfstep);
		j += sprintf_header_datafile_aux(header + j, smoothing_method, param);
		}
	if(param->d_coolrepeat > 0)
		{
		sprintf(smoothing_method, "%d coolrepeat each ncool = %d", param->d_coolrepeat, (int) param->d_coolsteps);
		j += sprintf_header_datafile_aux(header + j, smoothing_method, param);
		}

	#ifdef MULTICANONICAL_MODE
	j += sprintf(header + j, "mc_topcharge mc_weight");
	#endif
	j += sprintf(header + j, "\n");
	}


// open data file
FILE *open_file_with_header_replica(char const *const name, char const *const header,
                                    int const replica_index, GParam const *const param,
                                    int const binary_flag)
	{
	FILE *fp;
	char name_aux[STD_STRING_LENGTH];

	strcpy(name_aux, name);

	#ifdef REPLICA_MEAS_MODE
	if(param->d_N_replica_pt > 1)
		{
		char aux[STD_STRING_LENGTH];
		sprintf(aux, "_replica_%d", replica_index);
		REQUIRE(strlen(name) + strlen(aux) < STD_STRING_LENGTH, "filename too long");
		strcat(name_aux, aux);
		}
	#else
	(void) replica_index;
	#endif

	if(param->d_start == 2)
		{
		// open file in append mode
		fp = fopen(name_aux, "r");
		if(fp != NULL)
			{
			fclose(fp);
			fp = fopen(name_aux, "a");
			REQUIRE(fp != NULL, "failed to open %s for writing", name_aux);
			}
		else
			{
			fp = fopen(name_aux, "w");
			REQUIRE(fp != NULL, "failed to open %s for writing", name_aux);
			fputs(header, fp);
			}
		}
	else
		{
		// open file in write mode
		fp = fopen(name_aux, "w");
		REQUIRE(fp != NULL, "failed to open %s for writing", name_aux);
		fputs(header, fp);
		}

	fflush(fp);
	if(binary_flag == 1)
		{
		fclose(fp);
		fp = fopen(name_aux, "ab");
		REQUIRE(fp != NULL, "failed to open %s for writing in binary mode", name_aux);
		}
	return fp;
	}

// open data files
void open_data_files(Meas_Utils *meas_aux, int const replica_index, GParam const *const param)
	{
	char header[10 * STD_STRING_LENGTH];

	// data file
	header_datafile(header, param);
	meas_aux->datafilep = open_file_with_header_replica(param->d_data_file, header, replica_index, param, 0);

	// header for other files
	int j = sprintf(header, "# %d ", STDIM);
	for(int i = 0; i < STDIM; i++) j += sprintf(header + j, "%d ", param->d_size[i]);
	sprintf(header + j, "\n");

	// chiprime file
	if(param->d_chi_prime_meas == 1)
		meas_aux->chiprimefilep = open_file_with_header_replica(param->d_chiprime_file, header, replica_index, param, 0);

	// energy_slices file
	if(param->d_energy_slices_meas == 1)
		meas_aux->e_slices_filep = open_file_with_header_replica(param->d_energy_slices_file, header, replica_index, param, 0);

	// topcharge_tcorr file
	if(param->d_charge_slices_meas == 1 || param->d_charge_p_slices_meas == 1)
		meas_aux->q_slices_filep = open_file_with_header_replica(param->d_charge_slices_file, header, replica_index, param, 0);

	// energy_density file
	if(param->d_energy_density_meas == 1)
		meas_aux->energydensityfilep = open_file_with_header_replica(param->d_energydensity_file, header, replica_index, param, 1);

	// charge_density file
	if(param->d_charge_density_meas == 1)
		meas_aux->chargedensityfilep = open_file_with_header_replica(param->d_chargedensity_file, header, replica_index, param, 1);

	// polyakov_density files
	if(param->d_polyakov_density_meas == 1)
		{
		char filename[2 * STD_STRING_LENGTH];
		for(int i = 0; i < STDIM; i++)
			{
			sprintf(filename, "%s_dir%d", param->d_polyakovdensity_file, i);
			sprintf(header + j, "%d \n", i);
			meas_aux->polyakovdensityfilep[i] = open_file_with_header_replica(filename, header, replica_index, param, 1);
			}
		}
	}

// close data files
void close_data_files(Meas_Utils meas_aux, int const replica_index, GParam const *const param)
	{
	#ifndef REPLICA_MEAS_MODE
	if(replica_index != 0) return;
	#else
	(void) replica_index;
	#endif

	// data file
	fclose(meas_aux.datafilep);

	// chiprime file
	if(param->d_chi_prime_meas == 1)
		fclose(meas_aux.chiprimefilep);

	// energy_slices file
	if(param->d_energy_slices_meas == 1)
		fclose(meas_aux.e_slices_filep);

	// topcharge_tcorr file
	if(param->d_charge_slices_meas == 1)
		fclose(meas_aux.q_slices_filep);

	// clover_energy_density file
	if(param->d_energy_density_meas == 1)
		fclose(meas_aux.energydensityfilep);

	// charge_density file
	if(param->d_charge_density_meas == 1)
		fclose(meas_aux.chargedensityfilep);

	// polyakov_density files
	if(param->d_polyakov_density_meas == 1)
		for(int i = 0; i < STDIM; i++)
			fclose(meas_aux.polyakovdensityfilep[i]);
	}


void init_meas_utils(Meas_Utils *meas_aux, GParam const *const param, int const replica_index)
	{
	// max number of measures needed using any smoothing method
	int num_meas = param->d_agf_num_meas;
	if(num_meas < param->d_gf_num_meas)
		num_meas = param->d_gf_num_meas;
	if(num_meas < param->d_coolrepeat)
		num_meas = param->d_coolrepeat;

	if(num_meas > 0)
		{
		// allocate meas arrays
		if(param->d_plaquette_meas == 1)
			allocate_array_double(&(meas_aux->meanplaq), num_meas, __FILE__, __LINE__);

		if(param->d_clover_energy_meas == 1)
			allocate_array_double(&(meas_aux->clover_energy), num_meas, __FILE__, __LINE__);

		if(param->d_charge_meas == 1)
			allocate_array_double(&(meas_aux->charge), num_meas, __FILE__, __LINE__);

		if(param->d_polyakov_meas == 1)
			{
			allocate_array_double_pointer(&(meas_aux->polyre), num_meas, __FILE__, __LINE__);
			allocate_array_double_pointer(&(meas_aux->polyim), num_meas, __FILE__, __LINE__);
			for(int i = 0; i < num_meas; i++)
				{
				allocate_array_double(&(meas_aux->polyre[i]), STDIM, __FILE__, __LINE__);
				allocate_array_double(&(meas_aux->polyim[i]), STDIM, __FILE__, __LINE__);
				}
			}

		if(param->d_multipolyakov_order >= 1)
			{
			allocate_array_double(&(meas_aux->multipolyre), num_meas, __FILE__, __LINE__);
			allocate_array_double(&(meas_aux->multipolyim), num_meas, __FILE__, __LINE__);
			}

		if(param->d_chi_prime_meas == 1)
			allocate_array_double(&(meas_aux->chi_prime), num_meas, __FILE__, __LINE__);

		if(param->d_charge_prime_meas == 1)
			{
			allocate_array_double_pointer(&(meas_aux->charge_prime), num_meas, __FILE__, __LINE__);
			for(int i = 0; i < num_meas; i++)
				allocate_array_double(&(meas_aux->charge_prime[i]), STDIM, __FILE__, __LINE__);
			}

		if(param->d_action_meas == 1)
			{
			allocate_array_double(&(meas_aux->action1), num_meas, __FILE__, __LINE__);
			allocate_array_double(&(meas_aux->action2), num_meas, __FILE__, __LINE__);
			allocate_array_double(&(meas_aux->action3), num_meas, __FILE__, __LINE__);
			allocate_array_double(&(meas_aux->potential), num_meas, __FILE__, __LINE__);
			}

		// allocate auxiliary lattices
		for(int i = 0; i < 4; i++)
			{
			allocate_array_GAUGE_GROUP_pointer(&(meas_aux->lattice_aux[i]), param->d_volume, __FILE__, __LINE__);
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < (param->d_volume); r++)
				{
				allocate_array_GAUGE_GROUP(&(meas_aux->lattice_aux[i][r]), STDIM, __FILE__, __LINE__);
				}
			}
		}

	// allocate arrays for density profiles
	if(param->d_energy_slices_meas == 1 || param->d_charge_slices_meas == 1 || param->d_charge_p_slices_meas == 1)
		allocate_array_double(&(meas_aux->real_slices), param->d_max_size, __FILE__, __LINE__);

	if(param->d_charge_p_slices_meas == 1)
		allocate_array_double(&(meas_aux->imag_slices), param->d_max_size, __FILE__, __LINE__);

	if(param->d_energy_density_meas == 1 || param->d_charge_density_meas == 1 || param->d_charge_p_slices_meas == 1)
		allocate_array_double(&(meas_aux->scalar_density), param->d_volume, __FILE__, __LINE__);

	if(param->d_polyakov_density_meas == 1)
		{
		long max_space_vol = 0;
		for(int i = 0; i < STDIM; i++) if(param->d_space_vol[i] > max_space_vol) max_space_vol = param->d_space_vol[i];
		allocate_array_double(&(meas_aux->polyre_density), max_space_vol, __FILE__, __LINE__);
		allocate_array_double(&(meas_aux->polyim_density), max_space_vol, __FILE__, __LINE__);
		}

	// open data files
	open_data_files(meas_aux, replica_index, param);
	}

void init_meas_utils_replica(Meas_Utils **meas_aux, GParam const *const param)
	{
	allocate_array_Meas_Utils(meas_aux, param->d_N_replica_pt, __FILE__, __LINE__);

	// init meas utils and data files for physical replica
	init_meas_utils(&((*meas_aux)[0]), param, 0);

	// init meas utils for other replicas if using replica meas mode
	#ifdef REPLICA_MEAS_MODE
	for(int i = 1; i < param->d_N_replica_pt; i++)
		init_meas_utils(&((*meas_aux)[i]), param, i);
	#endif
	}

void free_meas_utils(Meas_Utils meas_aux, GParam const *const param, int const replica_index)
	{
	int num_meas = param->d_agf_num_meas;
	if(num_meas < param->d_gf_num_meas)
		num_meas = param->d_gf_num_meas;
	if(num_meas < param->d_coolrepeat)
		num_meas = param->d_coolrepeat;

	if(num_meas > 0)
		{
		// free meas arrays
		if(param->d_plaquette_meas == 1)
			free(meas_aux.meanplaq);

		if(param->d_clover_energy_meas == 1)
			free(meas_aux.clover_energy);

		if(param->d_charge_meas == 1)
			free(meas_aux.charge);

		if(param->d_polyakov_meas == 1)
			{
			for(int i = 0; i < num_meas; i++)
				{
				free(meas_aux.polyre[i]);
				free(meas_aux.polyim[i]);
				}
			free(meas_aux.polyre);
			free(meas_aux.polyim);
			}

		if(param->d_multipolyakov_order >= 1)
			{
			free(meas_aux.multipolyre);
			free(meas_aux.multipolyim);
			}

		if(param->d_chi_prime_meas == 1)
			free(meas_aux.chi_prime);

		if(param->d_charge_prime_meas == 1)
			{
			for(int i = 0; i < num_meas; i++)
				free(meas_aux.charge_prime[i]);
			free(meas_aux.charge_prime);
			}

		if(param->d_action_meas == 1)
			{
			free(meas_aux.action1);
			free(meas_aux.action2);
			free(meas_aux.action3);
			free(meas_aux.potential);
			}

		// free auxiliary lattices
		for(int i = 0; i < 4; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r = 0; r < (param->d_volume); r++)
				{
				free(meas_aux.lattice_aux[i][r]);
				}
			free(meas_aux.lattice_aux[i]);
			}
		}

	// free arrays for density profiles
	if(param->d_energy_slices_meas == 1 || param->d_charge_slices_meas == 1 || param->d_charge_p_slices_meas == 1)
		free(meas_aux.real_slices);

	if(param->d_charge_p_slices_meas == 1)
		free(meas_aux.imag_slices);

	if(param->d_energy_density_meas == 1 || param->d_charge_density_meas == 1 || param->d_charge_p_slices_meas == 1)
		free(meas_aux.scalar_density);

	if(param->d_polyakov_density_meas == 1)
		{
		free(meas_aux.polyre_density);
		free(meas_aux.polyim_density);
		}

	// close data files
	close_data_files(meas_aux, replica_index, param);
	}

void free_meas_utils_replica(Meas_Utils *meas_aux, GParam const *const param)
	{
	// free meas utils for physical replica
	free_meas_utils(meas_aux[0], param, 0);

	// free meas utils for other replicas if using replica meas mode
	#ifdef REPLICA_MEAS_MODE
	for(int i = 1; i < param->d_N_replica_pt; i++)
		free_meas_utils(meas_aux[i], param, i);
	#endif

	free(meas_aux);
	}

void print_measures_aux(int const num_meas, long const update_index, GParam const *const param, Meas_Utils const *const meas_aux)
	{
	double time_step;
	if(param->d_agf_meas_each > 0.0) time_step = (param->d_agf_meas_each);
	else time_step = param->d_ngfsteps;

	for(int i = 0; i < num_meas; i++)
		{
		if(param->d_plaquette_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->meanplaq[i]);
		if(param->d_clover_energy_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->clover_energy[i]);
		if(param->d_charge_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->charge[i]);
		if(param->d_polyakov_meas == 1) for(int j = 0; j < STDIM; j++) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", meas_aux->polyre[i][j], meas_aux->polyim[i][j]);
		if(param->d_multipolyakov_order >= 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", meas_aux->multipolyre[i], meas_aux->multipolyim[i]);
		if(param->d_chi_prime_meas == 1) fprintf(meas_aux->chiprimefilep, "%ld % 18.12e % 18.12e\n", update_index, (i + 1) * time_step, meas_aux->chi_prime[i]);
		if(param->d_charge_prime_meas == 1) for(int j = 0; j < STDIM; j++) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->charge_prime[i][j]);
		if(param->d_action_meas == 1)
			{
			fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->action1[i]);
			fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->action2[i]);
			fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->action3[i]);
			fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->potential[i]);
			}
		}
	}

#endif
