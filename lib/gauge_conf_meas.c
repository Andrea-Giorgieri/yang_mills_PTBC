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

// computation of the plaquette (1/NCOLOR the trace of, with twist factor) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param,
					long r, int i, int j)
	{
	GAUGE_GROUP matrix;
	
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(j >= STDIM || i >= STDIM)
		{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//		 ^ i
//		 |	(2)
//		 +---<---+
//		 |		 |
//	(3) V		 ^ (1)
//		 |		 |
//		 +--->---+---> j
//		 r	(4)
//

	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\nu\mu(x) = conj(Z_\mu\nu(x))
	
	return retr(&matrix);
	}


// computation of the plaquette (1/NCOLOR the trace of, with twist factor) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
											 Geometry const * const geo,
											 GParam const * const param,
											 long r,
											 int i,
											 int j)
	{
	GAUGE_GROUP matrix;
	
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(j >= STDIM || i >= STDIM)
		{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//		 ^ i
//		 |	(2)
//		 +---<---+
//		 |		 |
//	(3) V		 ^ (1)
//		 |		 |
//		 +--->---+---> j
//		 r	(4)
//

	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\mu\nu(x) = conj(Z_\nu\mu(x))
	
	return retr(&matrix)+I*imtr(&matrix);
	}


// computation of the plaquette with twist factor (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						long r,
						int i,
						int j,
						GAUGE_GROUP *matrix)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(j >= STDIM || i >= STDIM)
		{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//		 ^ j
//		 |	(3)
//		 +---<---+
//		 |		 |
//	(4) V		 ^ (2)
//		 |		 |
//		 +--->---+---> i
//		 r	(1)
//

	equal(matrix, &(GC->lattice[r][i]));
	times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\mu\nu(x) = conj(Z_\nu\mu(x))
	}


// compute the four-leaf clover in position r, in the plane i,j and save it in M 
void clover(Gauge_Conf const * const GC,
				Geometry const * const geo,
				GParam const * const param,
				long r,
				int i,
				int j,
				GAUGE_GROUP *M)
	{
	GAUGE_GROUP aux;
	long k, p;
	
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM || j >= STDIM)
		{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#else
	(void) param;
	#endif
	
	zero(M);

//
//					  j ^
//						|
//				 (14)	|		(3)
//			+-----<-----++-----<-----+
//			|			||			 |
//			|			||			 |
//	   (15) V (13)      ^V (4)		 ^ (2)
//			|			||			 |
//			|	 (16)	|| r  (1)	 |
//	     p  +----->-----++----->-----+------> i
//			+-----<-----++-----<-----+
//			|	 (9)	||	  (8)	 |
//			|			||			 |
//	   (10) V (12)	    ^V (5)		 ^ (7)
//			|			||			 |
//			|			||			 |
//			+------>----++----->-----+
//				 (11)	 k	  (6)
//
	// avanti-avanti
	equal(&aux, &(GC->lattice[r][i]) );							// 1
	times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );		// 2
	times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );	// 3
	times_equal_dag(&aux, &(GC->lattice[r][j]) );				// 4
	times_equal_complex(&aux, GC->Z[r][dirs_to_si(i,j)]);		// twist anticlockwise
	plus_equal(M, &aux);

	k=nnm(geo, r, j);

	// avanti-indietro
	equal_dag(&aux, &(GC->lattice[k][j]) );						// 5
	times_equal(&aux, &(GC->lattice[k][i]) );					// 6
	times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );		// 7
	times_equal_dag(&aux, &(GC->lattice[r][i]) );				// 8
	times_equal_complex(&aux, GC->Z[k][dirs_to_si(i,j)]);		// twist anticlockwise
	plus_equal(M, &aux);

	p=nnm(geo, r, i);

	// indietro-indietro
	equal_dag(&aux, &(GC->lattice[p][i]) );								// 9
	times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );			// 10
	times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );				// 11
	times_equal(&aux, &(GC->lattice[k][j]) );							// 12
	times_equal_complex(&aux, GC->Z[nnm(geo, k, i)][dirs_to_si(i,j)]);	// twist anticlockwise
	plus_equal(M, &aux);

	// indietro-avanti
	equal(&aux, &(GC->lattice[r][j]) );							// 13
	times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );	// 14
	times_equal_dag(&aux, &(GC->lattice[p][j]) );				// 15
	times_equal(&aux, &(GC->lattice[p][i]) );					// 16
	times_equal_complex(&aux, GC->Z[p][dirs_to_si(i,j)]);		// twist anticlockwise
	plus_equal(M, &aux);
	}


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					double *plaqs,
					double *plaqt)
	{
	long r;
	double ps, pt;
	
	ps=0.0;
	pt=0.0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		int i, j;
		i=0;
		for(j=1; j<STDIM; j++)
			{
			pt+=plaquettep(GC, geo, param, r, i, j);
			}
		
		for(i=1; i<STDIM; i++)
			{
			for(j=i+1; j<STDIM; j++)
				{
				ps+=plaquettep(GC, geo, param, r, i, j);
				}
			}
		}
	
	if(STDIM>2)
		{
		ps*=param->d_inv_vol;
		ps/=((double) (STDIM-1)*(STDIM-2)/2);
		}
	else
		{
		ps=0.0;
		}
	
	pt*=param->d_inv_vol;
	pt/=((double) STDIM-1);
	
	*plaqs=ps;
	*plaqt=pt;
	}


// compute the clover discretization of
// sum_{\mu\nu}	Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
								Geometry const * const geo,
								GParam const * const param,
								double *energy)
	{
	long r;
	double ris;
	
	ris=0.0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		int i, j;
		GAUGE_GROUP aux1, aux2;
		
		for(i=0; i<STDIM; i++)
			{
			for(j=i+1; j<STDIM; j++)
				{
				clover(GC, geo, param, r, i, j, &aux1);
				ta(&aux1);
				equal(&aux2, &aux1);
				times_equal(&aux1, &aux2);
				ris+=-NCOLOR*retr(&aux1)/16.0;
				}
			}
		}
	*energy=ris*param->d_inv_vol;
	}


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					double *repoly,
					double *impoly)
	{
	long rsp;
	double rep, imp;

	rep=0.0;
	imp=0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
	#endif
	for(rsp=0; rsp<param->d_space_vol; rsp++)
		{
		long r;
		int i;
		GAUGE_GROUP matrix;

		r=sisp_and_t_to_si(geo, rsp, 0);

		one(&matrix);
		for(i=0; i<param->d_size[0]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][0]));
			r=nnp(geo, r, 0);
			}

		rep+=retr(&matrix);
		imp+=imtr(&matrix);
		}

	*repoly=rep*param->d_inv_space_vol;
	*impoly=imp*param->d_inv_space_vol;
	}


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						double *repoly,
						double *impoly)
	{
	long rsp;
	double rep, imp;
	double complex tr;

	rep=0.0;
	imp=0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
	#endif
	for(rsp=0; rsp<param->d_space_vol; rsp++)
		{
		long r;
		int i;
		GAUGE_GROUP matrix;

		r=sisp_and_t_to_si(geo, rsp, 0);

		one(&matrix);
		for(i=0; i<param->d_size[0]; i++)
			{
			times_equal(&matrix, &(GC->lattice[r][0]));
			r=nnp(geo, r, 0);
			}
		tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

		#if NCOLOR==1
			rep+=0.0;
		#else
			rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
		#endif

		imp+=0.0;
		}

	*repoly=rep*param->d_inv_space_vol;
	*impoly=imp*param->d_inv_space_vol;
	}



// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_with_tracedef(Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							double *repoly,
							double *impoly)
	{
	long rsp;
	double **rep, **imp;
	int j;
	long i;

	for(j=0;j<(int)floor(NCOLOR/2);j++)
		{
		repoly[j]=0.0;
		impoly[j]=0.0;
		}
	
	allocate_array_double_pointer(&rep, param->d_space_vol, __FILE__, __LINE__);
	allocate_array_double_pointer(&imp, param->d_space_vol, __FILE__, __LINE__);
	for(i=0; i<param->d_space_vol; i++)
		{
		allocate_array_double(&(rep[i]), (int)floor(NCOLOR/2), __FILE__, __LINE__);
		allocate_array_double(&(imp[i]), (int)floor(NCOLOR/2), __FILE__, __LINE__);
		}

	for(i=0; i<param->d_space_vol; i++)
		{
		for(j=0; j<(int)floor(NCOLOR/2); j++)
			{
			rep[i][j] = 0.0;
			imp[i][j] = 0.0;
			}
		}

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(rsp)
	#endif
	for(rsp=0; rsp<param->d_space_vol; rsp++)
		{
		long r;
		int k;
		GAUGE_GROUP matrix, matrix2;

		r=sisp_and_t_to_si(geo, rsp, 0);

		one(&matrix);
		for(k=0; k<param->d_size[0]; k++)
			{
			times_equal(&matrix, &(GC->lattice[r][0]));
			r=nnp(geo, r, 0);
			}

		rep[rsp][0] = retr(&matrix);
		imp[rsp][0] = imtr(&matrix);

		equal(&matrix2, &matrix);

		for(k=1; k<(int)floor(NCOLOR/2.0); k++)
			{
			times_equal(&matrix2, &matrix);
			rep[rsp][k] = retr(&matrix2);
			imp[rsp][k] = imtr(&matrix2);
			}
		}

	for(j=0; j<(int)floor(NCOLOR/2); j++)
		{
		for(i=0; i<param->d_space_vol; i++)
			{
			repoly[j] += rep[i][j];
			impoly[j] += imp[i][j];
			}
		}

	for(j=0; j<(int)floor(NCOLOR/2.0); j++)
		{
		repoly[j] *= param->d_inv_space_vol;
		impoly[j] *= param->d_inv_space_vol;
		}

	for(i=0; i<param->d_space_vol; i++)
		{
		free(rep[i]);
		free(imp[i]);
		}
	free(rep);
	free(imp);
	}


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					long r)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		(void) GC;
		(void) geo;
		(void) param;
		(void) r;
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	double ris;
	
	#if (STDIM==4 && NCOLOR>1)
	GAUGE_GROUP aux1, aux2, aux3;
	double real1, real2, loc_charge;
	const double chnorm=1.0/(128.0*PI*PI);
	int i, dir[4][3], sign;
	
	dir[0][0] = 0;
	dir[0][1] = 0;
	dir[0][2] = 0;
	
	dir[1][0] = 1;
	dir[1][1] = 2;
	dir[1][2] = 3;
	
	dir[2][0] = 2;
	dir[2][1] = 1;
	dir[2][2] = 1;
	
	dir[3][0] = 3;
	dir[3][1] = 3;
	dir[3][2] = 2;
	
	sign=-1;
	loc_charge=0.0;
	
	for(i=0; i<3; i++)
		{
		clover(GC, geo, param, r, dir[0][i], dir[1][i], &aux1);
		clover(GC, geo, param, r, dir[2][i], dir[3][i], &aux2);
		
		times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
		real1=retr(&aux3)*NCOLOR;
		
		times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
		real2=retr(&aux3)*NCOLOR;
		
		loc_charge+=((double) sign*(real1-real2));
		sign=-sign;
		}
	ris=(loc_charge*chnorm);
	#endif
	
	#if (STDIM==2 && NCOLOR==1)
	GAUGE_GROUP u1matrix;
	double angle;
	
	plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
	angle=atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;
	
	ris=angle;
	#endif
	
	return ris;
	}

// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	double ris;
	long r;

	ris=0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		ris+=loc_topcharge(GC, geo, param, r);
		}

	return ris;
	}

// sum loc_topcharge over STDIM-1 dirs and then the abs value of the result over the remaining dir
double topcharge_prime(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					int const dir)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	int i, cartcoord[STDIM];
	double ris, *tmp;
	long r;
	
	ris=0.0;
	allocate_array_double(&tmp, param->d_size[dir], __FILE__, __LINE__);
	for(i=0; i<param->d_size[dir]; i++) tmp[i] = 0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel num_threads(NTHREADS)
	#endif
		{
		int j;
		double *tmp_private;
		
		allocate_array_double(&tmp_private, param->d_size[dir], __FILE__, __LINE__);
		for(j=0; j<param->d_size[dir]; j++) tmp_private[j] = 0;
		
		#ifdef OPENMP_MODE
		#pragma omp for private(r, cartcoord) 
		#endif
		for(r=0; r<(param->d_volume); r++)
			{
			si_to_cart(cartcoord,r,param);
			tmp_private[cartcoord[dir]]+=loc_topcharge(GC, geo, param, r);
			}
		#ifdef OPENMP_MODE
		#pragma omp critical 
		#endif
			{
			for(j=0; j<param->d_size[dir]; j++) tmp[j] += tmp_private[j];
			free(tmp_private);
			}
		}

	for(i=0; i<param->d_size[dir]; i++) ris += fabs(tmp[i]);
	free(tmp);
	return ris;
	}

// chi^\prime = (1/8) int d^4x |x|^2 <q(x)q(0)> = < (1/8) int d^4x |x|^2 q(x) q(0) > = < G2 >
// This function computes the quantity (q(0)/8) sum_{x} d(x,0)^2 q(x) = a^2 G2, whose mean over the ensamble is a^2 chi^\prime
// d(x,y) = lattice distance between sites x and y keeping periodic boundary conditions into account (i.e., the shortest distance between x and y)
double topo_chi_prime(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	double ris=0.0, factor=0.125; // factor = 1/(2D) = 1/8
	long r;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		double d2 = square_distance(r, 0, param); // d(r,0)^2
		ris += d2 * loc_topcharge(GC, geo, param, r);
		}
	ris *=	loc_topcharge(GC, geo, param, 0) * factor; // ris *= q(0) / 8
	
	return ris;
	}


void topcharge_timeslices(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param, double *ris, int ncool, FILE *topchar_tcorr_filep)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	long r;
	int N_t = param->d_size[0];

	for (int i=0; i<N_t; i++) ris[i]=0.0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		int t = geo->d_timeslice[r];
		ris[t] += loc_topcharge(GC, geo, param, r);
		}

	fprintf(topchar_tcorr_filep, "%ld %d ", GC->update_index, ncool);
	for (int i=0; i<param->d_size[0]; i++) fprintf(topchar_tcorr_filep, " % 18.12e", ris[i]);
	fprintf(topchar_tcorr_filep, "\n");
	}


void topcharge_timeslices_cooling(Gauge_Conf *const GC,
						Geometry const * const geo,
						GParam const * const param,
						Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	// measure no cooling
	topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, 0, meas_aux->topchar_tcorr_filep); 
	
	if(param->d_coolrepeat>0)	// if using cooling
		{	
		int iter;

		// measure with cooling
		for(iter=0; iter<(param->d_coolrepeat); iter++)
			{
			cooling(GC, geo, param, param->d_coolsteps);
			topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, (iter+1)*param->d_coolsteps, meas_aux->topchar_tcorr_filep);
			}
		restore_gauge_conf(GC, param);
		}

	fflush(meas_aux->topchar_tcorr_filep);
	}


void topcharge_timeslices_gradflow(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	// measure no gradient flow
	topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, 0, meas_aux->topchar_tcorr_filep);
	
	if(param->d_gf_num_meas > 0)	// if using gradient flow
		{
		int count;

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			
			if ( (count % param->d_gf_meas_each) == 0)
				{
				topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, count, meas_aux->topchar_tcorr_filep);
				}
			}
		restore_gauge_conf(GC, param);		
		}

	fflush(meas_aux->topchar_tcorr_filep);
	}


// compute topological observables (Q, chi_prime) after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topo_obs_cooling(Gauge_Conf * const GC,
					Geometry const * const geo,
					GParam const * const param,
					Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_coolsteps>0)	// if using cooling
		{	
		double plaqs, plaqt;
		int iter;
	
		for(iter=0; iter<(param->d_coolrepeat); iter++)
			{
			cooling(GC, geo, param, param->d_coolsteps);
			
			meas_aux->charge[iter] = topcharge(GC, geo, param);
			meas_aux->chi_prime[iter] = topo_chi_prime(GC, geo, param);
			
			plaquette(GC, geo, param, &plaqs, &plaqt);
			#if(STDIM==4)
			meas_aux->meanplaq[iter]=0.5*(plaqs+plaqt);
			#else
			meas_aux->meanplaq[iter]=plaqt;
			#endif
			}
		restore_gauge_conf(GC, param);
		}
	else	// no cooling
		{
		double ris, ris2, plaqs, plaqt; 
		int iter;

		ris=topcharge(GC, geo, param);
		ris2=topo_chi_prime(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
	
		for(iter=0; iter<(param->d_coolrepeat); iter++)
			{
			meas_aux->charge[iter]=ris;
			meas_aux->chi_prime[iter]=ris2;
			#if(STDIM==4)
			meas_aux->meanplaq[iter]=0.5*(plaqs+plaqt);
			#else
			meas_aux->meanplaq[iter]=plaqt;
			#endif
			}
		}
	}
	

void topo_obs_gradflow(Gauge_Conf * const GC,
						Geometry const * const geo,
						GParam const * const param,
						Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_ngfsteps>0)	// if using gradient flow
		{	
		double plaqs, plaqt;
		int count, meas_count;
		
		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			
			if ( (count % param->d_gf_meas_each) == 0)
				{
				meas_count = count/param->d_gf_meas_each-1;
				meas_aux->charge[meas_count]=topcharge(GC, geo, param);
				meas_aux->chi_prime[meas_count]=topo_chi_prime(GC, geo, param);
				plaquette(GC, geo, param, &plaqs, &plaqt);
				#if(STDIM==4)
					meas_aux->meanplaq[meas_count]=0.5*(plaqs+plaqt);
				#else
					meas_aux->meanplaq[meas_count]=plaqt;
				#endif
				}
			}
		restore_gauge_conf(GC, param);
		}
	else	// no gradient flow
		{
		double plaqs, plaqt; 
		
		meas_aux->charge[0]=topcharge(GC, geo, param);
		meas_aux->chi_prime[0]=topo_chi_prime(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
		#if(STDIM==4)
			meas_aux->meanplaq[0]=0.5*(plaqs+plaqt);
		#else
			meas_aux->meanplaq[0]=plaqt;
		#endif
		}
	}

	
void topo_obs_clover_energy_gradflow(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_ngfsteps>0)	// if using gradient flow
		{
		double tmp_energy;
		int count, meas_count;
		
		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			
			if ( (count % param->d_gf_meas_each) == 0)
				{
				meas_count = count/param->d_gf_meas_each-1;
				meas_aux->charge[meas_count]=topcharge(GC, geo, param);
				meas_aux->chi_prime[meas_count]=topo_chi_prime(GC, geo, param);
				clover_disc_energy(GC, geo, param, &tmp_energy);
				meas_aux->clover_energy[meas_count] = tmp_energy;
				}
			}
		restore_gauge_conf(GC, param);
		}
	else	// no gradient flow
		{
		meas_aux->charge[0]=topcharge(GC, geo, param);
		meas_aux->chi_prime[0]=topo_chi_prime(GC, geo, param);
		clover_disc_energy(GC, geo, param, &(meas_aux->clover_energy[0]));
		}
	}


/*---------------------------------------------*/
// OBSERVABLE NEEDED JUST TO CHECK HOW COOLING DESTROYS TOPOLOGICAL CORRELATIONS
void check_correlation_decay_cooling(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, double *ratio)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_coolsteps>0)	// if using cooling
		{	
		Gauge_Conf helperconf; 
		double Q, satd;
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		for(int iter=0; iter<(param->d_coolrepeat); iter++)
			{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			Q = fabs(topcharge(&helperconf, geo, param));
			satd = sum_abs_topcharge_dens(&helperconf, geo, param);
			ratio[iter] = (satd-Q)/satd;
			}
		free_gauge_conf(&helperconf, param);
		}
	}


double sum_abs_topcharge_dens(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param)
	{
	double sum=0.0;
	for (long r=0; r<(param->d_volume); r++)
		{
		sum += fabs(loc_topcharge(GC, geo, param, r));
		}
	return sum;
	}

/*---------------------------------------------*/

// compute the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf * const GC,
						Geometry const * const geo,
						GParam const * const param,
						Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	if(param->d_coolsteps>0)	// if using cooling
		{ 
		double plaqs, plaqt;
		int iter;
	
		for(iter=0; iter<(param->d_coolrepeat); iter++)
			{
			cooling(GC, geo, param, param->d_coolsteps);
			
			meas_aux->charge[iter] = topcharge(GC, geo, param);
			
			plaquette(GC, geo, param, &plaqs, &plaqt);
			#if(STDIM==4)
				meas_aux->meanplaq[iter]=0.5*(plaqs+plaqt);
			#else
				meas_aux->meanplaq[iter]=plaqt;
			#endif
			}

		restore_gauge_conf(GC, param);
		}
	else	// no cooling
		{
		double ris, plaqs, plaqt; 
		int iter;

		ris=topcharge(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
	
		for(iter=0; iter<(param->d_coolrepeat); iter++)
			{
			meas_aux->charge[iter]=ris;
			#if(STDIM==4)
				meas_aux->meanplaq[iter]=0.5*(plaqs+plaqt);
			#else
				meas_aux->meanplaq[iter]=plaqt;
			#endif
			}
		} 
	}


void topcharge_gradflow(Gauge_Conf * const GC,
						Geometry const * const geo,
						GParam const * const param,
						Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_ngfsteps>0)	// if using gradient flow
		{ 
		double plaqs, plaqt;
		int count, meas_count;

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			
			if ( (count % param->d_gf_meas_each) == 0)
				{
				meas_count = count/param->d_gf_meas_each-1;
				meas_aux->charge[meas_count]=topcharge(GC, geo, param);
				plaquette(GC, geo, param, &plaqs, &plaqt);
				#if(STDIM==4)
					meas_aux->meanplaq[meas_count]=0.5*(plaqs+plaqt);
				#else
					meas_aux->meanplaq[meas_count]=plaqt;
				#endif
				}
			}
		restore_gauge_conf(GC, param);
		}
	else	// no gradient flow
		{
		double plaqs, plaqt; 
		
		meas_aux->charge[0]=topcharge(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
		#if(STDIM==4)
			meas_aux->meanplaq[0]=0.5*(plaqs+plaqt);
		#else
			meas_aux->meanplaq[0]=plaqt;
		#endif
		}
	}

	
void topcharge_clover_energy_gradflow(Gauge_Conf * const GC,
								Geometry const * const geo,
								GParam const * const param,
								Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	if(param->d_ngfsteps>0)	// if using gradient flow
		{
		double tmp_energy;
		int count, meas_count;

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			
			if ( (count % param->d_gf_meas_each) == 0)
				{
				meas_count = count/param->d_gf_meas_each-1;
				meas_aux->charge[meas_count]=topcharge(GC, geo, param);
				clover_disc_energy(GC, geo, param, &tmp_energy);
				meas_aux->clover_energy[meas_count] = tmp_energy;
				}
			}
		restore_gauge_conf(GC, param);
		}
	else	// no gradient flow
		{
		double tmp_energy;
		
		meas_aux->charge[0]=topcharge(GC, geo, param);
		clover_disc_energy(GC, geo, param, &tmp_energy);
		meas_aux->clover_energy[0] = tmp_energy;
		}
	}


// compute the correlator of the local topological charge
// after "ncool" cooling steps up to spatial distance "dist"
void loc_topcharge_corr(Gauge_Conf * const GC,
						Geometry const * const geo,
						GParam const * const param,
						int ncool,
						int dist,
						double *ris)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	double *topch;
	long r;
	int i;
	
	// TO DO: refactor with meas_aux
	allocate_array_double(&topch, param->d_volume, __FILE__, __LINE__);
	
	// compute the local topological charge
	if(ncool>0)
		{		
		// cooling
		cooling(GC, geo, param, ncool);
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			topch[r]=loc_topcharge(GC, geo, param, r);
			}
		
		restore_gauge_conf(GC, param);
		}
	else
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			topch[r]=loc_topcharge(GC, geo, param, r);
			}
		}
	
	// compute correlators
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif
	for(i=0; i<dist; i++)
		{
		double aux;
		long r1, r2;
		int j, dir;
		
		ris[i]=0.0;
		
		for(r1=0; r1<param->d_volume; r1++)
			{
			aux=0.0;
			
			for(dir=1; dir<STDIM; dir++)
				{
				r2=r1;
				for(j=0; j<i; j++)
					{
					r2=nnp(geo, r2, dir);
					}
				aux+=topch[r2];
				}
			aux/=(double)(STDIM-1);
			
			ris[i]+=aux*topch[r1];
			}
		ris[i]*=param->d_inv_vol;
		}
	
	// free memory
	free(topch);
	}


void perform_measures_aux(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
							int const meas_count, Meas_Utils *meas_aux)
	{
	if (param->d_plaquette_meas == 1 ) 
		{
		double plaqs, plaqt;
		plaquette(GC, geo, param, &plaqs, &plaqt);
		#if(STDIM==4)
		meas_aux->meanplaq[meas_count]=0.5*(plaqs+plaqt);
		#else
		meas_aux->meanplaq[meas_count]=plaqt;
		#endif
		}
	if (param->d_clover_energy_meas   == 1) clover_disc_energy(GC, geo, param, &(meas_aux->clover_energy[meas_count]));
	if (param->d_charge_meas          == 1) meas_aux->charge[meas_count] = topcharge(GC, geo, param);
	if (param->d_topcharge_tcorr_meas == 1) topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, meas_count+1, meas_aux->topchar_tcorr_filep);				
	if (param->d_chi_prime_meas       == 1) meas_aux->chi_prime[meas_count] = topo_chi_prime(GC, geo, param);
	if (param->d_charge_prime_meas    == 1) for (int i=0; i<STDIM; i++) meas_aux->charge_prime[meas_count][i] = topcharge_prime(GC, geo, param, i);
	}


void perform_measures_localobs(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
								Meas_Utils *meas_aux)
	{
	int i;
	double plaqs=0.0, plaqt=0.0, polyre=0.0, polyim=0.0, clover_energy=0.0, charge=0.0, chi_prime=0.0, charge_prime[STDIM]; // =0.0 to suppress gcc warning
	
	// perform meas
	if (param->d_plaquette_meas       == 1) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas   == 1) clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_charge_meas          == 1) charge = topcharge(GC, geo, param);
	if (param->d_polyakov_meas        == 1) polyakov(GC, geo, param, &polyre, &polyim);
	if (param->d_chi_prime_meas       == 1) chi_prime = topo_chi_prime(GC, geo, param);
	if (param->d_charge_prime_meas    == 1) for (i=0; i<STDIM; i++) charge_prime[i] = topcharge_prime(GC, geo, param, i);
	if (param->d_topcharge_tcorr_meas == 1) topcharge_timeslices(GC, geo, param, meas_aux->sum_q_timeslices, 0, meas_aux->topchar_tcorr_filep);
	
	//print meas (topcharge_tcorr_timeslices already printed by topcharge_timeslices())
	fprintf(meas_aux->datafilep, "%ld ", GC->update_index);
	if (param->d_plaquette_meas     == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", clover_energy);
	if (param->d_charge_meas        == 1) fprintf(meas_aux->datafilep, "% 18.12e ", charge);
	if (param->d_polyakov_meas      == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", polyre, polyim);
	if (param->d_chi_prime_meas     == 1) fprintf(meas_aux->chiprimefilep, "%ld 0 % 18.12e\n", GC->update_index, chi_prime);
	if (param->d_charge_prime_meas  == 1) for (i=0; i<STDIM; i++) fprintf(meas_aux->datafilep, "% 18.12e ", charge_prime[i]);
	
	// flush data files
	fflush(meas_aux->datafilep);
	if (param->d_topcharge_tcorr_meas == 1) fflush(meas_aux->topchar_tcorr_filep);
	if (param->d_chi_prime_meas       == 1) fflush(meas_aux->chiprimefilep);
	}


void perform_measures_localobs_notopo(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
										Meas_Utils *meas_aux)
	{
	double plaqs=0.0, plaqt=0.0, clover_energy=0.0, polyre=0.0, polyim=0.0;
	
	if (param->d_plaquette_meas     == 1) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1) clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_polyakov_meas      == 1) polyakov(GC, geo, param, &polyre, &polyim);
	
	if (param->d_plaquette_meas     == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", clover_energy);
	if (param->d_polyakov_meas      == 1) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", polyre, polyim);
	
	fflush(meas_aux->datafilep);
	}


void perform_measures_localobs_cooling(Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )

	// meas no cooling
	perform_measures_localobs(GC, geo, param, meas_aux);
	
	// meas cooling
	if (param->d_topcharge_tcorr_meas == 1) topcharge_timeslices_cooling(GC, geo, param, meas_aux);
	if (param->d_chi_prime_meas       == 1) topo_obs_cooling(GC, geo, param, meas_aux);
	else                                    topcharge_cooling(GC, geo, param, meas_aux);
	
	// print meas cooling
	for(int i=0; i<param->d_coolrepeat; i++)
		{
		fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", meas_aux->charge[i], meas_aux->meanplaq[i]);
		if (param->d_chi_prime_meas == 1 )
			fprintf(meas_aux->chiprimefilep, "%ld %d % 18.12e\n", GC->update_index, (i+1)*param->d_coolsteps, meas_aux->chi_prime[i]);
		}
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if (param->d_chi_prime_meas == 1) fflush(meas_aux->chiprimefilep);
	
	#else
	perform_measures_localobs_notopo(GC, geo, param, meas_aux);
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	#endif
	}


void perform_measures_localobs_with_gradflow(Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											Meas_Utils *meas_aux)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int gradflowrepeat;
	
	// meas no gradflow
	perform_measures_localobs(GC, geo, param, meas_aux);
	
	// meas gradflow
	gradflowrepeat = (int)(param->d_ngfsteps/param->d_gf_meas_each);
	if (gradflowrepeat > 0)
		{
		int count, meas_count;

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
			{
			gradflow_RKstep(GC, geo, param, param->d_gfstep, meas_aux);
			if (count % param->d_gf_meas_each == 0)
				{
				meas_count = count/param->d_gf_meas_each-1;
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				}
			}
		
		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);
		
		// print meas gradflow
		print_measures_arrays(gradflowrepeat, GC->update_index, param, meas_aux);
		}
	
	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	fprintf(meas_aux->datafilep, "% 18.12e ", GC->stored_topcharge);
	#endif
	
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if (param->d_topcharge_tcorr_meas == 1) fflush(meas_aux->topchar_tcorr_filep);
	if (param->d_chi_prime_meas       == 1) fflush(meas_aux->chiprimefilep);
	
	#else
	perform_measures_localobs_notopo(GC, geo, param, meas_aux);
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	#endif
	}


void perform_measures_localobs_with_adaptive_gradflow(Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											Meas_Utils *meas_aux)
	{
	// meas no gradflow
	perform_measures_localobs(GC, geo, param, meas_aux);
	
	// meas gradflow
	if (param->d_agf_num_meas > 0)
		{
		int meas_count, accepted;
		double gftime, gftime_step;

		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		meas_count = 0;
		while(meas_count < param->d_agf_num_meas)
			{
			gradflow_RKstep_adaptive(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux);

			// if step accepted, perform measures
			if (accepted == 1 && fabs(gftime - param->d_agf_meas_each*(meas_count+1)) - param->d_agf_time_bin < MIN_VALUE )
				{
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}

			// adapt step to the time of next if this would be skipped
			if ((gftime + gftime_step - param->d_agf_meas_each*(meas_count+1)) > param->d_agf_time_bin )
				{
				gftime_step = param->d_agf_meas_each*(meas_count+1) - gftime;
				}
			}

		// restore gauge conf before gradflow from GC->lattice_copy
		restore_gauge_conf(GC, param);
		
		// print meas gradflow
		print_measures_arrays(param->d_agf_num_meas, GC->update_index, param, meas_aux);
		}
	
	// cold topcharge used to evaluate the multicanonical potential
	#ifdef MULTICANONICAL_MODE
	fprintf(meas_aux->datafilep, "% 18.12e ", GC->stored_topcharge);
	#endif
	
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	if (param->d_topcharge_tcorr_meas == 1) fflush(meas_aux->topchar_tcorr_filep);
	if (param->d_chi_prime_meas       == 1) fflush(meas_aux->chiprimefilep);
	}


void perform_measures_localobs_with_adaptive_gradflow_debug(Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											Meas_Utils *meas_aux,
											FILE *step_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int gradflowrepeat;
	
	// meas no gradflow
	perform_measures_localobs(GC, geo, param, meas_aux);
	
	// meas gradflow
	gradflowrepeat = (int)floor((param->d_agf_length+MIN_VALUE)/param->d_agf_meas_each);
	if (gradflowrepeat > 0)
		{
		int meas_count, accepted;
		double gftime, gftime_step, total_error;

		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		meas_count = 0;
		total_error = 0.0;
		fprintf(step_filep, "%ld % 18.12e % 18.12e % 18.12e\n", GC->update_index, gftime, gftime_step, total_error);
		while(meas_count < gradflowrepeat)
			{
			gradflow_RKstep_adaptive_debug(GC, geo, param, &gftime, &gftime_step, &accepted, &total_error, meas_aux);
			fprintf(step_filep, "%ld % 18.12e % 18.12e % 18.12e\n", GC->update_index, gftime, gftime_step, total_error);
			// step accepted, perform measures
			if (accepted == 1 && fabs(gftime - param->d_agf_meas_each*(meas_count+1)) - param->d_agf_time_bin < MIN_VALUE)
				{
				perform_measures_aux(GC, geo, param, meas_count, meas_aux);
				meas_count = meas_count + 1;
				}
			}
		//fprintf(step_filep, "\n");
		restore_gauge_conf(GC, param);

		// print meas gradflow
		print_measures_arrays(gradflowrepeat, GC->update_index, param, meas_aux);
		}
	
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	fflush(step_filep);
	if (param->d_topcharge_tcorr_meas == 1) fflush(meas_aux->topchar_tcorr_filep);
	if (param->d_chi_prime_meas       == 1) fflush(meas_aux->chiprimefilep);
	
	#else
	perform_measures_localobs_notopo(GC, geo, param, meas_aux);
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	#endif
	}


void perform_measures_localobs_with_adaptive_gradflow_debug2(Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											Meas_Utils *meas_aux,
											FILE *step_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	// meas no gradflow
	perform_measures_localobs(GC, geo, param, meas_aux);
	
	// meas gradflow
	if (param->d_agf_length > 0.0)
		{
		Gauge_Conf conf_accepted, conf_accepted_old;
		int accepted;
		double gftime, gftime_step, total_error, gftime_step_accepted;		
		
		// allocate memory
		init_gauge_conf_from_gauge_conf(&conf_accepted, GC, param);
		init_gauge_conf_from_gauge_conf(&conf_accepted_old, GC, param);

		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		total_error = 0.0;
		fprintf(step_filep, "%ld % 18.12e % 18.12e % 18.12e\n", GC->update_index, gftime, gftime_step, total_error);
		fflush(step_filep);
		while(gftime < param->d_agf_length)
			{
			gradflow_RKstep_adaptive_debug2(GC, geo, param, &gftime, &gftime_step, &accepted, &total_error, meas_aux);
			fprintf(step_filep, "%ld % 18.12e % 18.12e % 18.12e\n", GC->update_index, gftime, gftime_step, total_error);
			if (accepted == 1) //&& fabs(gftime - param->d_agf_meas_each*(meas_count+1)) - param->d_agf_time_bin < MIN_VALUE)
				{
				// save adaptive gradflow status
				equal_gauge_conf(&conf_accepted, GC, param);
				equal_gauge_conf(&conf_accepted_old, GC, param);
				gftime_step_accepted = gftime_step;
				
				// keep reducing gftime_step without advancing
				while (gftime_step > 1.01*param->d_agf_meas_each)
					{
					gftime_step = gftime_step-param->d_agf_meas_each;
					equal_gauge_conf(GC, &conf_accepted_old, param);
					gradflow_RKstep_adaptive_debug(GC, geo, param, &gftime, &gftime_step, &accepted, &total_error, meas_aux);
					fprintf(step_filep, "%ld % 18.12e % 18.12e % 18.12e\n", GC->update_index, gftime, gftime_step, total_error);
					}
				
				// restore adaptive gradflow status
				gftime = gftime + gftime_step_accepted;
				gftime_step = param->d_agf_step;
				equal_gauge_conf(GC, &conf_accepted, param);
				
				// meas gradflow
				fprintf(meas_aux->datafilep, "% 18.12e ", gftime);
				perform_measures_aux(GC, geo, param, 0, meas_aux);
				}
			else
				{
				if (gftime_step < 1.01*param->d_agf_meas_each) gftime_step = gftime_step/2.0;
				else gftime_step = gftime_step-param->d_agf_meas_each;
				}
			}
		restore_gauge_conf(GC, param);
		// free memory
		free_gauge_conf(&conf_accepted, param);
		free_gauge_conf(&conf_accepted_old, param);
		}
	
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	fflush(step_filep);
	if (param->d_topcharge_tcorr_meas == 1) fflush(meas_aux->topchar_tcorr_filep);
	if (param->d_chi_prime_meas       == 1) fflush(meas_aux->chiprimefilep);
	
	#else
	perform_measures_localobs_notopo(GC, geo, param, meas_aux);
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	#endif
	}


// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf * const GC,
												Geometry const * const geo,
												GParam const * const param,
												Meas_Utils *meas_aux)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	
	int i;
	double plaqs, plaqt, charge_nocooling, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1
	
	// meas no cooling
	plaquette(GC, geo, param, &plaqs, &plaqt);
	polyakov_with_tracedef(GC, geo, param, polyre, polyim);
	charge_nocooling=topcharge(GC, geo, param);

	// print meas no cooling
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", plaqs, plaqt);
	for(i=0; i<(int)floor(NCOLOR/2); i++) fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", polyre[i], polyim[i]);
	fprintf(meas_aux->datafilep, "% 18.12e ", charge_nocooling);
	
	// meas cooling
	topcharge_cooling(GC, geo, param, meas_aux);
	for(i=0; i<param->d_coolrepeat; i++)
		fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", meas_aux->charge[i], meas_aux->meanplaq[i]);
	
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	
	#else
	
	double plaqs, plaqt, polyre, polyim;
	
	plaquette(GC, geo, param, &plaqs, &plaqt);
	polyakov(GC, geo, param, &polyre, &polyim);
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e % 18.12e % 18.12e ", plaqs, plaqt, polyre, polyim);
	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	
	#endif
	}


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf * const GC,
								Geometry const * const geo,
								GParam const * const param,
								FILE *datafilep)
	{
	const int max_hit=50;
	const int dir=1;

	int i, mh, t_tmp;
	long r, r1, r2;
	double complex poly_corr;
	double poly_corr_abs, poly_corr_fluct, diff_sec;
	double complex *poly_array;
	time_t time1, time2;
	GAUGE_GROUP matrix, tmp;
	
	allocate_array_double_complex(&poly_array, param->d_space_vol, __FILE__, __LINE__);
	
	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0);
	#endif

	fprintf(datafilep, "Multihit optimization: \n");
	fprintf(datafilep, "the smaller the value the better the multihit\n");

	for(mh=1; mh<max_hit; mh++)
		{
		time(&time1);

		// polyakov loop computation
		for(r=0; r<param->d_space_vol; r++)
			{
			one(&matrix);
			for(i=0; i<param->d_size[0]; i++)
				{
				multihit(GC,
							geo,
							param,
							sisp_and_t_to_si(geo, r, i),
							0,
							mh,
							&tmp);
				times_equal(&matrix, &tmp);
				}
			poly_array[r]=retr(&matrix)+I*imtr(&matrix);
			}

		// average correlator computation
		poly_corr=0.0+I*0.0;
		poly_corr_abs=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			r1=sisp_and_t_to_si(geo, r, 0);
			for(i=0; i<param->d_dist_poly; i++)
				{
				r1=nnp(geo, r1, dir);
				}
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

			poly_corr += poly_array[r]*conj(poly_array[r2]);
			poly_corr_abs += cabs(poly_array[r]*conj(poly_array[r2]));
			}
		poly_corr*=param->d_inv_space_vol;
		poly_corr_abs*=param->d_inv_space_vol;

		// fluctuation of the average correlator computation
		poly_corr_fluct=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			r1=sisp_and_t_to_si(geo, r, 0);
			for(i=0; i<param->d_dist_poly; i++)
				{
				r1=nnp(geo, r1, dir);
				}
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
			poly_corr_fluct+=cabs( poly_array[r]*conj(poly_array[r2]) - poly_corr );
			}
		poly_corr_fluct*=param->d_inv_space_vol;


		time(&time2);
		diff_sec = difftime(time2, time1);

		fprintf(datafilep, "%d	% 18.12e	% 18.12e (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

		fflush(datafilep);
		}

	free(poly_array);
	}


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									FILE *datafilep)
	{
	int i;
	long r;
	double complex poly_corr;
	double poly_corr_abs, poly_corr_fluct;
	double complex *poly_array;
	
	allocate_array_double_complex(&poly_array, param->d_space_vol, __FILE__, __LINE__);
	
	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	multilevel_polycorr(GC,
								geo,
								param,
								param->d_size[0]);
	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
			}
		}

	// averages
	poly_corr=0.0+I*0.0;
	poly_corr_abs=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

		poly_corr+=poly_array[r];
		poly_corr_abs+=cabs(poly_array[r]);
		}
	poly_corr*=param->d_inv_space_vol;
	poly_corr_abs*=param->d_inv_space_vol;

	// fluctuations
	poly_corr_fluct=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_corr_fluct += cabs(poly_array[r]-poly_corr);
		}
	poly_corr_fluct*=param->d_inv_space_vol;

	// normalizations
	for(i=0; i<NLEVELS; i++)
		{
		poly_corr_abs*= sqrt(param->d_ml_upd[i]);
		poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs*=sqrt(param->d_multihit);
	poly_corr_fluct*=sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e ", poly_corr_abs);
	for(i=0; i<NLEVELS; i++)
		{
		fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
		}
	fprintf(datafilep, "(1, %d) ", param->d_multihit);
	fprintf(datafilep, "\n");
	fflush(datafilep);

	free(poly_array);
	}


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf * const GC,
								Geometry const * const geo,
								GParam const * const param,
								Meas_Utils *meas_aux)
	{
	#ifndef OPT_MULTIHIT
	#ifndef OPT_MULTILEVEL
		double ris;
		long r;
		int i;

		multilevel_polycorr(GC,
					 geo,
					 param,
					 param->d_size[0]);

		for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_space_vol; r++)
				{
				times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
				}
			}

		ris=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
			}
		ris*=param->d_inv_space_vol;

		fprintf(meas_aux->datafilep, "% 18.12e\n", ris);
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
void optimize_multihit_polycorradj(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									FILE *datafilep)
	{
	const int max_hit=50;
	const int dir=1;

	int i, mh, t_tmp;
	long r, r1, r2;
	double poly_corr, poly_corr_abs, poly_corr_fluct, diff_sec;
	double complex tr;
	double *poly_array;
	time_t time1, time2;
	GAUGE_GROUP matrix, tmp;
	
	allocate_array_double(&poly_array, param->d_space_vol, __FILE__, __LINE__);
	
	#ifdef THETA_MODE
	compute_clovers(GC, geo, param, 0);
	#endif

	fprintf(datafilep, "Multihit optimization: \n");
	fprintf(datafilep, "the smaller the value the better the multihit\n");

	for(mh=1; mh<max_hit; mh++)
		{
		time(&time1);

		// polyakov loop in the adjoint representation computation
		for(r=0; r<param->d_space_vol; r++)
			{
			one(&matrix);
			for(i=0; i<param->d_size[0]; i++)
				{
				multihit(GC,
							geo,
							param,
							sisp_and_t_to_si(geo, r, i),
							0,
							mh,
							&tmp);
				times_equal(&matrix, &tmp);
				}

			//trace of the matix in the fundamental representation
			tr=NCOLOR*(retr(&matrix)+I*imtr(&matrix));

			//trace of the matrix in adjoint representation
			poly_array[r]=cabs(tr*conj(tr))-1.0;
			#if NCOLOR != 1
			 poly_array[r]/=(NCOLOR*NCOLOR-1);
			#endif
			}

		// average correlator computation
		poly_corr=0.0;
		poly_corr_abs=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			r1=sisp_and_t_to_si(geo, r, 0);
			for(i=0; i<param->d_dist_poly; i++)
				{
				r1=nnp(geo, r1, dir);
				}
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

			poly_corr+= poly_array[r]*poly_array[r2];
			poly_corr_abs+=fabs(poly_array[r]*poly_array[r2]);
			}
		poly_corr*=param->d_inv_space_vol;
		poly_corr_abs*=param->d_inv_space_vol;

		// fluctuation of the average correlator computation
		poly_corr_fluct=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			r1=sisp_and_t_to_si(geo, r, 0);
			for(i=0; i<param->d_dist_poly; i++)
				{
				r1=nnp(geo, r1, dir);
				}
			si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

			poly_corr_fluct+=fabs(poly_array[r]*poly_array[r2]-poly_corr);
			}
		poly_corr_fluct*=param->d_inv_space_vol;

		time(&time2);
		diff_sec = difftime(time2, time1);

		fprintf(datafilep, "%d	% 18.12e	% 18.12e (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

		fflush(datafilep);
		}

	free(poly_array);
	}


// to optimize the multilevel (adjoint representation)
void optimize_multilevel_polycorradj(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									FILE *datafilep)
	{
	int i;
	long r;
	double poly_corr, poly_corr_abs, poly_corr_fluct;
	double *poly_array;
	
	allocate_array_double(&poly_array, param->d_space_vol, __FILE__, __LINE__);
	
	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	multilevel_polycorradj(GC,
									geo,
									param,
									param->d_size[0]);

	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
			}
		}

	// averages
	poly_corr=0.0;
	poly_corr_abs=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_array[r]=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));

		poly_corr+=poly_array[r];
		poly_corr_abs+=fabs(poly_array[r]);
		}
	poly_corr*=param->d_inv_space_vol;
	poly_corr_abs*=param->d_inv_space_vol;

	// fluctuations
	poly_corr_fluct=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_corr_fluct+=fabs(poly_array[r]-poly_corr);
		}
	poly_corr_fluct*=param->d_inv_space_vol;

	// normalizations
	for(i=0; i<NLEVELS; i++)
		{
		poly_corr_abs*=sqrt(param->d_ml_upd[i]);
		poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs*=sqrt(param->d_multihit);
	poly_corr_fluct*=sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e % 18.12e ", poly_corr_abs, poly_corr_fluct);
	for(i=0; i<NLEVELS; i++)
		{
		fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
		}
	fprintf(datafilep, "(1, %d) ", param->d_multihit);
	fprintf(datafilep, "\n");
	fflush(datafilep);

	free(poly_array);
	}


// perform the computation of the polyakov loop correlator in the adjoint representation with the multilevel algorithm
void perform_measures_polycorradj(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	#ifndef OPT_MULTIHIT
	#ifndef OPT_MULTILEVEL
		double ris;
		long r;
		int i;

		multilevel_polycorradj(GC, geo, param, param->d_size[0]);

		for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_space_vol; r++)
				{
				times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
				}
			}

		ris=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			ris+=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));
			}
		ris*=param->d_inv_space_vol;

		fprintf(meas_aux->datafilep, "% 18.12e\n", ris);
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
void optimize_multilevel_polycorr_long(Gauge_Conf * const GC,
										GParam const * const param,
										FILE *datafilep)
	{
	int i;
	long r;
	double poly_corr_abs, poly_corr_fluct;
	double complex poly_corr;
	double complex *poly_array;
	
	allocate_array_double_complex(&poly_array, param->d_space_vol, __FILE__, __LINE__);
	
	fprintf(datafilep, "Multilevel optimization: ");
	fprintf(datafilep, "the smaller the value the better the update\n");

	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
			}
		}

	// average
	poly_corr=0.0+I*0.0;
	poly_corr_abs=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

		poly_corr+=poly_array[r];
		poly_corr_abs+=cabs(poly_array[r]);
		}
	poly_corr*=param->d_inv_space_vol;
	poly_corr_abs*=param->d_inv_space_vol;

	// fluctuation
	poly_corr_fluct=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		poly_corr_fluct+=cabs(poly_array[r]-poly_corr);
		}
	poly_corr_fluct*=param->d_inv_space_vol;

	// normalization
	for(i=0; i<NLEVELS; i++)
		{
		poly_corr_abs*=sqrt(param->d_ml_upd[i]);
		poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
		}
	poly_corr_abs*=sqrt(param->d_ml_level0_repeat);
	poly_corr_fluct*=sqrt(param->d_ml_level0_repeat);

	poly_corr_abs*=sqrt(param->d_multihit);
	poly_corr_fluct*=sqrt(param->d_multihit);

	fprintf(datafilep, "% 18.12e % 18.12e ", poly_corr_abs, poly_corr_fluct);
	for(i=0; i<NLEVELS; i++)
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
void perform_measures_polycorr_long(Gauge_Conf * const GC,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	#ifdef OPT_MULTILEVEL
		optimize_multilevel_polycorr_long(GC, param, meas_aux->datafilep);
	#else
		double ris;
		long r;
		int i;

		for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_space_vol; r++)
				{
				times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
				}
			}

		ris=0.0;
		for(r=0; r<param->d_space_vol; r++)
			{
			ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
			}
		ris*=param->d_inv_space_vol;

		fprintf(meas_aux->datafilep, "% 18.12e\n", ris);
		fflush(meas_aux->datafilep);
	#endif
	}


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf * const GC,
								Geometry const * const geo,
								GParam const * const param,
								Meas_Utils *meas_aux)
	{
	double risr, risi;
	long r;
	int i;

	multilevel_tube_disc(GC, geo, param, param->d_size[0]);

	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
			}
		}

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
		risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf * const GC,
								Geometry const * const geo,
								GParam const * const param,
								Meas_Utils *meas_aux)
	{
	double risr, risi;
	long r;
	int i;

	multilevel_tube_conn(GC, geo, param, param->d_size[0]);

	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
			times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
			}
		}

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
		risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf * const GC,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	double risr, risi;
	long r;
	int i;

	for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_space_vol; r++)
			{
			times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
			times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
			times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
			}
		}

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
		risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
		risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	risr=0.0;
	risi=0.0;
	for(r=0; r<param->d_space_vol; r++)
		{
		risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
		}
	risr*=param->d_inv_space_vol;
	risi*=param->d_inv_space_vol;
	fprintf(meas_aux->datafilep, "% 18.12e % 18.12e ", risr, risi);

	fprintf(meas_aux->datafilep, "\n");
	fflush(meas_aux->datafilep);
	}


void allocate_measures_arrays(int const num_meas, GParam const * const param, double **meanplaq,
								double **clover_energy, double **charge, double **sum_q_timeslices,
								double **chi_prime, double ***charge_prime)
	{
	if (param->d_plaquette_meas == 1)
		allocate_array_double(meanplaq, num_meas, __FILE__, __LINE__);
	
	if (param->d_clover_energy_meas == 1)
		allocate_array_double(clover_energy, num_meas, __FILE__, __LINE__);	
	
	if (param->d_charge_meas == 1)
		allocate_array_double(charge, num_meas, __FILE__, __LINE__);
	
	if (param->d_topcharge_tcorr_meas == 1 )
		allocate_array_double(sum_q_timeslices, param->d_size[0], __FILE__, __LINE__);
	
	if (param->d_chi_prime_meas == 1)
		allocate_array_double(chi_prime, num_meas, __FILE__, __LINE__);
	
	if (param->d_charge_prime_meas == 1)
		{
		allocate_array_double_pointer(charge_prime, num_meas, __FILE__, __LINE__);
		for(int i=0; i<num_meas; i++)
			allocate_array_double(&((*charge_prime)[i]), STDIM, __FILE__, __LINE__);
		}
	}


void allocate_measures_arrays_cooling(int const num_meas, GParam const * const param,
										double **const meanplaq, double **const charge, double **const chi_prime)
	{
	allocate_array_double(charge, num_meas, __FILE__, __LINE__);
	allocate_array_double(meanplaq, num_meas, __FILE__, __LINE__);
	if (param->d_chi_prime_meas == 1)
		allocate_array_double(chi_prime, num_meas, __FILE__, __LINE__);
	}

// open data files
void open_data_file(Meas_Utils *meas_aux, char const * const data, char const * const chiprime,
						char const * const topchar_tcorr, GParam const * const param)
	{
	int i;

	if(param->d_start==2)
		{
		// open std data file (plaquette, polyakov, topological charge)
		meas_aux->datafilep=fopen(data, "r");
		if(meas_aux->datafilep!=NULL) // file exists
			{
			fclose(meas_aux->datafilep);
			meas_aux->datafilep=fopen(data, "a");
			}
		else
			{
 			meas_aux->datafilep=fopen(data, "w");
			print_header_datafile(meas_aux->datafilep, param);
			}
		
		// open chi prime data file
		if (param->d_chi_prime_meas == 1)
			{
			meas_aux->chiprimefilep=fopen(chiprime, "r");
			if(meas_aux->chiprimefilep!=NULL) // file exists
				{
				fclose(meas_aux->chiprimefilep);
				meas_aux->chiprimefilep=fopen(chiprime, "a");
				}
			else
				{
 				meas_aux->chiprimefilep=fopen(chiprime, "w");
				fprintf(meas_aux->chiprimefilep, "# %d ", STDIM);
				for(i=0; i<STDIM; i++) fprintf(meas_aux->chiprimefilep, "%d ", param->d_size[i]);
				fprintf(meas_aux->chiprimefilep, "\n");
				}
			}
		
		// open topocharge_tcorr data file
		if (param->d_topcharge_tcorr_meas == 1)
			{
			meas_aux->topchar_tcorr_filep=fopen(topchar_tcorr, "r");
			if(meas_aux->topchar_tcorr_filep!=NULL) // file exists
				{
				fclose(meas_aux->topchar_tcorr_filep);
				meas_aux->topchar_tcorr_filep=fopen(topchar_tcorr, "a");
				}
			else
				{
 				meas_aux->topchar_tcorr_filep=fopen(topchar_tcorr, "w");
				fprintf(meas_aux->topchar_tcorr_filep, "# %d ", STDIM);
				for(i=0; i<STDIM; i++) fprintf(meas_aux->topchar_tcorr_filep, "%d ", param->d_size[i]);
				fprintf(meas_aux->topchar_tcorr_filep, "\n");
				}
			}
		}
	else
		{
		// open std data file
		meas_aux->datafilep=fopen(data, "w");
		print_header_datafile(meas_aux->datafilep, param);
		
		// open chi prime data file
		if (param->d_chi_prime_meas == 1)
			{
			meas_aux->chiprimefilep=fopen(chiprime, "w");
			fprintf(meas_aux->chiprimefilep, "# %d ", STDIM);
			for(i=0; i<STDIM; i++) fprintf(meas_aux->chiprimefilep, "%d ", param->d_size[i]);
			fprintf(meas_aux->chiprimefilep, "\n");
			}
		
		// open topocharge_tcorr data file
		if (param->d_topcharge_tcorr_meas == 1)
			{
			meas_aux->topchar_tcorr_filep=fopen(topchar_tcorr, "w");
			fprintf(meas_aux->topchar_tcorr_filep, "# %d ", STDIM);
			for(i=0; i<STDIM; i++) fprintf(meas_aux->topchar_tcorr_filep, "%d ", param->d_size[i]);
			fprintf(meas_aux->topchar_tcorr_filep, "\n");
			}
		}
	
	fflush(meas_aux->datafilep);
	if (param->d_chi_prime_meas == 1 ) fflush(meas_aux->chiprimefilep);
	if (param->d_topcharge_tcorr_meas == 1 ) fflush(meas_aux->topchar_tcorr_filep);
	}


void init_meas_utils(Meas_Utils *meas_aux, GParam const * const param, int const replica_index)
	{
	int num_meas;
	char data_filename[STD_STRING_LENGTH], chiprime_filename[STD_STRING_LENGTH], topchar_tcorr_filename[STD_STRING_LENGTH];
	
	// max number of measures neeeded using any smoothing method
	num_meas = param->d_agf_num_meas;
	if (num_meas < param->d_gf_num_meas)
		num_meas = param->d_gf_num_meas;
	if (num_meas < param->d_coolrepeat) 
		num_meas = param->d_coolrepeat;
	
	if (num_meas > 0)
		{
		// allocate meas arrays
		if (param->d_plaquette_meas == 1)
			allocate_array_double(&(meas_aux->meanplaq), num_meas, __FILE__, __LINE__);
		
		if (param->d_clover_energy_meas == 1)
			allocate_array_double(&(meas_aux->clover_energy), num_meas, __FILE__, __LINE__);	
		
		if (param->d_charge_meas == 1)
			allocate_array_double(&(meas_aux->charge), num_meas, __FILE__, __LINE__);
		
		if (param->d_topcharge_tcorr_meas == 1)
			allocate_array_double(&(meas_aux->sum_q_timeslices), param->d_size[0], __FILE__, __LINE__);
		
		if (param->d_chi_prime_meas == 1)
			allocate_array_double(&(meas_aux->chi_prime), num_meas, __FILE__, __LINE__);
		
		if (param->d_charge_prime_meas == 1)
			{
			allocate_array_double_pointer(&(meas_aux->charge_prime), num_meas, __FILE__, __LINE__);
			for(int i=0; i<num_meas; i++)
				allocate_array_double(&(meas_aux->charge_prime[i]), STDIM, __FILE__, __LINE__);
			}
		
		// allocate auxiliary lattices
		for (int i=0; i<4; i++)
			{
			allocate_array_GAUGE_GROUP_pointer(&(meas_aux->lattice_aux[i]), param->d_volume, __FILE__, __LINE__);
			// TO DO: is parallelization ok? is it useful?
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r=0; r<(param->d_volume); r++)
				{
				allocate_array_GAUGE_GROUP(&(meas_aux->lattice_aux[i][r]), STDIM, __FILE__, __LINE__);
				}
			}
		}
	
	// open data files
	strcpy(data_filename, param->d_data_file);
	strcpy(chiprime_filename, param->d_chiprime_file);
	strcpy(topchar_tcorr_filename, param->d_topcharge_tcorr_file);
		
	#ifdef REPLICA_MEAS_MODE
	if (param->d_N_replica_pt > 1)
		{
		char aux[STD_STRING_LENGTH];
		sprintf(aux, "_replica_%d", replica_index);
		strcat(data_filename, aux);
		strcat(chiprime_filename, aux);
		strcat(topchar_tcorr_filename, aux);		
		}
	open_data_file(meas_aux, data_filename, chiprime_filename, topchar_tcorr_filename, param);
	#else
	if (replica_index == 0)
		{
		open_data_file(meas_aux, data_filename, chiprime_filename, topchar_tcorr_filename, param);
		}
	#endif
	}

void init_meas_utils_replica(Meas_Utils **meas_aux, GParam const * const param)
	{
	allocate_array_Meas_Utils(meas_aux, param->d_N_replica_pt, __FILE__, __LINE__);
	
	// init meas utils and data files for physical replica
	init_meas_utils(&((*meas_aux)[0]), param, 0);
	
	// init meas utils for other replicas if using replica meas mode
	#ifdef REPLICA_MEAS_MODE
	for (int i=1; i<param->d_N_replica_pt; i++)
		init_meas_utils(&((*meas_aux)[i]), param, i);
	#endif
	}

void free_meas_utils(Meas_Utils meas_aux, GParam const * const param, int const replica_index)
	{
	int num_meas;
	num_meas = param->d_agf_num_meas;
	if (num_meas < param->d_gf_num_meas)
		num_meas = param->d_gf_num_meas;
	if (num_meas < param->d_coolrepeat) 
		num_meas = param->d_coolrepeat;

	if (num_meas > 0)
		{
		// free meas arrays
		if (param->d_plaquette_meas == 1)
			free(meas_aux.meanplaq);
		
		if (param->d_clover_energy_meas == 1)
			free(meas_aux.clover_energy);	
		
		if (param->d_charge_meas == 1)
			free(meas_aux.charge);
		
		if (param->d_topcharge_tcorr_meas == 1 )
			free(meas_aux.sum_q_timeslices);
		
		if (param->d_chi_prime_meas == 1)
			free(meas_aux.chi_prime);
		
		if (param->d_charge_prime_meas == 1)
			{
			for(int i=0; i<num_meas; i++)
				free(meas_aux.charge_prime[i]);
			free(meas_aux.charge_prime);
			}
		
		// free auxiliary lattices
		for (int i=0; i<4; i++)
			{
			// TO DO: is parallelization ok? is it useful?
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS)
			#endif
			for(long r=0; r<(param->d_volume); r++)
				{
				free(meas_aux.lattice_aux[i][r]);
				}
			free(meas_aux.lattice_aux[i]);
			}
		
		// close data files
		#ifdef REPLICA_MEAS_MODE
		(void) replica_index;
		fclose(meas_aux.datafilep);
		if (param->d_chi_prime_meas==1) fclose(meas_aux.chiprimefilep);
		if (param->d_topcharge_tcorr_meas==1) fclose(meas_aux.topchar_tcorr_filep);
		#else
		if (replica_index == 0)
			{
			fclose(meas_aux.datafilep);
			if (param->d_chi_prime_meas==1) fclose(meas_aux.chiprimefilep);
			if (param->d_topcharge_tcorr_meas==1) fclose(meas_aux.topchar_tcorr_filep);
			}
		#endif
		}
	}

void free_meas_utils_replica(Meas_Utils *meas_aux, GParam const * const param)
	{
	// free meas utils for physical replica
	free_meas_utils(meas_aux[0], param, 0);
	
	// free meas utils for other replicas if using replica meas mode
	#ifdef REPLICA_MEAS_MODE
	for (int i=1; i<param->d_N_replica_pt; i++)
		free_meas_utils(meas_aux[i], param, i);
	#endif
	
	free(meas_aux);
	}
	
void print_measures_arrays(int const num_meas, long const update_index, GParam const * const param,
							Meas_Utils *meas_aux)
	{
	double time_step;
	if (param->d_agf_meas_each > 0.0) time_step = (param->d_agf_meas_each);
	else time_step = param->d_ngfsteps;
	
	for(int i=0; i<num_meas; i++)
		{
		if (param->d_plaquette_meas     == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->meanplaq[i]);
		if (param->d_clover_energy_meas == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->clover_energy[i]);
		if (param->d_charge_meas        == 1) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->charge[i]);
		if (param->d_chi_prime_meas     == 1) fprintf(meas_aux->chiprimefilep, "%ld % 18.12e % 18.12e\n", update_index, (i+1)*time_step, meas_aux->chi_prime[i]);
		if (param->d_charge_prime_meas  == 1) for (int j=0; j<STDIM; j++) fprintf(meas_aux->datafilep, "% 18.12e ", meas_aux->charge_prime[i][j]);
		}
	}


void free_measures_arrays(int const num_meas, GParam const * const param, double *meanplaq, double *clover_energy,
							double *charge, double *sum_q_timeslices,
							double *chi_prime, double **charge_prime)
	{
	if (num_meas > 0)
		{
		if (param->d_plaquette_meas       == 1) free(meanplaq);
		if (param->d_clover_energy_meas   == 1) free(clover_energy);
		if (param->d_charge_meas          == 1) free(charge);
		if (param->d_topcharge_tcorr_meas == 1) free(sum_q_timeslices);
		if (param->d_chi_prime_meas       == 1) free(chi_prime);
		if (param->d_charge_prime_meas    == 1)
			{
			for(int i=0; i<num_meas; i++) free(charge_prime[i]);
			free(charge_prime);
			}
		}
	}

#endif
