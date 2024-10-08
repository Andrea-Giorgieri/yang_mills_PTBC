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
void calcstaples_wilson(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						long r,
						int i,
						GAUGE_GROUP *M)
	{
	int j, l;
	long k;
	double complex factor;
	GAUGE_GROUP link1, link2, link3, link12, stap;

	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#else
	(void) param; // just to avoid warnings
	#endif

	zero(M); // M=0

	for(l=i+1; l< i + STDIM; l++)
	{
		j = (l % STDIM);

//
//		i ^
//		|	(1)
//		+----->-----+
//		|			|
//					|
//		|			V (2)
//					|
//		|			|
//		+-----<-----+-->	j
//		r	(3)
//

		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));				// link3 = (3)

		times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);	// stap=link12*stap^{dag}
	
		//twist (clockwise plaquette) modification
		factor=GC->Z[r][dirs_to_si(i,j)];		//Z_\mu\nu(x)
	
		times_equal_complex(&stap, factor); // Z_\mu\nu(x) * staple
	
		plus_equal(M, &stap);

//
//		i ^
//		|	(1)
//		|----<------+
//		|			|
//		|
//	(2) V			|
//		|
//		|			|
//		+------>----+--->j
//		k	(3)		r
//

		k=nnm(geo, r, j);

		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));				// link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));				// link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3);		// stap=link12*link3
		
		//twist (anticlockwise plaquette) modification
		factor=GC->Z[k][dirs_to_si(j,i)];	//Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
	
		times_equal_complex(&stap, factor); // Z_\mu\nu(x-\nu) * staple

		plus_equal(M, &stap);
	}
	}

// compute the staple for the trace deformed theory:
// in practice a Polyakov loop without a link
void calcstaples_tracedef(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						long r,
						int i,
						GAUGE_GROUP * M)
	{
	if(i!=0)
	{
	zero(M);
	#ifdef DEBUG
	fprintf(stderr, "Using calcstaples_tracedef for a non-temporal link (%s, %d)\n", __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	#endif
	}
	else
	{
	int j;
	long int rnext;
	GAUGE_GROUP aux;

	one(&aux);

	rnext=r;
	for(j=1; j<param->d_size[0]; j++)
		{
		rnext=nnp(geo, rnext, 0);
		times_equal(&aux, &(GC->lattice[rnext][0]));
		}

	equal(M, &aux);
	}
	}


// compute all the clovers in directions ortogonal to "dir"
void compute_clovers(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					int dir)
	{
	long r;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
	{
	GAUGE_GROUP aux;
	int i, j;

	for(i=0; i<STDIM; i++)
		{
		for(j=i+1; j<STDIM; j++)
			{
			if(i!=dir && j!=dir)
			{
			clover(GC, geo, param, r, i, j, &aux);

			equal(&(GC->clover_array[r][i][j]), &aux);
			minus_equal_dag(&(GC->clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

			equal(&(GC->clover_array[r][j][i]), &(GC->clover_array[r][i][j]));
			times_equal_real(&(GC->clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
			}
			}
		}
	}
	}


// compute the staple in position r, direction i and save it in M
// when an imaginary theta term is present
void calcstaples_with_topo(Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							long r,
							int i,
							GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#endif

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, M);
	#else

	if(STDIM!=4)
		{
		fprintf(stderr, "Error: imaginary theta term can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, topo_M, aux;
	const double coeff=(param->d_theta)*((double) NCOLOR)/(param->d_beta*128.0*PI*PI);
	long k;
	int j, l;
	int i0, j0;
	int sood1[4][4], sood2[4][4]; // signed ordered orthogonal directions
	double complex factor;

	zero(M);			// M=0
	zero(&topo_M);		// topo_M=0

	// the theta term is written as
	// theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
	// the independent permutations are 0123 0231 0312

	sood1[0][1] = 2;
	sood2[0][1] = 3;
	sood1[1][0] = 3;
	sood2[1][0] = 2;

	sood1[0][2] = 3;
	sood2[0][2] = 1;
	sood1[2][0] = 1;
	sood2[2][0] = 3;

	sood1[0][3] = 1;
	sood2[0][3] = 2;
	sood1[3][0] = 2;
	sood2[3][0] = 1;

	sood1[1][2] = 0;
	sood2[1][2] = 3;
	sood1[2][1] = 3;
	sood2[2][1] = 0;

	sood1[1][3] = 2;
	sood2[1][3] = 0;
	sood1[3][1] = 0;
	sood2[3][1] = 2;

	sood1[2][3] = 0;
	sood2[2][3] = 1;
	sood1[3][2] = 1;
	sood2[3][2] = 0;

	for(l=i+1; l< i + STDIM; l++)
		{
		j = (l % STDIM);

		i0=sood1[i][j];
		j0=sood2[i][j];

	//
	//		i ^
	//		|	(1)
	//	(b) +----->-----+ (c)
	//		|			|
	//					|
	//		|			V (2)
	//					|
	//		|			|
	//	(a) +-----<-----+-->	j
	//		r	(3)	(d)
	//
		
		zero(&topo_stap);	// topo_stap=0
		
		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));	// link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));				// link3 = (3)

		times_dag2(&link12, &link1, &link2);	// link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);		// stap=link12*stap^{dag}

		// clover insertion in (a)
		times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
		plus_equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);  // clover*stap
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link1, &(GC->clover_array[nnp(geo, nnp(geo, r, i), j)][i0][j0]));  // link1*clover
		times_equal_dag(&aux, &link2);		// *=link2^{dag}
		times_equal_dag(&aux, &link3);		// *=link3^{dag}
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times(&aux, &link12, &(GC->clover_array[nnp(geo, r, j)][i0][j0]));  // link1*link2*quadri
		times_equal_dag(&aux, &link3);		// *=link3^{dag}
		plus_equal(&topo_stap, &aux);
			
		//twist modification (clockwise plaquette)
		factor=GC->Z[r][dirs_to_si(i,j)];			//Z_\mu\nu(x)
		times_equal_complex(&topo_stap, factor);	//Z_\mu\nu(x) * staple
		times_equal_complex(&stap, factor);			//Z_\mu\nu(x) * staple
		
		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);

	//
	//		i ^
	//		|	(1)
	//	(d) +----<------+ (a)
	//		|			|
	//		|
	//	(2) V			|
	//		|
	//		|			| (b)
	//	(c) +------>----+--->j
	//		k	(3)	r
	//

		k=nnm(geo, r, j);
		zero(&topo_stap);	// topo_stap=0

		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));				// link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));				// link3 = (3)

		times_dag12(&link12, &link1, &link2);	// link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3);			// stap=link12*link3

		// clover insertion in (a)
		times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);	// clover*stap
		minus_equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &stap, &(GC->clover_array[r][i0][j0]));	// stap*clover
		minus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link12, &(GC->clover_array[k][i0][j0]));	// link1^{dag}*link2^{dag}*clover
		times_equal(&aux, &link3);								// *=link3
		minus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times_dag1(&aux, &link1, &(GC->clover_array[nnp(geo, k, i)][i0][j0]));  // link1^{dag}*clover
		times_equal_dag(&aux, &link2);			// *=link2^{dag}
		times_equal(&aux, &link3);				// *=link3
		minus_equal(&topo_stap, &aux);
			
		// twist modification (anticlockwise plaquette)
		factor=GC->Z[k][dirs_to_si(j,i)];			//Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		times_equal_complex(&topo_stap, factor);	//Z_\nu\mu(x-\nu) * staple
		times_equal_complex(&stap, factor);			//Z_\nu\mu(x-\nu) * staple

		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);	
		}

	times_equal_real(&topo_M, coeff);
	plus_equal(M, &topo_M);
	
	#endif
	}
	
// update functions for parallel tempering on defect

// compute all the clovers in directions ortogonal to "dir" for all replica
void compute_clovers_replica(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					int dir)
	{
	long s;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s)
	#endif
	for(s=0; s<((param->d_N_replica_pt)*param->d_volume); s++)
	{
	GAUGE_GROUP aux;
	int i, j;
		
		// s = i_r * volume + r
		long r = s % (param->d_volume); // space-time index
		int i_r = (int) (  (s-r) / (param->d_volume) ) ; // replica index

	for(i=0; i<STDIM; i++)
		{
		for(j=i+1; j<STDIM; j++)
			{
			if(i!=dir && j!=dir)
			{
			clover(&(GC[i_r]), geo, param, r, i, j, &aux);

			equal(&(GC[i_r].clover_array[r][i][j]), &aux);
			minus_equal_dag(&(GC[i_r].clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

			equal(&(GC[i_r].clover_array[r][j][i]), &(GC[i_r].clover_array[r][i][j]));
			times_equal_real(&(GC[i_r].clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
			}
			}
		}
	}
	}

// compute all the clovers in directions ortogonal to "dir" for all replica on a given rectangle
void compute_clovers_replica_rect(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					int dir, Rectangle const * const clover_rectangle)
	{
	long s;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s)
	#endif
	for(s=0; s<((param->d_N_replica_pt)*(clover_rectangle->d_vol_rect)); s++)
	{
	GAUGE_GROUP aux;
	int i, j;
		
		// s = i_r * volume_rect + n
		long n = s % (clover_rectangle->d_vol_rect); // space-time index on the rectangle
		long r = clover_rectangle->rect_sites[n];	// space-time index on the lattice
		int i_r = (int) ( (s-n) / (clover_rectangle->d_vol_rect) ) ; // replica index

	for(i=0; i<STDIM; i++)
		{
		for(j=i+1; j<STDIM; j++)
			{
			if(i!=dir && j!=dir)
			{
			clover(&(GC[i_r]), geo, param, r, i, j, &aux);

			equal(&(GC[i_r].clover_array[r][i][j]), &aux);
			minus_equal_dag(&(GC[i_r].clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

			equal(&(GC[i_r].clover_array[r][j][i]), &(GC[i_r].clover_array[r][i][j]));
			times_equal_real(&(GC[i_r].clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
			}
			}
		}
	}
	}

// evaluate non-topo staples with defect and twist factors in position r and direction i and save it in M
void calcstaples_wilson_with_defect(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						long r,
						int i,
						GAUGE_GROUP *M)
	{
	int j, l;
	long k;
	double complex factor;
	GAUGE_GROUP link1, link2, link3, link12, stap;

	#ifdef DEBUG
		if(r >= param->d_volume)
		{
			fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if(i >= STDIM)
		{
			fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	#else
		(void) param; // just to avoid warnings
	#endif

	zero(M); // M=0

	for(l=i+1; l< i + STDIM; l++)
	{
		j = (l % STDIM);

//
//	i	^
//		|	(1)
//		+----->-----+
//		|			|
//					|
//		|			V (2)
//					|
//		|			|
//		+-----<-----+-->	j
//		r	(3)
//

		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));				// link3 = (3)

		times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);	// stap=link12*stap^{dag}
		
		// boundary condition and twist (clockwise plaquette) modification
		factor=(GC->C[r][i])*(GC->C[nnp(geo, r, i)][j])*(GC->C[nnp(geo, r, j)][i])*(GC->C[r][j]); //K_\mu\nu(x)
		factor*=GC->Z[r][dirs_to_si(i,j)];		//Z_\mu\nu(x)

		times_equal_complex(&stap, factor); //K_\mu\nu(x) * Z_\mu\nu(x) * staple

		plus_equal(M, &stap);

//
//	i	^
//		|	(1)
//		|----<------+
//		|			|
//		|
//	(2) V			|
//		|
//		|			|
//		+------>----+--->j
//		k	(3)	r
//

		k=nnm(geo, r, j);

		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));				// link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));				// link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3);		// stap=link12*link3

		// boundary condition and twist (anticlockwise plaquette) modification
		factor=(GC->C[k][i])*(GC->C[nnp(geo, k, i)][j])*(GC->C[nnp(geo, k, j)][i])*(GC->C[k][j]); // K_\mu\nu(x-\nu)
		factor*=GC->Z[k][dirs_to_si(j,i)];	//Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))

		times_equal_complex(&stap, factor); //K_\mu\nu(x-\nu) * Z_\mu\nu(x-\nu) * staple

		plus_equal(M, &stap);
	}
	}

// evauate topo staples with defect on the non-topo contributions
// in psition r and direction i and save it in M
void calcstaples_with_topo_with_defect(Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							long r,
							int i,
							GAUGE_GROUP *M)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#endif

	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, M);
	#else

	if(STDIM!=4)
		{
		fprintf(stderr, "Error: imaginary theta term can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, topo_M, aux;
	const double coeff=(param->d_theta)*((double) NCOLOR)/(param->d_beta*128.0*PI*PI);
	long k;
	int j, l;
	int i0, j0;
	int sood1[4][4], sood2[4][4]; // signed ordered orthogonal directions
	double complex factor;

	zero(M);			// M=0
	zero(&topo_M);		// topo_M=0

	// the theta term is written as
	// theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
	// the independent permutations are 0123 0231 0312

	sood1[0][1] = 2;
	sood2[0][1] = 3;
	sood1[1][0] = 3;
	sood2[1][0] = 2;

	sood1[0][2] = 3;
	sood2[0][2] = 1;
	sood1[2][0] = 1;
	sood2[2][0] = 3;

	sood1[0][3] = 1;
	sood2[0][3] = 2;
	sood1[3][0] = 2;
	sood2[3][0] = 1;

	sood1[1][2] = 0;
	sood2[1][2] = 3;
	sood1[2][1] = 3;
	sood2[2][1] = 0;

	sood1[1][3] = 2;
	sood2[1][3] = 0;
	sood1[3][1] = 0;
	sood2[3][1] = 2;

	sood1[2][3] = 0;
	sood2[2][3] = 1;
	sood1[3][2] = 1;
	sood2[3][2] = 0;

	for(l=i+1; l< i + STDIM; l++)
		{
		j = (l % STDIM);

		i0=sood1[i][j];
		j0=sood2[i][j];

	//
	//		i ^
	//		|	(1)
	//	(b) +----->-----+ (c)
	//		|			|
	//					|
	//		|			V (2)
	//					|
	//		|			|
	//	(a) +-----<-----+-->	j
	//		r	(3)	(d)
	//
		
		zero(&topo_stap);	// topo_stap=0
		
		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));	// link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));				// link3 = (3)

		times_dag2(&link12, &link1, &link2);	// link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);		// stap=link12*stap^{dag}

		// clover insertion in (a)
		times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
		plus_equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);  // clover*stap
		plus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link1, &(GC->clover_array[nnp(geo, nnp(geo, r, i), j)][i0][j0]));  // link1*clover
		times_equal_dag(&aux, &link2);		// *=link2^{dag}
		times_equal_dag(&aux, &link3);		// *=link3^{dag}
		plus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times(&aux, &link12, &(GC->clover_array[nnp(geo, r, j)][i0][j0]));  // link1*link2*quadri
		times_equal_dag(&aux, &link3);		// *=link3^{dag}
		plus_equal(&topo_stap, &aux);
			
		// boundary condition modification (only affects non-topo staple) and twist (clockwise plaquette)
		factor=GC->Z[r][dirs_to_si(i,j)];			//Z_\mu\nu(x)
		times_equal_complex(&topo_stap, factor);	//Z_\mu\nu(x) * staple
		
		factor*=(GC->C[r][i])*(GC->C[nnp(geo, r, i)][j])*(GC->C[nnp(geo, r, j)][i])*(GC->C[r][j]); //K_\mu\nu(x)
		times_equal_complex(&stap, factor);			//K_\mu\nu(x) * Z_\mu\nu(x) * staple
		
		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);

	//
	//		i ^
	//		|	(1)
	//	(d) +----<------+ (a)
	//		|			|
	//		|
	//	(2) V			|
	//		|
	//		|			| (b)
	//	(c) +------>----+--->j
	//		k	(3)	r
	//

		k=nnm(geo, r, j);
		zero(&topo_stap);	// topo_stap=0

		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));				// link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));				// link3 = (3)

		times_dag12(&link12, &link1, &link2);	// link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3);			// stap=link12*link3

		// clover insertion in (a)
		times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);	// clover*stap
		minus_equal(&topo_stap, &aux);

		// clover insertion in (b)
		times(&aux, &stap, &(GC->clover_array[r][i0][j0]));	// stap*clover
		minus_equal(&topo_stap, &aux);

		// clover insertion in (c)
		times(&aux, &link12, &(GC->clover_array[k][i0][j0]));	// link1^{dag}*link2^{dag}*clover
		times_equal(&aux, &link3);								// *=link3
		minus_equal(&topo_stap, &aux);

		// clover insertion in (d)
		times_dag1(&aux, &link1, &(GC->clover_array[nnp(geo, k, i)][i0][j0]));  // link1^{dag}*clover
		times_equal_dag(&aux, &link2);			// *=link2^{dag}
		times_equal(&aux, &link3);				// *=link3
		minus_equal(&topo_stap, &aux);
			
		// boundary condition modification (only affects non-topo staple) and twist (anticlockwise plaquette)
		factor=GC->Z[k][dirs_to_si(j,i)];			//Z_\nu\mu(x-\nu) = conj(Z_\mu\nu(x-\nu))
		times_equal_complex(&topo_stap, factor);	//K_\mu\nu(x-\nu) * Z_\nu\mu(x-\nu) * staple
		
		factor*=(GC->C[k][i])*(GC->C[nnp(geo, k, i)][j])*(GC->C[nnp(geo, k, j)][i])*(GC->C[k][j]); // K_\mu\nu(x-\nu)
		times_equal_complex(&stap, factor);			//K_\mu\nu(x-\nu) * Z_\nu\mu(x-\nu) * staple

		plus_equal(M, &stap);
		plus_equal(&topo_M, &topo_stap);	
		}

	times_equal_real(&topo_M, coeff);
	plus_equal(M, &topo_M);
	
	#endif
	}


// perform an update with heatbath
void heatbath(Gauge_Conf *GC,
			Geometry const * const geo,
			GParam const * const param,
			long r,
			int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
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
void overrelaxation(Gauge_Conf *GC,
					Geometry const * const geo,
					GParam const * const param,
					long r,
					int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#endif
	(void) param; // just to avoid wornings

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap);
	#endif

	single_overrelaxation(&(GC->lattice[r][i]), &stap);
	}


// perform an update with metropolis
// return 1 if the proposed update is accepted
int metropolis(Gauge_Conf *GC,
				Geometry const * const geo,
				GParam const * const param,
				long r,
				int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#endif

	GAUGE_GROUP stap, new_link, tmp_matrix, rnd_matrix;
	double action_new, action_old;
	int acc;

	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap);
	#endif

	// compute old action
	times(&tmp_matrix, &(GC->lattice[r][i]), &stap);
	action_old=param->d_beta*(1.0-retr(&tmp_matrix));

	// compute the new link
	one(&tmp_matrix);
	rand_matrix(&rnd_matrix);
	times_equal_real(&rnd_matrix, param->d_epsilon_metro);
	plus_equal(&rnd_matrix, &tmp_matrix);
	unitarize(&rnd_matrix);	// rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
	if(casuale()<0.5)
	{
	times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
	}
	else
	{
	times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
	}

	// new action
	times(&tmp_matrix, &new_link, &stap);
	action_new=param->d_beta*(1.0-retr(&tmp_matrix));

	if(casuale()< exp(action_old-action_new))
	{
	equal(&(GC->lattice[r][i]), &new_link);
	acc=1;
	}
	else
	{
	acc=0;
	}

	return acc;
	}


// perform an update with metropolis with trace deformations
// return 1 if the proposed update is accepted
int metropolis_with_tracedef(Gauge_Conf *GC,
							Geometry const * const geo,
							GParam const * const param,
							long r,
							int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#endif

	GAUGE_GROUP stap_w, stap_td, new_link, tmp_matrix, rnd_matrix, poly;
	double action_new, action_old;
	double rpart, ipart;
	int j, acc;

	// compute old action
	#ifndef THETA_MODE
	calcstaples_wilson(GC, geo, param, r, i, &stap_w);
	#else
	calcstaples_with_topo(GC, geo, param, r, i, &stap_w);
	#endif
	times(&tmp_matrix, &(GC->lattice[r][i]), &stap_w);
	action_old=param->d_beta*(1.0-retr(&tmp_matrix));
	if(i==0) // just if we are updating a temporal link
	{
	// "staple" for trace deformation
	calcstaples_tracedef(GC, geo, param, r, i, &stap_td);

	// trace deformation contribution to action_old
	times(&poly, &(GC->lattice[r][i]), &stap_td);
	one(&tmp_matrix);
	for(j=0; j<(int)floor(NCOLOR/2.0); j++)
		{
		times_equal(&tmp_matrix, &poly);
		rpart=NCOLOR*retr(&tmp_matrix);
		ipart=NCOLOR*imtr(&tmp_matrix);
		action_old += param->d_h[j]*(rpart*rpart+ipart*ipart);
		}
	}

	// compute the update to be proposed
	one(&tmp_matrix);
	rand_matrix(&rnd_matrix);
	times_equal_real(&rnd_matrix, param->d_epsilon_metro);
	plus_equal(&rnd_matrix, &tmp_matrix);
	unitarize(&rnd_matrix);	// rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
	if(casuale()<0.5)
	{
	times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
	}
	else
	{
	times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
	}

	// compute the new action
	times(&tmp_matrix, &new_link, &stap_w);
	action_new=param->d_beta*(1.0-retr(&tmp_matrix));
	if(i==0) // just if we are updating a temporal link
	{
	// trace deformation contribution to action_new
	times(&poly, &new_link, &stap_td);
	one(&tmp_matrix);
	for(j=0; j<(int)floor(NCOLOR/2.0); j++)
		{
		times_equal(&tmp_matrix, &poly);
		rpart=NCOLOR*retr(&tmp_matrix);
		ipart=NCOLOR*imtr(&tmp_matrix);
		action_new += param->d_h[j]*(rpart*rpart+ipart*ipart);
		}
	}

	if(casuale()< exp(action_old-action_new))
	{
	equal(&(GC->lattice[r][i]), &new_link);
	acc=1;
	}
	else
	{
	acc=0;
	}

	return acc;
	}

	
// perform an update with heatbath in the presence of a defect
void heatbath_with_defect(Gauge_Conf *GC,
			Geometry const * const geo,
			GParam const * const param,
			long r,
			int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
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
void overrelaxation_with_defect(Gauge_Conf *GC,
					Geometry const * const geo,
					GParam const * const param,
					long r,
					int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
	fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
	fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#endif
	(void) param; // just to avoid wornings

	GAUGE_GROUP stap;

	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
	#endif

	single_overrelaxation(&(GC->lattice[r][i]), &stap);
	}
	

// perform a complete update
void update(Gauge_Conf * GC,
			Geometry const * const geo,
			GParam const * const param)
	{
	for(int i=0; i<STDIM; i++)
		{
		if(param->d_size[i]==1)
			{
			fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		}

	long r, num_even;
	int j, dir;
	
	/* Check if there's at least one even dimension of the lattice, i.e. check if d_volume is even.
	If there's at least one even dimension: d_volume/2 even sites and d_volume/2 odd sites.
	Otherwise: (d_volume+1)/2 even sites and (d_volume-1)/2 odd sites. */
	long is_even = ( param->d_volume ) % 2;
	num_even = ( param->d_volume + is_even ) / 2; // number of even sites
	
	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif 
		for(r=0; r<num_even; r++)
			{
			heatbath(GC, geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif 
		for(r=num_even; r<(param->d_volume); r++)
			{
			heatbath(GC, geo, param, r, dir);
			} 
		}

	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif 
			for(r=0; r<num_even; r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif 
			for(r=num_even; r<(param->d_volume); r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}
			}
		}
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r, dir)
	#endif 
	for(r=0; r<(param->d_volume); r++)
		{
		for(dir=0; dir<STDIM; dir++)
			{
			unitarize(&(GC->lattice[r][dir]));
			} 
		}

	GC->update_index++;
	}
	
// update all replica in the presence of a defect
void update_with_defect(Gauge_Conf * GC, Geometry const * const geo, GParam const * const param)
	{
	for(int i=0; i<STDIM; i++)
		{
		if(param->d_size[i]==1)
			{
			fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		}

	long s, num_even, num_odd;
	int j, dir;
	
	/* Check if there's at least one even dimension of the lattice, i.e. check if d_volume is even.
	If there's at least one even dimension: d_volume/2 even sites and d_volume/2 odd sites.
	Otherwise: (d_volume+1)/2 even sites and (d_volume-1)/2 odd sites. */
	long is_even = ( param->d_volume ) % 2;
	num_even = ( param->d_volume + is_even ) / 2; // number of even sites
	num_odd  = ( param->d_volume - is_even ) / 2; // number of odd sites

	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif
		for(s=0; s<((param->d_N_replica_pt)*num_even); s++)
			{
			// s = i * num_even + r
			long r = s % num_even; // site index
			int i = (int) ( (s-r) / num_even ); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<((param->d_N_replica_pt)*num_odd); s++)
			{
			// s = i * num_odd + aux ; aux = r - num_even
			long aux = s % num_odd;
			long r = num_even + aux; // site index
			int i = (int) ( (s-aux) / num_odd ); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		}

	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif

		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif
			for(s=0; s<((param->d_N_replica_pt)*num_even); s++)
				{
				// s = i * num_even + r
				long r = s % num_even; // site index
				int i = (int) ( (s-r) / num_even ); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif
			for(s=0; s<((param->d_N_replica_pt)*num_odd); s++)
				{
				// s = i * num_odd + aux ; aux = r - num_even
				long aux = s % num_odd;
				long r = num_even + aux; // site index
				int i = (int) ( (s-aux) / num_odd ); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<(param->d_N_replica_pt)*(param->d_volume); s++)
		{
		// s = i * volume + r
		long r = s % param->d_volume;
		int i = (int) ( (s-r) / (param->d_volume) );
		for(dir=0; dir<STDIM; dir++)
			{
			unitarize(&(GC[i].lattice[r][dir]));
			} 
		}
	}

// hierarchical update functions

// update all replica only on a given rectangle in the presence of a defect
void update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
											Rectangle const * const most_update, Rectangle const * const clover_rectangle)
	{
	long s, num_even, num_odd;
	int j,dir;
	
	#ifndef THETA_MODE
	(void) clover_rectangle; // to avoid compiler warning of unused variable
	#endif	
		
	/* Check if there's at least one even dimension of the rectangle, i.e. check if d_vol_rect is even.
		If there's at least one even dimension: d_vol_rect/2 even sites and d_vol_rect/2 odd sites.
		Otherwise: (d_vol_rect+1)/2 even sites and (d_vol_rect-1)/2 odd sites. */
	long is_even = ( most_update->d_vol_rect ) % 2;
	num_even = ( most_update->d_vol_rect + is_even ) / 2; // number of even sites
	num_odd  = ( most_update->d_vol_rect - is_even ) / 2; // number of odd sites

	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_even + n
			long n = s % num_even;					// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-n) / num_even );		// replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_odd + aux; aux = n - num_even
			long aux = s % num_odd;
			long n = aux + num_even;				// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-aux) / num_odd );	// replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		}

	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif

		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_even + n
				long n = s % num_even;					// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-n) / num_even );		// replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_odd + aux; aux = n - num_even
				long aux = s % num_odd; 
				long n = aux + num_even;				// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-aux) / num_odd );	// replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((most_update->d_vol_rect)*(param->d_N_replica_pt)); s++)
		{
		for(dir=0; dir<STDIM; dir++)
			{
			// s = i * volume_rect + n
			long n = s % (most_update->d_vol_rect);					// site index on rectangle
			long r = most_update->rect_sites[n];					// site index on lattice
			int i = (int) ( (s-n) / (most_update->d_vol_rect) );	// replica index
			unitarize(&(GC[i].lattice[r][dir]));
			} 
		}
	}

// perform a hierarchical update on all rectangles
void hierarchical_update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
														int const hierarc_level, 
														Rectangle const * const most_update,
														Rectangle const * const clover_rectangle,
														Rectangle const * const swap_rectangle,
														Acc_Utils *acc_counters)
	{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
			{
			update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]) );
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			}
		} // end if
	else
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
			{
			update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]) );	
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			hierarchical_update_rectangle_with_defect(GC,geo,param,(hierarc_level+1),most_update,clover_rectangle,swap_rectangle,acc_counters);
			}
		} // end else

	}
	
// perform a single step of parallel tempering with hierarchical update
void parallel_tempering_with_hierarchical_update(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
												Rectangle const * const most_update, Rectangle const * const clover_rectangle,
												Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)
	{
	int i;
	int start_hierarc=0; // first hierarc level is  0
	
	// full update + hierarchical update + swaps and translations after every sweep
	update_with_defect(GC, geo, param);
	if(param->d_N_replica_pt>1)
	{
		swap(GC, geo, param, swap_rectangle, acc_counters);
		conf_translation(&(GC[0]), geo, param);
		if(param->d_N_hierarc_levels>0)
			hierarchical_update_rectangle_with_defect(GC,geo,param,start_hierarc,most_update,clover_rectangle,swap_rectangle,acc_counters);
	}

	// increase update index of all replica
	for(i=0;i<param->d_N_replica_pt; i++)
		GC[i].update_index++;
	}

// perform a complete update with trace deformation
void update_with_trace_def(Gauge_Conf * GC,
							Geometry const * const geo,
							GParam const * const param,
							double *acc_td)
	{
	int *a;
	long r, asum, num_even, num_sp_even;
	int j, dir, t;
	
	allocate_array_int(&a, param->d_space_vol, __FILE__, __LINE__);
	for(r=0; r<param->d_space_vol; r++) a[r]=0;

	num_even    = (param->d_volume + (param->d_volume % 2)) / 2;
	num_sp_even = (param->d_space_vol + (param->d_space_vol % 2)) / 2;

	// heatbath on spatial links
	for(dir=1; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<num_even; r++)
			{
			heatbath(GC, geo, param, r, dir);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=num_even; r<(param->d_volume); r++)
			{
			heatbath(GC, geo, param, r, dir);
			}
		}

	// metropolis on temporal links
	for(t=0; t<param->d_size[0]; t++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, 0);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<num_sp_even; r++)
			{
			long r4=sisp_and_t_to_si(geo, r, t);
			a[r]+=metropolis_with_tracedef(GC, geo, param, r4, 0);
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=num_sp_even; r<(param->d_space_vol); r++)
			{
			long r4=sisp_and_t_to_si(geo, r, t);
			a[r]+=metropolis_with_tracedef(GC, geo, param, r4, 0);
			}
		}

	asum=0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for reduction(+:asum) private(r)
	#endif
	for(r=0; r<param->d_space_vol; r++)
		{
		asum+=(long)a[r];
		}

	*acc_td=((double)asum)*param->d_inv_vol;

	// overrelax spatial links
	for(dir=1; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers(GC, geo, param, dir);
		#endif

		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<num_even; r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=num_even; r<(param->d_volume); r++)
				{
				overrelaxation(GC, geo, param, r, dir);
				}
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r, dir)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		for(dir=0; dir<STDIM; dir++)
			{
			unitarize(&(GC->lattice[r][dir]));
			}
		}

	free(a);
	GC->update_index++;
	}


// perform n cooling steps minimizing the action at theta=0
void cooling(Gauge_Conf *GC,
			Geometry const * const geo,
			GParam const * const param,
			int n)
	{
	long r, num_even;
	int i, k;
	
	num_even = (param->d_volume + (param->d_volume % 2)) / 2;

	for(k=0; k<n; k++)
		{
		// cooling
		for(i=0; i<STDIM; i++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<num_even; r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson(GC, geo, param, r, i, &staple);
				cool(&(GC->lattice[r][i]), &staple);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=num_even; r<(param->d_volume); r++)
				{
				GAUGE_GROUP staple;
				calcstaples_wilson(GC, geo, param, r, i, &staple);
				cool(&(GC->lattice[r][i]), &staple);
				}
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r, i)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		for(i=0; i<STDIM; i++)
			{
			unitarize(&(GC->lattice[r][i]));
			}
		}
	}


// perform a single step of the Runge Kutta integrator for the Wilson flow
// as described in Luscher arXiv:1006.4518 app. C
void gradflow_RKstep(Gauge_Conf *GC,
					Gauge_Conf *helper1,
					Gauge_Conf *helper2,
					Geometry const * const geo,
					GParam const *const param,
					double dt)
	{
	long r;
	int dir;

	// initialize
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
			equal(&(helper2->lattice[r][dir]), &(GC->lattice[r][dir]));
			}
		}


	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, link;

			calcstaples_wilson(helper1, geo, param, r, dir, &staple);
			equal(&link, &(helper1->lattice[r][dir]));
			times(&aux, &link, &staple);				// aux=link*staple
			times_equal_real(&aux, -dt/4.0);
			equal(&(helper2->lattice[r][dir]), &aux);	// helper2=aux
			taexp(&aux);
			times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
			}
		}

	// now helper1=W_0, helper2=(1/4)Z_0 and GC=W_1

	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, link;

			calcstaples_wilson(GC, geo, param, r, dir, &staple);
			equal(&link, &(GC->lattice[r][dir]));
			times(&aux, &link, &staple);				// aux=link*staple
			times_equal_real(&aux, -dt*8.0/9.0);
			minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/9.0); // 1/4 was in Z_0
			equal(&(helper2->lattice[r][dir]), &aux);
			taexp(&aux);
			times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
			}
		}

	// now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0 and GC=W_1

	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, link;

			calcstaples_wilson(helper1, geo, param, r, dir, &staple);
			equal(&link, &(helper1->lattice[r][dir]));
			times(&aux, &link, &staple);					// aux=link*staple
			times_equal_real(&aux, -dt*3.0/4.0);
			minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
			taexp(&aux);
			times(&(GC->lattice[r][dir]), &aux, &link);	// GC=aux*link
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		for(int i=0; i<STDIM; i++)
			{
			unitarize(&(GC->lattice[r][i]));
			}
		}
	}
	
// perform a single step of the Runge Kutta integrator for the Wilson flow
// with adaptive integration step as described in Fritzsch-Ramos arXiv:1301.4388 app. D
void gradflow_RKstep_adaptive(Gauge_Conf *GC,
					Gauge_Conf *GC_old,
					Gauge_Conf *helper1,
					Gauge_Conf *helper2,
					Gauge_Conf *helper3,
					Geometry const * const geo,
					GParam const *const param,
					double *t,
					double *dt,
					int *accepted)
	{
	long r;
	int dir, j;
	double max_dist, *local_max_dist;
	
	// initialize
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			equal(&(GC_old->lattice[r][dir]), &(GC->lattice[r][dir]));
			equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
			}
		}

	// now helper1 = GC_old = GC = W_0
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, link;

			calcstaples_wilson(helper1, geo, param, r, dir, &staple);
			equal(&link, &(helper1->lattice[r][dir]));
			times(&aux, &link, &staple);				// aux=link*staple
			times_equal_real(&aux, -(*dt));
			equal(&(helper2->lattice[r][dir]), &aux);	// helper2=aux
			times_equal_real(&aux, 1.0/4.0);
			taexp(&aux);
			times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
			}
		}

	// now helper1=W_0, helper2=Z_0 and GC=W_1
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, aux2, link;

			calcstaples_wilson(GC, geo, param, r, dir, &staple);
			equal(&link, &(GC->lattice[r][dir]));
			times(&aux, &link, &staple);				// aux=link*staple
			times_equal_real(&aux, -(*dt));
			equal(&aux2, &aux);
			
			times_equal_real(&aux2, 2.0);
			minus_equal(&aux2, &(helper2->lattice[r][dir]));
			taexp(&aux2);
			times(&(helper3->lattice[r][dir]), &aux2, &(helper1->lattice[r][dir])); // helper3=aux2*helper1
			
			times_equal_real(&aux, 8.0/9.0);
			minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/36.0);
			equal(&(helper2->lattice[r][dir]), &aux);
			taexp(&aux);
			times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
			}
		}

	// now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0, helper3=W'_2, and GC=W_1
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++)
			{
			GAUGE_GROUP staple, aux, link;

			calcstaples_wilson(helper1, geo, param, r, dir, &staple);
			equal(&link, &(helper1->lattice[r][dir]));
			times(&aux, &link, &staple);					// aux=link*staple
			times_equal_real(&aux, -(*dt)*3.0/4.0);
			minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
			taexp(&aux);
			times(&(GC->lattice[r][dir]), &aux, &link);	// GC=aux*link
			}
		}
	// now helper3 = W'_2 and GC = W_3
	
	// error calculation and final unitarization
	allocate_array_double(&local_max_dist, NTHREADS, __FILE__, __LINE__);
	for (j=0; j<NTHREADS; j++) local_max_dist[j] = 0.0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		int i, thread_num;
		double dist;
		thread_num = 0;
		#ifdef OPENMP_MODE
		thread_num = omp_get_thread_num();
		#endif
		for(i=0; i<STDIM; i++)
			{
			// Note: unitarize() after dist calculation, otherwise the integration
			//       error can be greatly over-estimated for some configurations
			minus_equal(&(helper3->lattice[r][i]), &(GC->lattice[r][i]));
			dist = norm(&(helper3->lattice[r][i]))/((double)NCOLOR*(double)NCOLOR);
			unitarize(&(GC->lattice[r][i]));
			if (dist > local_max_dist[thread_num]) local_max_dist[thread_num] = dist;
			}
		}
	max_dist = MIN_VALUE;
	for (j=0; j<NTHREADS; j++)
		{
		if (local_max_dist[j] > max_dist) max_dist = local_max_dist[j];
		}
	
	if (max_dist < param->d_agf_delta) //if the integration step is accepted, advance t
		{
		*accepted = 1;
		*t = *t + *dt;
		}
	else //if the integration step is rejected, reset gauge conf
		{
		*accepted = 0;
		for(dir=0; dir<STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_volume; r++)
				{
				equal(&(GC->lattice[r][dir]), &(GC_old->lattice[r][dir]));
				}
			}
		}

	// new integration step
	*dt = *dt * 0.95 * pow(param->d_agf_delta/max_dist, 1.0/3.0);
	free(local_max_dist);
	}

void gradflow_RKstep_adaptive_debug(Gauge_Conf *GC,
					Gauge_Conf *GC_old,
					Gauge_Conf *helper1,
					Gauge_Conf *helper2,
					Gauge_Conf *helper3,
					Geometry const * const geo,
					GParam const *const param,
					double *t,
					double *dt,
					int *accepted, double *total_error)
	{
	long r;
	int dir, j;
	double max_dist, *local_max_dist;
	//double mean_dist, clover1, clover2;
	
	// initialize
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		equal(&(GC_old->lattice[r][dir]), &(GC->lattice[r][dir]));
		equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
		}
	}
	// now helper1 = GC_old = GC = W_0
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, link;

		calcstaples_wilson(helper1, geo, param, r, dir, &staple);
		equal(&link, &(helper1->lattice[r][dir]));
		times(&aux, &link, &staple);				// aux=link*staple
		times_equal_real(&aux, -(*dt));
		equal(&(helper2->lattice[r][dir]), &aux);	// helper2=aux
		times_equal_real(&aux, 1.0/4.0);
		taexp(&aux);
		times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		unitarize(&(helper1->lattice[r][dir]));
		}
	}
	// now helper1=W_0, helper2=Z_0 and GC=W_1
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, aux2, link;

		calcstaples_wilson(GC, geo, param, r, dir, &staple);
		equal(&link, &(GC->lattice[r][dir]));
		times(&aux, &link, &staple);				// aux=link*staple
		times_equal_real(&aux, -(*dt));
		equal(&aux2, &aux);
		
		times_equal_real(&aux2, 2.0);
		minus_equal(&aux2, &(helper2->lattice[r][dir]));
		taexp(&aux2);
		times(&(helper3->lattice[r][dir]), &aux2, &(helper1->lattice[r][dir])); // helper3=aux2*helper1
		
		times_equal_real(&aux, 8.0/9.0);
		minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/36.0);
		equal(&(helper2->lattice[r][dir]), &aux);
		taexp(&aux);
		times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		unitarize(&(helper1->lattice[r][dir]));
		unitarize(&(helper3->lattice[r][dir]));
		}
	}
	// now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0, helper3=W'_2, and GC=W_1
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, link;

		calcstaples_wilson(helper1, geo, param, r, dir, &staple);
		equal(&link, &(helper1->lattice[r][dir]));
		times(&aux, &link, &staple);					// aux=link*staple
		times_equal_real(&aux, -(*dt)*3.0/4.0);
		minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
		taexp(&aux);
		times(&(GC->lattice[r][dir]), &aux, &link);	// GC=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		}
	}
	// now helper3 = W'_2 and GC = W_3
	
	// error calculation
	allocate_array_double(&local_max_dist, NTHREADS, __FILE__, __LINE__);
	
	for (j=0; j<NTHREADS; j++) local_max_dist[j] = 0.0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		int i, thread_num;
		double dist;
		thread_num = 0;
		#ifdef OPENMP_MODE
		thread_num = omp_get_thread_num();
		#endif
		for(i=0; i<STDIM; i++)
			{
			minus_equal(&(helper3->lattice[r][i]), &(GC->lattice[r][i]));
			dist = norm(&(helper3->lattice[r][i]))/((double)NCOLOR*(double)NCOLOR);
			if (dist > local_max_dist[thread_num]) local_max_dist[thread_num] = dist;
			}
		}
	max_dist = MIN_VALUE;
	for (j=0; j<NTHREADS; j++)
		{
		if (local_max_dist[j] > max_dist) max_dist = local_max_dist[j];
		}
	
	// error calculation debug
	//clover_disc_energy(GC, geo, param, &clover1);
	//clover_disc_energy(helper3, geo, param, &clover2);
	//mean_dist = fabs((clover1-clover2)/clover1);
	
	*total_error = max_dist;
	if (max_dist < param->d_agf_delta && *dt < 1.01*param->d_agf_meas_each) //if the integration step is accepted, advance t
		{
		*accepted = 1;
		*t = *t + *dt;
		//*total_error += max_dist;
		*dt = param->d_agf_step;
		}
	else //if the integration step is rejected, reset gauge conf
		{
		*accepted = 0;
		for(dir=0; dir<STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_volume; r++)
				{
				equal(&(GC->lattice[r][dir]), &(GC_old->lattice[r][dir]));
				}
			}
		*dt = *dt-param->d_agf_meas_each;
		}
	// new integration step
	//*dt = *dt * 0.95 * pow(param->d_agf_delta/max_dist, 1.0/3.0);
	}

void gradflow_RKstep_adaptive_debug2(Gauge_Conf *GC,
					Gauge_Conf *GC_old,
					Gauge_Conf *helper1,
					Gauge_Conf *helper2,
					Gauge_Conf *helper3,
					Geometry const * const geo,
					GParam const *const param,
					double *t,
					double *dt,
					int *accepted, double *total_error)
	{
	long r;
	int dir, j;
	double max_dist, *local_max_dist;
	//double mean_dist, clover1, clover2;
	// initialize
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		equal(&(GC_old->lattice[r][dir]), &(GC->lattice[r][dir]));
		equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
		}
	}
	// now helper1 = GC_old = GC = W_0
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, link;

		calcstaples_wilson(helper1, geo, param, r, dir, &staple);
		equal(&link, &(helper1->lattice[r][dir]));
		times(&aux, &link, &staple);				// aux=link*staple
		times_equal_real(&aux, -(*dt));
		equal(&(helper2->lattice[r][dir]), &aux);	// helper2=aux
		times_equal_real(&aux, 1.0/4.0);
		taexp(&aux);
		times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		unitarize(&(helper1->lattice[r][dir]));
		}
	}
	// now helper1=W_0, helper2=Z_0 and GC=W_1
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, aux2, link;

		calcstaples_wilson(GC, geo, param, r, dir, &staple);
		equal(&link, &(GC->lattice[r][dir]));
		times(&aux, &link, &staple);				// aux=link*staple
		times_equal_real(&aux, -(*dt));
		equal(&aux2, &aux);
		
		times_equal_real(&aux2, 2.0);
		minus_equal(&aux2, &(helper2->lattice[r][dir]));
		taexp(&aux2);
		times(&(helper3->lattice[r][dir]), &aux2, &(helper1->lattice[r][dir])); // helper3=aux2*helper1
		
		times_equal_real(&aux, 8.0/9.0);
		minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/36.0);
		equal(&(helper2->lattice[r][dir]), &aux);
		taexp(&aux);
		times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		unitarize(&(helper1->lattice[r][dir]));
		unitarize(&(helper3->lattice[r][dir]));
		}
	}
	// now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0, helper3=W'_2, and GC=W_1
	for(dir=0; dir<STDIM; dir++)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<param->d_volume; r++)
		{
		GAUGE_GROUP staple, aux, link;

		calcstaples_wilson(helper1, geo, param, r, dir, &staple);
		equal(&link, &(helper1->lattice[r][dir]));
		times(&aux, &link, &staple);					// aux=link*staple
		times_equal_real(&aux, -(*dt)*3.0/4.0);
		minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
		taexp(&aux);
		times(&(GC->lattice[r][dir]), &aux, &link);	// GC=aux*link
		
		unitarize(&(GC->lattice[r][dir]));
		}
	}
	// now helper3 = W'_2 and GC = W_3
	
	// error calculation
	allocate_array_double(&local_max_dist, NTHREADS, __FILE__, __LINE__);
	
	for (j=0; j<NTHREADS; j++) local_max_dist[j] = 0.0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r=0; r<(param->d_volume); r++)
		{
		int i, thread_num;
		double dist;
		thread_num = 0;
		#ifdef OPENMP_MODE
		thread_num = omp_get_thread_num();
		#endif
		for(i=0; i<STDIM; i++)
			{
			minus_equal(&(helper3->lattice[r][i]), &(GC->lattice[r][i]));
			dist = norm(&(helper3->lattice[r][i]))/((double)NCOLOR*(double)NCOLOR);
			if (dist > local_max_dist[thread_num]) local_max_dist[thread_num] = dist;
			}
		}
	max_dist = MIN_VALUE;
	for (j=0; j<NTHREADS; j++)
		{
		if (local_max_dist[j] > max_dist) max_dist = local_max_dist[j];
		}
	
	// error calculation debug
	//clover_disc_energy(GC, geo, param, &clover1);
	//clover_disc_energy(helper3, geo, param, &clover2);
	//mean_dist = fabs((clover1-clover2)/clover1);
	
	*total_error = max_dist;
	if (max_dist < param->d_agf_delta)// && *dt < 1.01*param->d_agf_meas_each) //if the integration step is accepted, advance t
		{
		*accepted = 1;
		(void)t; 
	//	*t = *t + *dt;
	//	//*total_error += max_dist;
	//	*dt = param->d_agf_step;
		}
	else //if the integration step is rejected, reset gauge conf
		{
		*accepted = 0;
		for(dir=0; dir<STDIM; dir++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(r)
			#endif
			for(r=0; r<param->d_volume; r++)
				{
				equal(&(GC->lattice[r][dir]), &(GC_old->lattice[r][dir]));
				}
			}
		//*dt = *dt-param->d_agf_meas_each;
		}
	// new integration step
	//*dt = *dt * 0.95 * pow(param->d_agf_delta/max_dist, 1.0/3.0);
	}
	
// n step of ape smearing with parameter alpha
// new=Proj[old + alpha *staple ]
void ape_smearing(Gauge_Conf *GC,
				Geometry const * const geo,
				GParam const *const param,
				double alpha,
				int n)
	{
	Gauge_Conf helper1;
	long r;
	int dir, count;

	init_gauge_conf_from_gauge_conf(&helper1, GC, param); //helper1=GC

	for(count=0; count<n; count++)
		{
		if(count%2==0) // smear(helper1)->GC
			{
			for(dir=0; dir<STDIM; dir++)
				{
				#ifdef OPENMP_MODE
				#pragma omp parallel for num_threads(NTHREADS) private(r)
				#endif
				for(r=0; r<param->d_volume; r++)
					{
					GAUGE_GROUP staple, link;
					
					calcstaples_wilson(&helper1, geo, param, r, dir, &staple);
					equal(&link, &(helper1.lattice[r][dir]));
					times_equal_real(&link, 1-alpha);
					times_equal_real(&staple, alpha/6.0);
					plus_equal_dag(&link, &staple);
					unitarize(&link);
					equal(&(GC->lattice[r][dir]), &link);
					}
				}
			}
		else // smear(GC)->helper1
			{
			for(dir=0; dir<STDIM; dir++)
				{
				#ifdef OPENMP_MODE
				#pragma omp parallel for num_threads(NTHREADS) private(r)
				#endif
				for(r=0; r<param->d_volume; r++)
					{
					GAUGE_GROUP staple, link;
					
					calcstaples_wilson(GC, geo, param, r, dir, &staple);
					equal(&link, &(GC->lattice[r][dir]));
					times_equal_real(&link, 1-alpha);
					times_equal_real(&staple, alpha/6.0);
					plus_equal_dag(&link, &staple);
					unitarize(&link);
					equal(&(helper1.lattice[r][dir]), &link);
					}
				}
			}
		}

	if(n>0 && n%2==0) equal_gauge_conf(GC, &helper1, param); // GC=helper1
	free_gauge_conf(&helper1, param);
	}




#endif
