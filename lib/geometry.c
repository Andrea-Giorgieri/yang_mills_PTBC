#ifndef GEOMETRY_C
#define GEOMETRY_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>

#include"../include/memalign.h"
#include"../include/geometry.h"
#include"../include/gparam.h"

#if STDIM == 4

// signed ordered orthogonal directions, only needed in 4D to write the theta term as
// theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood[\mu][\nu][0] sood[\mu][\nu][1]} )
// the independent permutations are 0123 0231 0312
const int g_signed_ord_orth_dir[4][4][2] = {
		[0][1] = {2, 3},
		[1][0] = {3, 2},

		[0][2] = {3, 1},
		[2][0] = {1, 3},

		[0][3] = {1, 2},
		[3][0] = {2, 1},

		[1][2] = {0, 3},
		[2][1] = {3, 0},

		[1][3] = {2, 0},
		[3][1] = {0, 2},

		[2][3] = {0, 1},
		[3][2] = {1, 0}
		};

// the i-th independent permutation among (0123), (0213) and (0312) is given by
// (ipd[0][i], ipd[1][i], ipd[2][i], ipd[3][i])
const int g_indep_perm_dir[4][3] = {
		[0] = {0, 0, 0},
		[1] = {1, 2, 3},
		[2] = {2, 1, 1},
		[3] = {3, 3, 2}
		};

#endif

// single index 4d = even/odd lexicographic index 4d
// single index 3d = even/odd lexicographic index 3d
void init_indexing_lexeo(void)
	{
	cart_to_si = &cart_to_lexeo;
	si_to_cart = &lexeo_to_cart;
	lex_to_si = &lex_to_lexeo;
	si_to_lex = &lexeo_to_lex;
	sisp_and_t_to_si_compute = &lexeosp_and_t_to_lexeo;
	si_to_sisp_and_t_compute = &lexeo_to_lexeosp_and_t;
	si_to_sisp_and_mu_compute = &lexeo_to_lexeosp_and_mu;

	// for rectangles
	cart_to_si_rect = &cart_to_lexeo_rect;
	}


// initialize geometry
void init_geometry(Geometry *geo, GParam const *const param)
	{
	// allocate memory
	allocate_array_long_pointer(&(geo->d_nnp), param->d_volume, __FILE__, __LINE__);
	allocate_array_long_pointer(&(geo->d_nnm), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_long(&(geo->d_nnp[r]), STDIM, __FILE__, __LINE__);
		allocate_array_long(&(geo->d_nnm[r]), STDIM, __FILE__, __LINE__);
		}

	allocate_array_int_pointer(&(geo->d_mucomp), STDIM, __FILE__, __LINE__);
	allocate_array_long_pointer(&(geo->d_muorth), STDIM, __FILE__, __LINE__);
	allocate_array_long_pointer_pointer(&(geo->d_mutimespace), STDIM, __FILE__, __LINE__);
	for(int i = 0; i < STDIM; i++)
		{
		allocate_array_int(&(geo->d_mucomp[i]), param->d_volume, __FILE__, __LINE__);
		allocate_array_long(&(geo->d_muorth[i]), param->d_volume, __FILE__, __LINE__);
		allocate_array_long_pointer(&(geo->d_mutimespace[i]), param->d_size[i], __FILE__, __LINE__);
		for(long r = 0; r < param->d_size[i]; r++)
			allocate_array_long(&(geo->d_mutimespace[i][r]), (param->d_volume) / (param->d_size[i]), __FILE__, __LINE__);
		}

	// initialize nearest neighbors
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		int cartcoord[STDIM];
		si_to_cart(cartcoord, r, param);

		for(int i = 0; i < STDIM; i++)
			{
			int const value = cartcoord[i];

			cartcoord[i] = (value + 1) % param->d_size[i];
			geo->d_nnp[r][i] = cart_to_si(cartcoord, param);

			cartcoord[i] = (value - 1 + param->d_size[i]) % param->d_size[i];
			geo->d_nnm[r][i] = cart_to_si(cartcoord, param);

			cartcoord[i] = value;
			}
		}

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		int x_mu;
		long r_sp;

		for(int mu = 0; mu < STDIM; mu++)
			{
			si_to_sisp_and_mu_compute(&r_sp, &x_mu, r, mu, param);
			geo->d_mutimespace[mu][x_mu][r_sp] = r;
			geo->d_mucomp[mu][r] = x_mu;
			geo->d_muorth[mu][r] = r_sp;
			}
		}

	for(int mu = 0; mu < STDIM; mu++)
		for(int i = 0; i < STDIM - 1; i++)
			geo->d_orth_dir[mu][i] = orthogonal_dir(mu, i);

	#ifdef DEBUG
	test_geometry(geo, param);
	#endif
	}

// free memory
void free_geometry(Geometry *geo, GParam const *const param)
	{
	for(long r = 0; r < param->d_volume; r++)
		{
		free(geo->d_nnp[r]);
		free(geo->d_nnm[r]);
		}
	free(geo->d_nnp);
	free(geo->d_nnm);

	for(int i = 0; i < STDIM; i++)
		{
		for(long r = 0; r < param->d_size[i]; r++)
			{
			free(geo->d_mutimespace[i][r]);
			}
		free(geo->d_mutimespace[i]);
		free(geo->d_muorth[i]);
		free(geo->d_mucomp[i]);
		}
	free(geo->d_mutimespace);
	free(geo->d_muorth);
	free(geo->d_mucomp);
	}


long nnp(Geometry const *const geo, long r, int i);


long nnm(Geometry const *const geo, long r, int i);


long sisp_and_t_to_si(Geometry const *const geo, long sisp, int t);


void si_to_sisp_and_t(long *sisp, int *t, Geometry const *const geo, long si);


void test_geometry(Geometry const *const geo, GParam const *const param)
	{
	int cart[STDIM], cartsp[STDIM - 1], t;

	// test of lex_to_cart <-> cart_to_lex
	for(long si = 0; si < param->d_volume; si++)
		{
		lex_to_cart(cart, si, param);
		long const test_res = cart_to_lex(cart, param);

		REQUIRE(si == test_res, "geometry test failed: %ld != %ld", si, test_res);
		}

	// test of lexeo_to_cart <-> cart_to_lexeo
	for(long si = 0; si < param->d_volume; si++)
		{
		lexeo_to_cart(cart, si, param);
		long const test_res = cart_to_lexeo(cart, param);

		REQUIRE(si == test_res, "geometry test failed: %ld != %ld", si, test_res);
		}

	// test of nnp <-> nnm
	for(long si = 0; si < param->d_volume; si++)
		{
		for(int dir = 0; dir < STDIM; dir++)
			{
			long const si_bis = nnp(geo, si, dir);
			long const test_res = nnm(geo, si_bis, dir);

			REQUIRE(si == test_res, "geometry test failed: %ld != %ld", si, test_res);
			}
		}

	// test of lexsp_to_cartsp <-> cartsp_to_lexsp
	for(long sisp = 0; sisp < param->d_space_vol[0]; sisp++)
		{
		lexsp_to_cartsp(cartsp, sisp, param);
		long const test_res = cartsp_to_lexsp(cartsp, param);

		REQUIRE(sisp == test_res, "geometry test failed: %ld != %ld", sisp, test_res);
		}

	// test of lexeosp_to_cartsp <-> cartsp_to_lexeosp
	for(long sisp = 0; sisp < param->d_space_vol[0]; sisp++)
		{
		lexeosp_to_cartsp(cartsp, sisp, param);
		long const test_res = cartsp_to_lexeosp(cartsp, param);

		REQUIRE(sisp == test_res, "geometry test failed: %ld != %ld", sisp, test_res);
		}

	// test of lexeosp_and_t_to_lexeo <-> lexeo_to_lexeosp_and_t
	for(long si = 0; si < param->d_volume; si++)
		{
		long sisp;
		lexeo_to_lexeosp_and_t(&sisp, &t, si, param);
		long const test_res = lexeosp_and_t_to_lexeo(sisp, t, param);

		REQUIRE(si == test_res, "geometry test failed: %ld != %ld", si, test_res);
		}
	}


//plane i-j -> single index, for twist factors
int dirs_to_si(int const i, int const j)
	{
	#ifdef DEBUG
	ASSERT(i != j, "directions cannot be equal (%d)", i);
	#endif
	if(i < j) return i * (2 * STDIM - 3 - i) / 2 + j - 1;                 //clockwise
	return j * (2 * STDIM - 3 - j) / 2 + i - 1 + STDIM * (STDIM - 1) / 2; //anticlockwise
	}


// check if lexeo index r is on the defect
int is_on_defect(long const r, GParam const *const param)
	{
	int cartcoord[STDIM];
	si_to_cart(cartcoord, r, param);
	for(int i = 0; i < STDIM; i++)
		{
		if(i == param->d_defect_dir)
			{
			if(cartcoord[i] != param->d_size[i] - 1) return 0;
			}
		else
			{
			int i_defect = i;
			if(i > param->d_defect_dir) i_defect -= 1;
			if(cartcoord[i] >= param->d_L_defect[i_defect]) return 0;
			}
		}
	return 1;
	}


//------------ these are not to be used outside geometry.c ----------------

// cartesian coordinates -> lexicographic index
long cart_to_lex(int const *const cartcoord, GParam const *const param)
	{
	long res = 0;
	long aux = 1;
	for(int i = 0; i < STDIM; i++)
		{
		res += cartcoord[i] * aux;
		aux *= param->d_size[i];
		}

	// res = cartcoord[0]
	//		+cartcoord[1]*size[0]
	//		+cartcoord[2]*size[0]*size[1]
	//		+...
	//		+cartcoord[STDIM-1]*size[0]*size[1]*...*size[STDIM-2]

	return res;

	/*
	int i;
	long res, aux;

	res=0;
	aux=1;
	for(i=STDIM-1; i>=0; i--)
		{
		res+=cartcoord[i]*aux;
		aux*=param->d_size[i];
		}

	// res = cartcoord[STDIM-1] +
	//		+cartcoord[STDIM-2]*size[STDIM-1]+
	//		+cartcoord[STDIM-3]*size[STDIM-1]*size[STDIM-2]+
	//		+...
	//		+cartcoord[0]*size[STDIM-1]*size[STDIM-2]*...*size[1]

	return res;
	*/
	}


// lexicographic index -> cartesian coordinates
void lex_to_cart(int *cartcoord, long lex, GParam const *const param)
	{
	long aux[STDIM];

	aux[0] = 1;
	for(int i = 1; i < STDIM; i++)
		{
		aux[i] = aux[i - 1] * param->d_size[i - 1];
		}
	// aux[0]=1
	// aux[1]=size[0]
	// aux[2]=size[0]*size[1]
	// ...
	// aux[STDIM-1]=size[0]*size[1]*...*size[STDIM-2]

	for(int i = STDIM - 1; i >= 0; i--)
		{
		cartcoord[i] = (int) (lex / aux[i]);
		lex -= aux[i] * cartcoord[i];
		}

	/*
	int i;
	long aux[STDIM];

	aux[STDIM-1]=1;
	for(i=STDIM-2; i>=0; i--)
		{
		aux[i]=aux[i+1]*param->d_size[i+1];
		}
	// aux[STDIM-1] = 1
	// aux[STDIM-2] = size[STDIM-1]
	// aux[STDIM-3] = size[STDIM-1]*size[STDIM-2]
	// ...
	// aux[0]		 = size[STDIM-1]*size[STDIM-2]*...*size[1]

	for(i=0; i<STDIM; i++)
		{
		cartcoord[i]=(int) (lex/aux[i]);
		lex-=aux[i]*cartcoord[i];
		}
	*/
	}


// cartesian coordinates -> lexicographic eo index
long cart_to_lexeo(int const *const cartcoord, GParam const *const param)
	{
	long const lex = cart_to_lex(cartcoord, param);

	int eo = 0;
	for(int i = 0; i < STDIM; i++)
		{
		eo += cartcoord[i];
		}

	if(eo % 2 == 0)
		{
		return lex / 2;
		}
	return (lex + param->d_volume) / 2;
	// even sites are written first
	}


// lexicographic eo index -> cartesian coordinates
void lexeo_to_cart(int *cartcoord, long lexeo, GParam const *const param)
	{
	long lex;

	if(param->d_volume % 2 == 0)
		{
		if(lexeo < param->d_volume / 2)
			{
			lex = 2 * lexeo;
			}
		else
			{
			lex = 2 * (lexeo - param->d_volume / 2);
			}
		lex_to_cart(cartcoord, lex, param);

		int eo = 0;
		for(int i = 0; i < STDIM; i++)
			{
			eo += cartcoord[i];
			}
		eo = eo % 2;

		// this is to take care of the case d_volume is even but not
		// all the lattice extents are even
		if((eo == 0 && lexeo >= param->d_volume / 2) ||
		   (eo == 1 && lexeo < param->d_volume / 2))
			{
			lex += 1;
			lex_to_cart(cartcoord, lex, param);
			}
		}
	else
		{
		if(lexeo <= param->d_volume / 2)
			{
			lex = 2 * lexeo;
			}
		else
			{
			lex = 2 * (lexeo - param->d_volume / 2) - 1;
			}
		lex_to_cart(cartcoord, lex, param);
		}
	}


//	lexicographic index -> lexicographic eo index
long lex_to_lexeo(long lex, GParam const *const param)
	{
	int cartcoord[STDIM];

	lex_to_cart(cartcoord, lex, param);

	return cart_to_lexeo(cartcoord, param);
	}


//	lexicographic eo index -> lexicographic index
long lexeo_to_lex(long lexeo, GParam const *const param)
	{
	int cartcoord[STDIM];

	lexeo_to_cart(cartcoord, lexeo, param);

	return cart_to_lex(cartcoord, param);
	}


// spatial cartesian coordinates -> spatial lexicographic index
long cartsp_to_lexsp(int const *const ccsp, GParam const *const param)
	{
	// the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
	// cc	= t x1 x2 ... x_{STDIM-1}
	// ccsp =	x1 x2		x_{STDIM-1}

	long res = 0;
	long aux = 1;
	for(int i = 0; i < STDIM - 1; i++)
		{
		res += ccsp[i] * aux;
		aux *= param->d_size[i + 1];
		}

	// res = ccsp[0]
	//		+ccsp[1]*size[1]
	//		+ccsp[2]*size[1]*size[2]
	//		+...
	//		+ccsp[STDIM-2]*size[1]*size[2]*...*size[STDIM-2]

	return res;

	/*
	int i;
	long res, aux;

	res=0;
	aux=1;
	for(i=STDIM-2; i>=0; i--)
		{
		res+=ccsp[i]*aux;
		aux*=param->d_size[i+1];
		}

	// res = ccsp[STDIM-2] +
	//		+ccsp[STDIM-3]*size[STDIM-1]+
	//		+ccsp[STDIM-4]*size[STDIM-1]*size[STDIM-2]+
	//		+ccsp
	//		+ccsp[0]*size[STDIM-1]*size[STDIM-2]*...*size[2]

	return res;
	*/
	}

// spatial cartesian coordinates -> spatial lexicographic index
long cartsp_to_lexsp_mu(int const *const ccsp, int mu, GParam const *const param)
	{
	// the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
	// cc	= x0 x1 x_mu ... x_{STDIM-1}
	// ccsp = x0 x1	 -   ... x_{STDIM-1}

	long res = 0;
	long aux = 1;
	for(int i = 0; i < STDIM - 1; i++)
		{
		res += ccsp[i] * aux;
		int const j = orthogonal_dir(mu, i);
		aux *= param->d_size[j];
		}

	// res = ccsp[0]
	//		+ccsp[1]*size[0]
	//		+ccsp[2]*size[0]*size[1]
	//		+... (skipping size[mu])
	//		+ccsp[STDIM-2]*size[1]*size[2]*...*size[STDIM-2]

	return res;
	}


// spatial lexicographic index -> spatial cartesian coordinates
void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const *const param)
	{
	// the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
	// cc	= t x1 x2 ... x_{STDIM-1}
	// ccsp =	x1 x2		x_{STDIM-1}

	int i;
	long aux[STDIM - 1];

	aux[0] = 1;
	for(i = 1; i < STDIM - 1; i++)
		{
		aux[i] = aux[i - 1] * param->d_size[i];
		}
	// aux[0]=1
	// aux[1]=size[1]
	// aux[2]=size[1]*size[2]
	// ...
	// aux[STDIM-2]=size[1]*size[2]*...*size[STDIM-2]

	for(i = STDIM - 2; i >= 0; i--)
		{
		ccsp[i] = (int) (lexsp / aux[i]);
		lexsp -= aux[i] * ccsp[i];
		}

	/*
	int i;
	long aux[STDIM-1];

	aux[STDIM-2]=1;
	for(i=STDIM-3; i>=0; i--)
		{
		aux[i]=aux[i+1]*param->d_size[i+2];
		}
	// aux[STDIM-2] = 1
	// aux[STDIM-3] = size[STDIM-1]
	// aux[STDIM-4] = size[STDIM-1]*size[STDIM-2]
	// ...
	// aux[0]		 = size[STDIM-1]*size[STDIM-2]*...*size[2]

	for(i=0; i<STDIM-1; i++)
		{
		ccsp[i]=(int) (lexsp/aux[i]);
		lexsp-=aux[i]*ccsp[i];
		}
	*/
	}


// spatial cartesian coordinates -> spatial lexicographic eo index
long cartsp_to_lexeosp(int const *const ccsp, GParam const *const param)
	{
	long const lexsp = cartsp_to_lexsp(ccsp, param);

	int eo = 0;
	for(int i = 0; i < STDIM - 1; i++)
		{
		eo += ccsp[i];
		}

	if(eo % 2 == 0)
		{
		return lexsp / 2;
		}
	return (lexsp + param->d_space_vol[0]) / 2;
	}

// spatial cartesian coordinates -> spatial lexicographic eo index
long cartsp_to_lexeosp_mu(int const *const ccsp, int mu, GParam const *const param)
	{
	long const lexsp = cartsp_to_lexsp_mu(ccsp, mu, param);

	int eo = 0;
	for(int i = 0; i < STDIM - 1; i++)
		{
		eo += ccsp[i];
		}

	if(eo % 2 == 0)
		{
		return lexsp / 2;
		}
	return (lexsp + param->d_space_vol[mu]) / 2;
	}


// spatial lexicographic eo index -> spatial cartesian coordinates
void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const *const param)
	{
	long lexsp;

	if(param->d_space_vol[0] % 2 == 0)
		{
		if(lexeosp < param->d_space_vol[0] / 2)
			{
			lexsp = 2 * lexeosp;
			}
		else
			{
			lexsp = 2 * (lexeosp - param->d_space_vol[0] / 2);
			}
		lexsp_to_cartsp(ccsp, lexsp, param);

		int eo = 0;
		for(int i = 0; i < STDIM - 1; i++)
			{
			eo += ccsp[i];
			}
		eo = eo % 2;

		if((eo == 0 && lexeosp >= param->d_space_vol[0] / 2) ||
		   (eo == 1 && lexeosp < param->d_space_vol[0] / 2))
			{
			lexsp += 1;
			lexsp_to_cartsp(ccsp, lexsp, param);
			}
		}
	else
		{
		if(lexeosp <= param->d_space_vol[0] / 2)
			{
			lexsp = 2 * lexeosp;
			}
		else
			{
			lexsp = 2 * (lexeosp - param->d_space_vol[0] / 2) - 1;
			}
		lexsp_to_cartsp(ccsp, lexsp, param);
		}
	}


// spatial lexicographic index -> spatial lexicographic eo index
long lexsp_to_lexeosp(long lexsp, GParam const *const param)
	{
	int ccsp[STDIM];

	lexsp_to_cartsp(ccsp, lexsp, param);

	return cartsp_to_lexeosp(ccsp, param);
	}


//	spatial lexicographic eo index -> spatial lexicographic index
long lexeosp_to_lexsp(long lexeosp, GParam const *const param)
	{
	int ccsp[STDIM];

	lexeosp_to_cartsp(ccsp, lexeosp, param);

	return cartsp_to_lexsp(ccsp, param);
	}


// lexicographic eo spatial and time -> lexicographic eo index
long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const *const param)
	{
	int cc[STDIM];

	lexeosp_to_cartsp(cc + 1, lexeosp, param);
	cc[0] = t;

	return cart_to_lexeo(cc, param);
	}


// lexicographic eo index -> lexicographic eo spatial and time
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const *const param)
	{
	int cc[STDIM], ccsp[STDIM - 1];

	lexeo_to_cart(cc, lexeo, param);

	*t = cc[0];

	for(int i = 0; i < STDIM - 1; i++)
		{
		ccsp[i] = cc[i + 1];
		}

	*lexeosp = cartsp_to_lexeosp(ccsp, param);
	}

void lexeo_to_lexeosp_and_mu(long *lexeosp, int *t, long lexeo, int mu, GParam const *const param)
	{
	int cc[STDIM], ccsp[STDIM - 1];

	lexeo_to_cart(cc, lexeo, param);

	*t = cc[mu];

	for(int i = 0; i < STDIM - 1; i++)
		ccsp[i] = cc[(i < mu) ? i : i + 1];
	*lexeosp = cartsp_to_lexeosp_mu(ccsp, mu, param);
	}

// geometry for rectangles to be used for hierarchical update during parallel tempering

// reduce a generic integer component in the interval [0,L_max-1]
int periodic_condition(int const coord, int const L_max)
	{
	return (coord % L_max + L_max) % L_max;
	}


// i-th direction orthogonal to mu
int orthogonal_dir(int const mu, int const i)
	{
	return (i < mu) ? i : i + 1;
	}


// cartesian -> lexicographic eo (on a given rectangle)
// this function should not be used outside geometry.c
long cart_to_lexeo_rect(int const *const cartcoord, Rectangle const *const most_update)
	{
	int eo = 0;
	long res = 0, aux = 1;

	for(int i = 0; i < STDIM; i++)
		{
		res += cartcoord[i] * aux;
		aux *= most_update->d_size_rect[i];
		eo += cartcoord[i];
		}
	eo = eo % 2;
	return (eo * (most_update->d_vol_rect) + res) / 2; // even sites first
	}


void init_rect(Rectangle *most_update, int const L_R, GParam const *const param)
	{
	// sizes of rectangle and ranges of rect coordinates
	int aux_L[STDIM], size_min[STDIM];

	// rectangle sizes and min of rect coordinates along every direction
	aux_L[param->d_defect_dir] = 2 * L_R + 1;
	size_min[param->d_defect_dir] = param->d_size[param->d_defect_dir] - L_R - 1;
	for(int i = 0; i < STDIM - 1; i++)
		{
		int const j = orthogonal_dir(param->d_defect_dir, i);
		aux_L[j] = 2 * L_R + param->d_L_defect[i];
		size_min[j] = -L_R;
		}

	// d_size_rect[i] must not exceed d_size[i]
	// if so, the exceeding dimension of the rectangle is just the respective dimension of the whole lattice
	// and the i-th coordinate just ranges from 0 to d_size[i]-1
	for(int i = 0; i < STDIM; i++)
		{
		if(aux_L[i] >= param->d_size[i])
			{
			aux_L[i] = param->d_size[i];
			size_min[i] = 0;
			}
		}

	// volume of rectangle
	long V = 1;
	for(int i = 0; i < STDIM; i++)
		V *= aux_L[i];

	// allocate rectangle
	allocate_array_long(&(most_update->rect_sites), V, __FILE__, __LINE__);

	// save dimensions and volume of the rectangle
	for(int i = 0; i < STDIM; i++)
		most_update->d_size_rect[i] = aux_L[i];
	most_update->d_vol_rect = V;

	// save lexeo index of sites of the rectangle
	int coord[STDIM];      // cartesian coordinates on the whole lattice
	int real_coord[STDIM]; // cartesian coordinates after using periodic conditions
	int coord_rect[STDIM]; // cartesian coordinates on the rectangle
	long r, r_rect;        // lexeo index on the rectangle and lexeo index on the whole lattice

	for(int i = 0; i < STDIM; i++)
		coord_rect[i] = 0;

	// loop over each rectangle site
	while(1)
		{
		// compute rectangle coordinates coord_rect and lattice coordinates real_coord for current iteration
		for(int i = 0; i < STDIM; i++)
			{
			coord[i] = size_min[i] + coord_rect[i];
			real_coord[i] = periodic_condition(coord[i], param->d_size[i]);
			}

		// convert to lexicographical index on rectangle and lattice
		r_rect = cart_to_si_rect(coord_rect, most_update);
		r = cart_to_si(real_coord, param);

		most_update->rect_sites[r_rect] = r;

		// move to the next rectangle site
		int i = STDIM - 1;
		while(i >= 0)
			{
			if(coord_rect[i] < aux_L[i] - 1)
				{
				coord_rect[i]++;
				break;
				}
			coord_rect[i] = 0;
			i--;
			}

		// break when all indices have reached their max values
		if(i < 0)
			break;
		}
	}

void free_rect(Rectangle *most_update)
	{
	free(most_update->rect_sites);
	}

void init_rect_hierarc(Rectangle **most_update, Rectangle **clover_rect, GParam const *const param)
	{
	if(param->d_N_hierarc_levels == 0)
		{
		most_update = NULL;
		(void) clover_rect;
		}
	else
		{
		int i;
		allocate_array_Rectangle(most_update, param->d_N_hierarc_levels, __FILE__, __LINE__);
		for(i = 0; i < param->d_N_hierarc_levels; i++)
			init_rect(&((*most_update)[i]), param->d_L_rect[i], param);

		#ifdef THETA_MODE
		allocate_array_Rectangle(clover_rect, param->d_N_hierarc_levels, __FILE__, __LINE__);
		for(i = 0; i < param->d_N_hierarc_levels; i++)
			init_rect(&((*clover_rect)[i]), 2 + param->d_L_rect[i], param);
		#else
		(void) clover_rect;
		#endif
		}
	}

void free_rect_hierarc(Rectangle *most_update, Rectangle *clover_rect, GParam const *const param)
	{
	if(param->d_N_hierarc_levels == 0)
		{
		// just to avoid compiler warning of unused variables
		(void) most_update;
		(void) clover_rect;
		}
	else
		{
		int i;
		for(i = 0; i < param->d_N_hierarc_levels; i++)
			free_rect(&(most_update[i]));
		free(most_update);
		#ifdef THETA_MODE
		for(i = 0; i < param->d_N_hierarc_levels; i++)
			free_rect(&(clover_rect[i]));
		free(clover_rect);
		#else
		(void) clover_rect; // just to avoid compiler warning of unused variable
		#endif
		}
	}

void init_rect_utils(Rect_Utils *rect_aux, GParam const *const param)
	{
	// border = n -> up to n-th neighbors added to the rectangle
	// swap_rect: border = 1 (under swap, action changes only for plaquettes around the defect)
	int border = 1;
	init_rect(&(rect_aux->swap_rect), border, param);

	if(param->d_N_hierarc_levels == 0)
		{
		rect_aux->update_rect = NULL;
		#ifdef THETA_MODE
		rect_aux->clover_rect = NULL;
		#endif
		#ifdef MULTICANONICAL_MODE
		rect_aux->cooling_rect = NULL;
		rect_aux->topcharge_rect = NULL;
		#endif
		}
	else
		{
		// update_rect: border = d_L_rect
		allocate_array_Rectangle(&(rect_aux->update_rect), param->d_N_hierarc_levels, __FILE__, __LINE__);
		for(int i = 0; i < param->d_N_hierarc_levels; i++)
			{
			border = param->d_L_rect[i];
			init_rect(&(rect_aux->update_rect[i]), border, param);
			}

		#ifdef THETA_MODE
		// clover_rect: border = d_L_rect + 2
		allocate_array_Rectangle(&(rect_aux->clover_rect), param->d_N_hierarc_levels, __FILE__, __LINE__);
		for(int i = 0; i < param->d_N_hierarc_levels; i++)
			{
			border = param->d_L_rect[i] + 2;
			init_rect(&(rect_aux->clover_rect[i]), border, param);
			}
		#endif

		#ifdef MULTICANONICAL_MODE
		// topcharge_rect: border = d_L_rect + 4*coolsteps + 2
		allocate_array_Rectangle(&(rect_aux->topcharge_rect), param->d_N_hierarc_levels, __FILE__, __LINE__);
		for(int i = 0; i < param->d_N_hierarc_levels; i++)
			{
			if(param->d_topo_cooling == 0) border = param->d_L_rect[i] + 2;
			if(param->d_topo_cooling == 1) border = param->d_L_rect[i] + 2 + 4 * param->d_topo_coolsteps;
			init_rect(&(rect_aux->topcharge_rect[i]), border, param);
			}

		// cooling_rect: border = topcharge_rect + 4*(coolsteps - coolstep)
		if(param->d_topo_cooling == 1 && param->d_topo_coolsteps > 0)
			{
			allocate_array_Rectangle_pointer(&(rect_aux->cooling_rect), param->d_N_hierarc_levels, __FILE__, __LINE__);
			for(int i = 0; i < param->d_N_hierarc_levels; i++)
				{
				allocate_array_Rectangle(&(rect_aux->cooling_rect[i]), param->d_topo_coolsteps, __FILE__, __LINE__);
				for(int j = 0; j < param->d_topo_coolsteps; j++)
					{
					border = param->d_L_rect[i] + 2 + 8 * param->d_topo_coolsteps - 4 * (j + 1);
					init_rect(&(rect_aux->cooling_rect[i][j]), border, param);
					}
				}
			}
		else
			{
			rect_aux->cooling_rect = NULL;
			}
		#endif
		}
	}

void free_rect_utils(Rect_Utils *rect_aux, GParam const *const param)
	{
	free_rect(&(rect_aux->swap_rect));

	if(param->d_N_hierarc_levels > 0)
		{
		int i;

		for(i = 0; i < param->d_N_hierarc_levels; i++)
			free_rect(&(rect_aux->update_rect[i]));
		free(rect_aux->update_rect);

		#ifdef THETA_MODE
		for(i = 0; i < param->d_N_hierarc_levels; i++)
			free_rect(&(rect_aux->clover_rect[i]));
		free(rect_aux->clover_rect);
		#endif

		#ifdef MULTICANONICAL_MODE
		if(param->d_topo_cooling == 1 && param->d_topo_coolsteps > 0)
			{
			for(i = 0; i < param->d_N_hierarc_levels; i++)
				{
				for(int j = 0; j < param->d_topo_coolsteps; j++)
					free_rect(&(rect_aux->cooling_rect[i][j]));
				free(rect_aux->cooling_rect[i]);
				}
			free(rect_aux->cooling_rect);
			}

		for(i = 0; i < param->d_N_hierarc_levels; i++)
			free_rect(&(rect_aux->topcharge_rect[i]));
		free(rect_aux->topcharge_rect);
		#endif
		}
	}

// distance between sites i and j on a 1D-ring with L sites
int ring_distance(int const i, int const j, int const L)
	{
	int d = abs(i - j);
	return (d < L - d) ? d : L - d;
	}

// distance of site i from the link j -> j+1 on a 1D-ring with L sites
int link_ring_distance(int const i, int const j, int const L)
	{
	int d0 = ring_distance(i, j, L);
	int d1 = ring_distance(i, periodic_condition(j + 1, L), L);
	return (d0 < d1) ? d0 : d1;
	}

// square distance between sites i and j on a periodic lattice
double square_distance(long const i, long const j, GParam const *const param)
	{
	int x[STDIM], y[STDIM];
	double res = 0.0;

	si_to_cart(x, i, param); // i --> x
	si_to_cart(y, j, param); // j --> y

	for(int mu = 0; mu < STDIM; ++mu)
		{
		const double L = (double) param->d_size[mu];
		const double half_L = 0.5 * L;
		double d = (double) labs(x[mu] - y[mu]);

		// periodic boundary conditions
		if(d > half_L) d = L - d;

		res += d * d;
		}

	return res;
	}


#endif
