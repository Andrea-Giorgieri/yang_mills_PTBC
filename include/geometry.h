#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "macro.h"

#include "gparam.h"

#if STDIM == 4
// 4d-only global lookup tables for topological charge and theta-term
extern const int g_signed_ord_orth_dir[4][4][2]; // (mu, nu, sood[mu][nu][0], sood[mu][nu][1]) is an even permutation of (0,1,2,3)
extern const int g_indep_perm_dir[4][3];         // (ipd[0][i], ipd[1][i], ipd[2][i], ipd[3][i]) is the i-th independent permutation of (0,1,2,3)
#endif

typedef struct Geometry
	{
	long *d_sip_to_si;                // d_sip_to_si[i]     = lexeo index of site with lexeop index i
	long **d_nnp;                     // d_nnp_loc[r][i]    = next neighbour (on the local lattice) in dir.  i of the site r
	long **d_nnm;                     // d_nnm_loc[r][i]    = next neighbour (on the local lattice) in dir. -i of the site r
	int **d_mucomp;                   // d_mucomp[mu][r]    = mu component of r
	long **d_muorth;                  // d_muorth[mu][r]    = mu-orthogonal-space component of r
	long ***d_mutimespace;            // d_mutimespace[mu][t][rsp] = r such that mu component = t, nu != mu component = rsp
	int d_orth_dir[STDIM][STDIM - 1]; // d_orth_dir[mu][i] = i-th direction orthogonal to mu
	} Geometry;


typedef struct Rectangle
	{
	long *rect_sites;
	int d_size_rect[STDIM];
	long d_vol_rect;
	} Rectangle;


// auxiliary rectangles for hierarchical updates and swaps
typedef struct Rect_Utils
	{
	Rectangle swap_rect;
	Rectangle *update_rect;    // update_rect[hierarc_levels]
	Rectangle *clover_rect;    // clover_rect[hierarc_levels]
	Rectangle *topcharge_rect; // topcharge_rect[hierarc_levels]
	Rectangle **cooling_rect;  // cooling_rect[hierarc_levels][coolsteps]
	} Rect_Utils;


// general functions
void init_geometry(Geometry *geo, GParam const *const param);

void free_geometry(Geometry *geo, GParam const *const param);

void test_geometry(Geometry const *const geo, GParam const *const param);


// functions for switching between indexing methods inside geometry.c
long cart_to_lex(int const *const cartcoord, GParam const *const param);                            // cartesian coordinates -> lexicographic index

void lex_to_cart(int *cartcoord, long lex, GParam const *const param);                              // lexicographic index -> cartesian coordinates

long cart_to_lexeo(int const *const cartcoord, GParam const *const param);                          // cartesian coordinates -> lexicographic eo index

void lexeo_to_cart(int *cartcoord, long lexeo, GParam const *const param);                          // lexicographic eo index -> cartesian coordinates

long cart_to_lexeop(int const *const cartcoord, GParam const *const param);                         // cartesian coordinates -> lexicographic eo index for pbc with odd sides

long lex_to_lexeo(long lex, GParam const *const param);                                             // lexicographic index -> lexicographic eo index

long lexeo_to_lex(long lexeo, GParam const *const param);                                           // lexicographic eo index -> lexicographic index

long cartsp_to_lexsp(int const *const ccsp, GParam const *const param);                             // spatial cartesian coordinates -> spatial lexicographic index

long cartsp_to_lexsp_mu(int const *const ccsp, int mu, GParam const *const param);                  // (nu != mu) cartesian coordinates -> (nu != mu) lexicographic index

void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const *const param);                             // spatial lexicographic index -> spatial cartesian coordinates

long cartsp_to_lexeosp(int const *const ccsp, GParam const *const param);                           // spatial cartesian coordinates -> spatial lexicographic eo index

long cartsp_to_lexeosp_mu(int const *const ccsp, int mu, GParam const *const param);                // (nu != mu) cartesian coordinates -> (nu != mu) lexicographic eo index

void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const *const param);                         // spatial lexicographic eo index -> spatial cartesian coordinates

long lexsp_to_lexeosp(long lexsp, GParam const *const param);                                       // spatial lexicographic index -> spatial lexicographic eo index

long lexeosp_to_lexsp(long lexeosp, GParam const *const param);                                     // spatial lexicographic eo index -> spatial lexicographic index

long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const *const param);                        // lexicographic eo spatial and time -> lexicographic eo index

void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const *const param);          // lex. eo index -> lex. eo spatial and t

void lexeo_to_lexeosp_and_mu(long *lexeosp, int *t, long lexeo, int mu, GParam const *const param); // lex. eo index -> lex. eo (nu != mu) and mu

long cart_to_lexeo_rect(int const *const cartcoord, Rectangle const *const rect);                   // cartesian -> lexicographic eo on a given rectangle

int dirs_to_si(int const i, int const j);                                                           //plane i-j -> single index, for twist factors


// utilities and distances with periodic boundary conditions
int is_on_defect(long const r, GParam const *const param);

int orthogonal_dir(int const mu, int const i);

int periodic_condition(int const coord, int const L_max);

int ring_distance(int const i, int const j, int const L);

int link_ring_distance(int const i, int const j, int const L);

double square_distance(long const i, long const j, GParam const *const param);


// geometry of rectangles used in the hierarchical update during parallel tempering
void init_rect(Rectangle *rect, int L_R, GParam const *const param);

void free_rect(Rectangle *rect);

void init_rect_hierarc(Rectangle **rect, Rectangle **clover_rect, GParam const *const param);

void free_rect_hierarc(Rectangle *rect, Rectangle *clover_rect, GParam const *const param);

void init_rect_utils(Rect_Utils *rect_aux, GParam const *const param);

void free_rect_utils(Rect_Utils *rect_aux, GParam const *const param);


// wrappers for switching between indexing methods outside of geometry.c

// cartesian coordinates -> single index
static inline long cart_to_si(int const *const cartcoord, GParam const *const param)
	{
	return cart_to_lexeo(cartcoord, param);
	}


// single index -> cartesian coordinates
static inline void si_to_cart(int *cartcoord, long const si, GParam const *const param)
	{
	lexeo_to_cart(cartcoord, si, param);
	}


// lexicographic -> single index
static inline long lex_to_si(long const lex, GParam const *const param)
	{
	return lex_to_lexeo(lex, param);
	}


// single index -> lexicographic
static inline long si_to_lex(long const si, GParam const *const param)
	{
	return lexeo_to_lex(si, param);
	}


// single index spatial and time -> single index tot
static inline long sisp_and_t_to_si_compute(long const sisp, int const t, GParam const *const param)
	{
	return lexeosp_and_t_to_lexeo(sisp, t, param);
	}


// single index tot -> single index spatial and time
static inline void si_to_sisp_and_t_compute(long *const sisp, int *const t, long const si, GParam const *const param)
	{
	lexeo_to_lexeosp_and_t(sisp, t, si, param);
	}


// single index tot -> single index spatial and mu
static inline void si_to_sisp_and_mu_compute(long *const sisp, int *const t, long const si, int const mu, GParam const *const param)
	{
	lexeo_to_lexeosp_and_mu(sisp, t, si, mu, param);
	}


// cartesian -> single index on rectangles
static inline long cart_to_si_rect(int const *const cartcoord, Rectangle const *const rect)
	{
	return cart_to_lexeo_rect(cartcoord, rect);
	}


// functions for accessing pre-computed indices

// next neighbour in + direction
static inline long nnp(Geometry const *const geo, long const r, int const i)
	{
	return geo->d_nnp[r][i];
	}


// next neighbour in - direction
static inline long nnm(Geometry const *const geo, long const r, int const i)
	{
	return geo->d_nnm[r][i];
	}


// single index spatial and time -> single index tot
static inline long sisp_and_t_to_si(Geometry const *const geo, long const sisp, int const t)
	{
	return geo->d_mutimespace[0][t][sisp];
	}


// single index nu != mu and mu -> single index tot
static inline long sisp_and_mu_to_si(Geometry const *const geo, long const sisp, int const t, int const mu)
	{
	return geo->d_mutimespace[mu][t][sisp];
	}


// single index tot -> single index spatial and time
static inline void si_to_sisp_and_t(long *sisp, int *t, Geometry const *const geo, long const si)
	{
	*sisp = geo->d_muorth[0][si];
	*t = geo->d_mucomp[0][si];
	}


#endif
