#ifndef GEOMETRY_H
#define GEOMETRY_H

#include"macro.h"
#include"gparam.h"

typedef struct Geometry {
	long  **d_nnp;                           // d_nnp_loc[r][i]    = next neighbour (on the local lattice) in dir.  i of the site r
	long  **d_nnm;                           // d_nnm_loc[r][i]    = next neighbour (on the local lattice) in dir. -i of the site r
	int   **d_mucomp;                        // d_mucomp[mu][r]    = mu component of r
	long  **d_muorth;                        // d_muorth[mu][r]    = mu-orthogonal-space component of r
	long ***d_mutimespace;                   // d_mutimespace[mu][t][rsp] = r such that mu component = t, nu != mu component = rsp
	int     d_orth_dir[STDIM][STDIM-1];      // d_orth_dir[mu][i] = i-th direction orthogonal to mu
	int     d_signed_ord_orth_dir[4][4][2];  // (mu, nu, sood[mu][nu][0], sood[mu][nu][1]) is an even permutation of (0,1,2,3)
	int     d_indep_perm_dir[4][3];          // (ipd[0][i], ipd[1][i], ipd[2][i], ipd[3][i]) is the i-th independent permutation of (0,1,2,3)
} Geometry;

typedef struct Rectangle {
	long *rect_sites;
	int d_size_rect[STDIM];
	long d_vol_rect;
} Rectangle;

// auxiliary rectangles for hierarchical updates and swaps
typedef struct Rect_Utils {
	Rectangle *update_rect;		// update_rect[hierarc_levels]
	Rectangle swap_rect;

	Rectangle *clover_rect;		// clover_rect[hierarc_levels]

	Rectangle **cooling_rect;	// cooling_rect[hierarc_levels][coolsteps]
	Rectangle *topcharge_rect;	// topcharge_rect[hierarc_levels]
	}	Rect_Utils;


// these are the functions to be used in shwitching between different indices
long (*cart_to_si)(int const * const cartcoord, GParam const * const param); // cartesian coordinates -> single index
void (*si_to_cart)(int *cartcoord, long si, GParam const * const param);     // single index -> cartesian coordinates
long (*lex_to_si)(long lex, GParam const * const param);          // lexicographic -> single index
long (*si_to_lex)(long si, GParam const * const param);           // lexicographic -> single index
long (*sisp_and_t_to_si_compute)(long sisp, int t, GParam const * const param);            // single index spatial and time -> single index tot
void (*si_to_sisp_and_t_compute)(long *sisp, int *t, long si, GParam const * const param); // single index tot -> single index spatial and time
void (*si_to_sisp_and_mu_compute)(long *sisp, int *t, long si, int mu, GParam const * const param); // single index tot -> single index spatial and mu
long (*cart_to_si_rect)(int const * const cartcoord, Rectangle const * const most_update); // cartesian -> single index on rectangles
int	dirs_to_si(int const i, int const j); //plane i-j -> single index, for twist factors

// general functions
int is_on_defect(long const r, GParam const * const param);
void init_indexing_lexeo(void); // has to be called before init_geometry
void init_geometry(Geometry *geo, GParam const * const param);
void free_geometry(Geometry *geo, GParam const * const param);

// next neighbour in + direction
inline long nnp(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnp[r][i];
  }

// next neighbour in - direction
inline long nnm(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnm[r][i];
  }

// single index spatial and time -> single index tot
inline long sisp_and_t_to_si(Geometry const * const geo, long sisp, int t)
  {
  return geo->d_mutimespace[0][t][sisp];
  }

// single index nu != mu and mu -> single index tot
inline long sisp_and_mu_to_si(Geometry const * const geo, long sisp, int t, int mu)
  {
  return geo->d_mutimespace[mu][t][sisp];
  }

// single index tot -> single index spatial and time
inline void si_to_sisp_and_t(long *sisp, int *t, Geometry const * const geo, long si)
  {
  *sisp=geo->d_muorth[0][si];
  *t=geo->d_mucomp[0][si];
  }

// for debug
void test_geometry(Geometry const * const geo, GParam const * const param);

//------------ these are not to be used outside geometry.c ----------------

long cart_to_lex(int const * const cartcoord, GParam const * const param);   // cartesian coordinates -> lexicographic index
void lex_to_cart(int *cartcoord, long lex, GParam const * const param);      // lexicographic index -> cartesian coordinates

long cart_to_lexeo(int const * const cartcoord, GParam const * const param); // cartesian coordinates -> lexicographic eo index
void lexeo_to_cart(int *cartcoord, long lexeo, GParam const * const param);  // lexicographic eo index -> cartesian coordinates

long lex_to_lexeo(long lex, GParam const * const param);                     //  lexicographic index -> lexicographic eo index
long lexeo_to_lex(long lexeo, GParam const * const param);                   //  lexicographic eo index -> lexicographic index

long cartsp_to_lexsp(int const * const ccsp, GParam const * const param);              // spatial cartesian coordinates -> spatial lexicographic index
long cartsp_to_lexsp_mu(int const * const ccsp, int mu, GParam const * const param);   // (nu != mu) cartesian coordinates -> (nu != mu) lexicographic index
void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const * const param);               // spatial lexicographic index -> spatial cartesian coordinates

long cartsp_to_lexeosp(int const * const ccsp, GParam const * const param);            // spatial cartesian coordinates -> spatial lexicographic eo index
long cartsp_to_lexeosp_mu(int const * const ccsp, int mu, GParam const * const param); // (nu != mu) cartesian coordinates -> (nu != mu) lexicographic eo index
void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const * const param);           // spatial lexicographic eo index -> spatial cartesian coordinates

long lexsp_to_lexeosp(long lexsp, GParam const * const param);   //  spatial lexicographic index -> spatial lexicographic eo index
long lexeosp_to_lexsp(long lexeosp, GParam const * const param); //  spatial lexicographic eo index -> spatial lexicographic index

long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const * const param);                        // lexicographic eo spatial and time -> lexicographic eo index
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const * const param);          // lex. eo index -> lex. eo spatial and t
void lexeo_to_lexeosp_and_mu(long *lexeosp, int *t, long lexeo, int mu, GParam const * const param); // lex. eo index -> lex. eo (nu != mu) and mu

long cart_to_lexeo_rect(int const * const cartcoord, Rectangle const * const most_update);

// geometry of rectangles used in the hierarchical update during parallel tempering

int periodic_condition(int const coord, int const L_max);
int orthogonal_dir(int const mu, int const i);
void init_rect(Rectangle *most_update, int const L_R, GParam const * const param);
void free_rect(Rectangle *most_update);
void init_rect_hierarc(Rectangle **most_update, Rectangle **clover_rect, GParam const * const param);
void free_rect_hierarc(Rectangle *most_update, Rectangle *clover_rect, GParam const * const param);
void init_rect_utils(Rect_Utils *rect_aux, GParam const * const param);
void free_rect_utils(Rect_Utils *rect_aux, GParam const * const param);

// for boundary conditions
int ring_distance(int const i, int const j, int const L);
int link_ring_distance(int const i, int const j, int const L);

// needed to compute chi'
double square_distance(long const i, long const j, GParam const * const param);

#endif
