#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include "macro.h"

#include <complex.h>
#include <stdio.h>

#include "gparam.h"
#include "u1.h"
#include "u1_upd.h"
#include "su2.h"
#include "su2_upd.h"
#include "sun.h"
#include "sun_upd.h"


// select GAUGE_GROUP and function implementations according to NCOLOR:
// GAUGE_IMPL(func) becomes func_U1, func_Su2 or func_SuN depending on NCOLOR
#if NCOLOR == 1
#define GAUGE_GROUP U1
#define GAUGE_IMPL_SUFFIX _U1
#elif NCOLOR == 2
#define GAUGE_GROUP Su2
#define GAUGE_IMPL_SUFFIX _Su2
#else
#define GAUGE_GROUP SuN
#define GAUGE_IMPL_SUFFIX _SuN
#endif

#define GAUGE_PASTE_IMPL(name, suffix) name##suffix
#define GAUGE_PASTE(name, suffix) GAUGE_PASTE_IMPL(name, suffix)
#define GAUGE_IMPL(name) GAUGE_PASTE(name, GAUGE_IMPL_SUFFIX)


// basic operations

ACC_ROUTINE_SEQ
static inline void one(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(one)(A);
	}

ACC_ROUTINE_SEQ
static inline void zero(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(zero)(A);
	}

ACC_ROUTINE_SEQ
static inline void equal(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(equal)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void equal_dag(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(equal_dag)(A, B);
	}


// additions and subtractions

ACC_ROUTINE_SEQ
static inline void plus_equal(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(plus_equal)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void plus_equal_dag(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(plus_equal_dag)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void minus_equal(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(minus_equal)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void minus_equal_times_real(GAUGE_GROUP *A, GAUGE_GROUP const *const B, double const r)
	{
	GAUGE_IMPL(minus_equal_times_real)(A, B, r);
	}

ACC_ROUTINE_SEQ
static inline void minus_equal_dag(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(minus_equal_dag)(A, B);
	}


// multiplications

ACC_ROUTINE_SEQ
static inline void times_equal_real(GAUGE_GROUP *A,	double const r)
	{
	GAUGE_IMPL(times_equal_real)(A, r);
	}

ACC_ROUTINE_SEQ
static inline void times_equal_complex(GAUGE_GROUP *A, double complex const r)
	{
	GAUGE_IMPL(times_equal_complex)(A, r);
	}

ACC_ROUTINE_SEQ
static inline void times_equal(GAUGE_GROUP *A, GAUGE_GROUP const *const B)
	{
	GAUGE_IMPL(times_equal)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void times_equal_dag(GAUGE_GROUP *A, GAUGE_GROUP const *B)
	{
	GAUGE_IMPL(times_equal_dag)(A, B);
	}

ACC_ROUTINE_SEQ
static inline void times(GAUGE_GROUP *A, GAUGE_GROUP const *const B, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(times)(A, B, C);
	}

ACC_ROUTINE_SEQ
static inline void times_dag1(GAUGE_GROUP *A, GAUGE_GROUP const *const B, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(times_dag1)(A, B, C);
	}

ACC_ROUTINE_SEQ
static inline void times_dag2(GAUGE_GROUP *A, GAUGE_GROUP const *const B, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(times_dag2)(A, B, C);
	}

ACC_ROUTINE_SEQ
static inline void times_dag12(GAUGE_GROUP *A, GAUGE_GROUP const *const B, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(times_dag12)(A, B, C);
	}


// linear combinations

ACC_ROUTINE_SEQ
static inline void lin_comb(GAUGE_GROUP *A, double const b,	GAUGE_GROUP const *const B,	double const c,	GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(lin_comb)(A, b, B, c, C);
	}

ACC_ROUTINE_SEQ
static inline void lin_comb_dag1(GAUGE_GROUP *A, double const b, GAUGE_GROUP const *const B, double const c, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(lin_comb_dag1)(A, b, B, c, C);
	}

ACC_ROUTINE_SEQ
static inline void lin_comb_dag2(GAUGE_GROUP *A, double const b, GAUGE_GROUP const *const B, double const c, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(lin_comb_dag2)(A, b, B, c, C);
	}

ACC_ROUTINE_SEQ
static inline void lin_comb_dag12(GAUGE_GROUP *A, double const b, GAUGE_GROUP const *const B, double const c, GAUGE_GROUP const *const C)
	{
	GAUGE_IMPL(lin_comb_dag12)(A, b, B, c, C);
	}


// random generation

static inline void rand_matrix(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(rand_matrix)(A);
	}


// norms and traces

ACC_ROUTINE_SEQ
static inline double norm(GAUGE_GROUP const *const A)
	{
	return GAUGE_IMPL(norm)(A);
	}

ACC_ROUTINE_SEQ
static inline double relative_dist(GAUGE_GROUP const *const A, GAUGE_GROUP const *const B)
	{
	return GAUGE_IMPL(relative_dist)(A, B);
	}

ACC_ROUTINE_SEQ
static inline double retr(GAUGE_GROUP const *const A)
	{
	return GAUGE_IMPL(retr)(A);
	}

ACC_ROUTINE_SEQ
static inline double imtr(GAUGE_GROUP const *const A)
	{
	return GAUGE_IMPL(imtr)(A);
	}


// unitarization and exponentiation

ACC_ROUTINE_SEQ
static inline void unitarize(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(unitarize)(A);
	}

ACC_ROUTINE_SEQ
static inline void ta(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(ta)(A);
	}

ACC_ROUTINE_SEQ
static inline void taexp(GAUGE_GROUP *A)
	{
	GAUGE_IMPL(taexp)(A);
	}


// I/O operations

static inline void print_on_screen(GAUGE_GROUP const *const A)
	{
	GAUGE_IMPL(print_on_screen)(A);
	}

static inline void print_on_file(FILE *fp, GAUGE_GROUP const *const A)
	{
	GAUGE_IMPL(print_on_file)(fp, A);
	}

static inline void print_on_binary_file_bigen(FILE *fp, GAUGE_GROUP const *const A)
	{
	GAUGE_IMPL(print_on_binary_file_bigen)(fp, A);
	}

static inline void read_from_file(FILE *fp, GAUGE_GROUP *A)
	{
	GAUGE_IMPL(read_from_file)(fp, A);
	}

static inline void read_from_binary_file_bigen(FILE *fp, GAUGE_GROUP *A)
	{
	GAUGE_IMPL(read_from_binary_file_bigen)(fp, A);
	}


// Monte Carlo updates and cooling

static inline void single_heatbath(GAUGE_GROUP *link, GAUGE_GROUP const *const staple, GParam const *const param)
	{
	GAUGE_IMPL(single_heatbath)(link, staple, param);
	}

ACC_ROUTINE_SEQ
static inline void single_overrelaxation(GAUGE_GROUP *link, GAUGE_GROUP const *const staple)
	{
	GAUGE_IMPL(single_overrelaxation)(link, staple);
	}

ACC_ROUTINE_SEQ
static inline void cool(GAUGE_GROUP *link, GAUGE_GROUP const *const staple)
	{
	GAUGE_IMPL(cool)(link, staple);
	}


#undef GAUGE_IMPL
#undef GAUGE_PASTE
#undef GAUGE_PASTE_IMPL
#undef GAUGE_IMPL_SUFFIX

#endif
