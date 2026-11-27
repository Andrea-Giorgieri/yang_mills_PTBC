#ifndef SUN_H
#define SUN_H

#include "macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "endianness.h"
#include "random.h"
#include "su2.h"


// the element [i][j] is obtained as matrix.comp[m(i,j)] with m(i,j) defined in macro.h
typedef struct SuN
	{
	double complex comp[NCOLOR * NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
	} SuN;


// basic operations

// A=1
ACC_ROUTINE_SEQ
static inline void one_SuN(SuN *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] = 0.0 + 0.0 * I;
		}

	for(int i = 0; i < NCOLOR; i++)
		{
		A->comp[m(i, i)] = 1.0 + 0.0 * I;
		}
	}


// A=0
ACC_ROUTINE_SEQ
static inline void zero_SuN(SuN *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] = 0.0 + 0.0 * I;
		}
	}


// A=B
ACC_ROUTINE_SEQ
static inline void equal_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] = B->comp[i];
		}
	}


// A=B^{dag}
ACC_ROUTINE_SEQ
static inline void equal_dag_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] = conj(B->comp[m(j, i)]);
			}
		}
	}


// additions and subtractions

// A+=B
ACC_ROUTINE_SEQ
static inline void plus_equal_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] += B->comp[i];
		}
	}


// A+=B^{dag}
ACC_ROUTINE_SEQ
static inline void plus_equal_dag_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] += conj(B->comp[m(j, i)]);
			}
		}
	}


// A-=B
ACC_ROUTINE_SEQ
static inline void minus_equal_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] -= B->comp[i];
		}
	}


// A-=(r*B)
ACC_ROUTINE_SEQ
static inline void minus_equal_times_real_SuN(SuN *restrict A, SuN const *const restrict B, double const r)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] -= r * B->comp[i];
		}
	}


// A-=B^{dag}
ACC_ROUTINE_SEQ
static inline void minus_equal_dag_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] -= conj(B->comp[m(j, i)]);
			}
		}
	}


// multiplications

// A*=r
ACC_ROUTINE_SEQ
static inline void times_equal_real_SuN(SuN *restrict A, double const r)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] *= r;
		}
	}


// A*=r
ACC_ROUTINE_SEQ
static inline void times_equal_complex_SuN(SuN *restrict A, double complex const r)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] *= r;
		}
	}


// A*=B
ACC_ROUTINE_SEQ
static inline void times_equal_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			aux[j] = A->comp[m(i, j)];
			}

		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += aux[k] * B->comp[m(k, j)];
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// A*=B^{dag}
ACC_ROUTINE_SEQ
static inline void times_equal_dag_SuN(SuN *restrict A, SuN const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			aux[j] = A->comp[m(i, j)];
			}

		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += aux[k] * conj(B->comp[m(j, k)]);
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// A=B*C
ACC_ROUTINE_SEQ
static inline void times_SuN(SuN *restrict A, SuN const *const restrict B, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += B->comp[m(i, k)] * C->comp[m(k, j)];
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// A=B^{dag}*C
ACC_ROUTINE_SEQ
static inline void times_dag1_SuN(SuN *restrict A, SuN const *const restrict B, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += conj(B->comp[m(k, i)]) * C->comp[m(k, j)];
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// A=B*C^{dag}
ACC_ROUTINE_SEQ
static inline void times_dag2_SuN(SuN *restrict A, SuN const *const restrict B, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += B->comp[m(i, k)] * conj(C->comp[m(j, k)]);
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// A=B^{dag}*C^{dag}
ACC_ROUTINE_SEQ
static inline void times_dag12_SuN(SuN *restrict A, SuN const *const restrict B, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += conj(B->comp[m(k, i)]) * conj(C->comp[m(j, k)]);
				}
			A->comp[m(i, j)] = sum;
			}
		}
	}


// linear combinations

// A=b*B+c*C
ACC_ROUTINE_SEQ
static inline void lin_comb_SuN(SuN *restrict A, double const b, SuN const *const restrict B, double const c, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		A->comp[i] = b * B->comp[i] + c * C->comp[i];
		}
	}


// A=b*B^{dag}+c*C
ACC_ROUTINE_SEQ
static inline void lin_comb_dag1_SuN(SuN *restrict A, double const b, SuN const *const restrict B, double const c, SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] = b * conj(B->comp[m(j, i)]) + c * C->comp[m(i, j)];
			}
		}
	}


// A=b*B+c*C^{dag}
ACC_ROUTINE_SEQ
static inline void lin_comb_dag2_SuN(SuN *restrict A, double const b, SuN const *const restrict B, double const c,
                                     SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] = b * B->comp[m(i, j)] + c * conj(C->comp[m(j, i)]);
			}
		}
	}


// A=b*B^{dag}+c*C^{dag}
ACC_ROUTINE_SEQ
static inline void lin_comb_dag12_SuN(SuN *restrict A, double const b, SuN const *const restrict B, double const c,
                                      SuN const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT_DEV(A != B, "the same pointer is used twice");
	ASSERT_DEV(A != C, "the same pointer is used twice");
	ASSERT_DEV(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			A->comp[m(i, j)] = b * conj(B->comp[m(j, i)]) + c * conj(C->comp[m(j, i)]);
			}
		}
	}


// random generation

// SU(N) random matrix generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
ACC_ROUTINE_SEQ
static inline void rand_matrix_SuN(SuN *A, rng_t *rng)
	{
	double p0, p1, p2, p3, r2;

	one_SuN(A);

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			// SU(2) random components
			do
				{
				p0 = 1.0 - 2.0 * casuale(rng);
				p1 = 1.0 - 2.0 * casuale(rng);
				p2 = 1.0 - 2.0 * casuale(rng);
				p3 = 1.0 - 2.0 * casuale(rng);

				r2 = p0 * p0 + p1 * p1 + p2 * p2 + p3 * p3;
				} while(r2 > 1.0);

			double const invr = 1.0 / sqrt(r2);

			p0 *= invr;
			p1 *= invr;
			p2 *= invr;
			p3 *= invr;

			double complex const aux00 = p0 + p3 * I;
			double complex const aux01 = p2 + p1 * I;
			double complex const aux10 = -p2 + p1 * I;
			double complex const aux11 = p0 - p3 * I;

			for(int k = 0; k < NCOLOR; k++)
				{
				double complex *Aki = &A->comp[m(k, i)];
				double complex *Akj = &A->comp[m(k, j)];
				double complex const temp0 = *Aki * aux00 + *Akj * aux10;
				double complex const temp1 = *Aki * aux01 + *Akj * aux11;
				*Aki = temp0;
				*Akj = temp1;
				}
			}
		}
	}


// generate a matrix in the algebra of SU(N) with gaussian
// random components in the base T_i such that Tr(T_iT_j)=delta_{ij}
ACC_ROUTINE_SEQ
static inline void rand_algebra_gauss_matrix_SuN(SuN *A, rng_t *rng)
	{
	#if NCOLOR == 1
	(void) A; // just to avoid warnings
	(void) rng;
	#else
	double d1, d2, dd[NCOLOR - 1];

	zero_SuN(A);

	// out of diagonal elements
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			gauss2(&d1, &d2, rng);
			A->comp[m(i, j)] = d1 - d2 * I;
			A->comp[m(j, i)] = d1 + d2 * I;
			}
		}

	// random numbers to be used in the diagonal
	for(int i = 0; i < NCOLOR - 1; i++)
		{
		dd[i] = gauss1(rng);
		}

	// diagonal
	#if NCOLOR == 2
	A->comp[m(0, 0)] = dd[0];
	A->comp[m(1, 1)] = -dd[0];
	#else
	double const factor = sqrt(2.0 / (double) (NCOLOR * NCOLOR - NCOLOR));
	for(int i = 0; i < NCOLOR - 2; i++)
		{
		A->comp[m(i, i)] += dd[i];
		A->comp[m(i + 1, i + 1)] -= dd[i];
		}
	for(int i = 0; i < NCOLOR - 1; i++)
		{
		A->comp[m(i, i)] += factor * dd[NCOLOR - 2];
		}
	A->comp[m(NCOLOR - 1, NCOLOR - 1)] = factor * (1.0 - (double) NCOLOR) * dd[NCOLOR - 2];
	#endif

	times_equal_real_SuN(A, 0.7071067811865475244008); // *= 1 / sqrt(2)
	#endif
	}


// norms and traces

// l2 norm
ACC_ROUTINE_SEQ
static inline double norm_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	double res = 0.0;
	for(int i = 0; i < NCOLOR * NCOLOR; i++)
		{
		double const aux = cabs(A->comp[i]);
		res += aux * aux;
		}
	return sqrt(res);
	}


// ReTr[A]/N
ACC_ROUTINE_SEQ
static inline double retr_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	double complex tr = 0.0 + 0.0 * I;
	for(int i = 0; i < NCOLOR; i++)
		{
		tr += A->comp[m(i, i)];
		}
	return creal(tr) / (double) NCOLOR;
	}


// ImTr[A]/N
ACC_ROUTINE_SEQ
static inline double imtr_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	double complex tr = 0.0 + 0.0 * I;
	for(int i = 0; i < NCOLOR; i++)
		{
		tr += A->comp[m(i, i)];
		}
	return cimag(tr) / (double) NCOLOR;
	}


// ArgTr[A]
ACC_ROUTINE_SEQ
static inline double argtr_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	double complex tr = 0.0 + 0.0 * I;
	for(int i = 0; i < NCOLOR; i++)
		{
		tr += A->comp[m(i, i)];
		}
	return atan2(cimag(tr), creal(tr)); // = carg(tr), but carg is not available on GPU
	}


// Tr[A * B^{dag}] / N
ACC_ROUTINE_SEQ
static inline double complex tr_times_dag_SuN(SuN const *const restrict A, SuN const *const restrict B)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	double complex tr = 0.0 + 0.0 * I;
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			tr += A->comp[m(i, j)] * conj(B->comp[m(i, j)]);
			}
		}
	return tr / (double) NCOLOR;
	}


// norm(A - B) / (1/2 * \sqrt{norm(A - 1)**2 + norm(B - 1)**2})
ACC_ROUTINE_SEQ
static inline double relative_dist_SuN(SuN const *const restrict A, SuN const *const restrict B)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	#ifndef DEBUG // use more efficient expression assuming A and B are exactly in SU(N)

	double complex tr_AxB = 0.0 + 0.0 * I;
	double complex tr_ApB = 0.0 + 0.0 * I;
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			tr_AxB += (A->comp[m(i, j)]) * conj(B->comp[m(i, j)]);
			}
		tr_ApB += (A->comp[m(i, i)]) + (B->comp[m(i, i)]);
		}
	double const aux_ApB = 2.0 - creal(tr_ApB) / (double) NCOLOR;
	double const aux_AxB = 1.0 - creal(tr_AxB) / (double) NCOLOR;

	#else

	double aux_AxB = 0.0;
	double aux_ApB = 0.0;
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double aux2;
			double aux1 = cabs(A->comp[m(i, j)] - B->comp[m(i, j)]);
			aux_AxB += aux1 * aux1;
			if(i == j)
				{
				aux1 = cabs(A->comp[m(i, j)] - 1.0);
				aux2 = cabs(B->comp[m(i, j)] - 1.0);
				}
			else
				{
				aux1 = cabs(A->comp[m(i, j)]);
				aux2 = cabs(B->comp[m(i, j)]);
				}
			aux_ApB += aux1 * aux1 + aux2 * aux2;
			}
		}

	#endif

	double const check = 2 * MIN_VALUE * (double) NCOLOR;
	if(aux_ApB <= check || aux_AxB <= check)
		{
		return 0.0;
		}
	return sqrt(2.0 * aux_AxB / aux_ApB);
	}


// LU decomposition with partial pivoting
// from Numerical Recipes in C, pag 46
ACC_ROUTINE_SEQ
static inline void LU_SuN(SuN const *const restrict A, SuN *restrict res, int *restrict sign)
	{
	int i, j, k;
	double big, temp;
	double complex sum, dum;
	double vv[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

	int imax = 0;
	equal_SuN(res, A);

	*sign = 1;
	for(i = 0; i < NCOLOR; i++)
		{
		big = 0.0;
		for(j = 0; j < NCOLOR; j++)
			{
			temp = cabs(res->comp[m(i, j)]);
			if(temp > big) big = temp;
			}
		vv[i] = 1.0 / big;
		}

	for(j = 0; j < NCOLOR; j++)
		{
		for(i = 0; i < j; i++)
			{
			sum = res->comp[m(i, j)];
			for(k = 0; k < i; k++)
				{
				sum -= res->comp[m(i, k)] * res->comp[m(k, j)];
				}
			res->comp[m(i, j)] = sum;
			}

		big = 0.0;
		for(i = j; i < NCOLOR; i++)
			{
			sum = res->comp[m(i, j)];
			for(k = 0; k < j; k++)
				{
				sum -= res->comp[m(i, k)] * res->comp[m(k, j)];
				}
			res->comp[m(i, j)] = sum;

			temp = vv[i] * cabs(sum);
			if(temp >= big)
				{
				big = temp;
				imax = i;
				}
			}

		if(j != imax)
			{
			for(k = 0; k < NCOLOR; k++)
				{
				dum = res->comp[m(imax, k)];
				res->comp[m(imax, k)] = res->comp[m(j, k)];
				res->comp[m(j, k)] = dum;
				}
			*sign *= -1;
			vv[imax] = vv[j];
			}

		if(j != NCOLOR - 1)
			{
			dum = (1.0 + 0.0 * I) / res->comp[m(j, j)];
			for(i = j + 1; i < NCOLOR; i++)
				{
				res->comp[m(i, j)] *= dum;
				}
			}
		}
	}


// Det[A]
ACC_ROUTINE_SEQ
static inline complex double det_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	#if NCOLOR == 3
	complex double res = 0.0 + 0.0 * I;

	res += A->comp[m(0, 0)] * A->comp[m(1, 1)] * A->comp[m(2, 2)];
	res += A->comp[m(1, 0)] * A->comp[m(2, 1)] * A->comp[m(0, 2)];
	res += A->comp[m(2, 0)] * A->comp[m(0, 1)] * A->comp[m(1, 2)];
	res -= A->comp[m(2, 0)] * A->comp[m(1, 1)] * A->comp[m(0, 2)];
	res -= A->comp[m(1, 0)] * A->comp[m(0, 1)] * A->comp[m(2, 2)];
	res -= A->comp[m(0, 0)] * A->comp[m(2, 1)] * A->comp[m(1, 2)];

	#else

	int i;
	double complex res;
	SuN lu;

	LU_SuN(A, &lu, &i);

	if(i > 0)
		{
		res = 1.0 + 0.0 * I;
		}
	else
		{
		res = -1.0 + 0.0 * I;
		}

	for(i = 0; i < NCOLOR; i++)
		{
		res *= lu.comp[m(i, i)];
		}

	#endif

	return res;
	}


// gives 0 if the matrix is in SU(N) and 1 otherwise
ACC_ROUTINE_SEQ
static inline int scheck_SuN(SuN const *const restrict A)
	{
	int res = 0;

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex aux = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				aux += A->comp[m(i, k)] * conj(A->comp[m(j, k)]);
				}
			if(i == j) aux -= 1.0 + 0.0 * I;
			if(cabs(aux) > MIN_VALUE) res = 1;
			}
		}

	if(res == 0)
		{
		if(cabs(det_SuN(A) - 1) > MIN_VALUE)
			{
			res = 1;
			}
		}

	return res;
	}


// fused products for the staples: no temporary matrix. Every element is computed with the same
// operations, in the same order, as the sequence of the separate functions named in the comment.

// A+=r*(B*C^{dag}) = times_dag2_SuN(aux, B, C); times_equal_complex_SuN(aux, r); plus_equal_SuN(A, aux)
ACC_ROUTINE_SEQ
static inline void plus_equal_times_dag2_complex_SuN(SuN *restrict A, SuN const *const restrict B,
                                                     SuN const *const restrict C, double complex const r)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += B->comp[m(i, k)] * conj(C->comp[m(j, k)]);
				}
			sum *= r;
			A->comp[m(i, j)] += sum;
			}
		}
	}


// A+=r*(B*C) = times_SuN(aux, B, C); times_equal_complex_SuN(aux, r); plus_equal_SuN(A, aux)
ACC_ROUTINE_SEQ
static inline void plus_equal_times_complex_SuN(SuN *restrict A, SuN const *const restrict B,
                                                SuN const *const restrict C, double complex const r)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex sum = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				sum += B->comp[m(i, k)] * C->comp[m(k, j)];
				}
			sum *= r;
			A->comp[m(i, j)] += sum;
			}
		}
	}


// SU(2) subgroups

// given the NxN matrix "in", extracts the i, j lines and column and
// gives the real number "xi" and the SU(2) matrix "u" such that
// 4 xi^2 = redet2[s-s^(dag)+1*tr(s^(dag))]
// u = [s-s^(dag)+1*tr(s^(dag))]/2/xi
// (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))
// (the same from the four elements in_ii, in_ij, in_ji, in_jj of "in", used by ennetodue and ennetodue_times)
ACC_ROUTINE_SEQ
static inline void ennetodue_elements(double complex const in_ii, double complex const in_ij, double complex const in_ji,
                                      double complex const in_jj, double *xi, Su2 *u)
	{
	double s[2][2][2], aux_re[2][2], aux_im[2][2];

	s[0][0][0] = creal(in_ii);
	s[0][0][1] = cimag(in_ii);

	s[0][1][0] = creal(in_ij);
	s[0][1][1] = cimag(in_ij);

	s[1][0][0] = creal(in_ji);
	s[1][0][1] = cimag(in_ji);

	s[1][1][0] = creal(in_jj);
	s[1][1][1] = cimag(in_jj);

	aux_re[0][0] = s[0][0][0] + s[1][1][0];
	aux_im[0][0] = s[0][0][1] - s[1][1][1];

	aux_re[0][1] = s[0][1][0] - s[1][0][0];
	aux_im[0][1] = s[0][1][1] + s[1][0][1];

	aux_re[1][0] = s[1][0][0] - s[0][1][0];
	aux_im[1][0] = s[1][0][1] + s[0][1][1];

	aux_re[1][1] = s[0][0][0] + s[1][1][0];
	aux_im[1][1] = s[1][1][1] - s[0][0][1];

	double const p = sqrt(aux_re[0][0] * aux_re[1][1] - aux_im[0][0] * aux_im[1][1] - aux_re[0][1] * aux_re[1][0] + aux_im[0][1] * aux_im[1][0]);

	*xi = p / 2.0;

	if(*xi > MIN_VALUE)
		{
		double const inv_p = 1 / p;
		aux_re[0][0] *= inv_p;
		aux_im[0][1] *= inv_p;
		aux_re[0][1] *= inv_p;
		aux_im[0][0] *= inv_p;
		}

	u->comp[0] = aux_re[0][0];
	u->comp[1] = aux_im[0][1];
	u->comp[2] = aux_re[0][1];
	u->comp[3] = aux_im[0][0];
	}


ACC_ROUTINE_SEQ
static inline void ennetodue(SuN const *const in, int const i, int const j, double *xi, Su2 *u)
	{
	ennetodue_elements(in->comp[m(i, i)], in->comp[m(i, j)], in->comp[m(j, i)], in->comp[m(j, j)], xi, u);
	}


// element (a, b) of B*C
ACC_ROUTINE_SEQ
static inline double complex times_element_SuN(SuN const *const B, SuN const *const C, int const a, int const b)
	{
	double complex sum = 0.0 + 0.0 * I;
	for(int k = 0; k < NCOLOR; k++)
		{
		sum += B->comp[m(a, k)] * C->comp[m(k, b)];
		}
	return sum;
	}


// ennetodue of the product B*C without forming it: only the four elements needed by the (i, j) subgroup are computed
// (4 N complex products instead of N^3; same values as times_equal_SuN followed by ennetodue)
ACC_ROUTINE_SEQ
static inline void ennetodue_times(SuN const *const B, SuN const *const C, int const i, int const j, double *xi, Su2 *u)
	{
	ennetodue_elements(times_element_SuN(B, C, i, i), times_element_SuN(B, C, i, j), times_element_SuN(B, C, j, i),
	                   times_element_SuN(B, C, j, j), xi, u);
	}


// given a 2x2 matrix, extend it to an NxN matrix with 1 on the diagonal
ACC_ROUTINE_SEQ
static inline void duetoenne(Su2 const *const in, int const i, int const j, SuN *out)
	{
	one_SuN(out);

	out->comp[m(i, i)] = in->comp[0] + in->comp[3] * I;
	out->comp[m(i, j)] = in->comp[2] + in->comp[1] * I;
	out->comp[m(j, i)] = -in->comp[2] + in->comp[1] * I;
	out->comp[m(j, j)] = in->comp[0] - in->comp[3] * I;
	}


// cooling, unitarization and exponentiation

// cooling
ACC_ROUTINE_SEQ
static inline void cool_SuN(SuN *link, SuN const *const staple)
	{
	SuN prod;
	Su2 u, udag;
	double complex temp0, temp1;
	double aux;

	equal_SuN(&prod, staple);     // prod=staple
	times_equal_SuN(&prod, link); // prod=staple*link

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			ennetodue(&prod, i, j, &aux, &u); // aux=xi unused
			equal_dag_Su2(&udag, &u);

			double complex const fii = udag.comp[0] + udag.comp[3] * I;
			double complex const fij = udag.comp[2] + udag.comp[1] * I;
			double complex const fji = -udag.comp[2] + udag.comp[1] * I;
			double complex const fjj = udag.comp[0] - udag.comp[3] * I;

			// link*=final
			for(int k = 0; k < NCOLOR; k++)
				{
				temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
				temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
				link->comp[m(k, i)] = temp0;
				link->comp[m(k, j)] = temp1;
				}

			// prod*=final
			for(int k = 0; k < NCOLOR; k++)
				{
				temp0 = prod.comp[m(k, i)] * fii + prod.comp[m(k, j)] * fji;
				temp1 = prod.comp[m(k, i)] * fij + prod.comp[m(k, j)] * fjj;
				prod.comp[m(k, i)] = temp0;
				prod.comp[m(k, j)] = temp1;
				}
			}
		}

	// Maximize ReTr[staple * C *link] for C \in Z(SU(N)) and update link *= C
	// \phi \equiv carg(Tr[staple * link]) => C = \exp{-i * 2\pi/N * round(\phi / (2*\pi/N))}
	#if NCOLOR > 3
	aux = argtr_SuN(&prod);                        // aux = phi
	aux = round(aux / PI2_N) * PI2_N;              // round aux to nearest center phase (PI2_N = 2*pi/N in marco.h)
	times_equal_complex_SuN(link, cexp(-I * aux)); // link *= exp(-i * aux)
	#endif
	}


// TODO: bugged cooling for testing, remove
ACC_ROUTINE_SEQ
static inline void bad_cool_SuN(SuN *link, SuN const *const staple)
	{
	SuN prod;
	Su2 u, udag;
	double complex temp0, temp1;
	double aux;

	equal_SuN(&prod, staple);     // prod=staple
	times_equal_SuN(&prod, link); // prod=staple*link

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			ennetodue(&prod, i, j, &aux, &u); // aux=xi unused
			equal_dag_Su2(&udag, &u);

			double complex const fii = udag.comp[0] + udag.comp[3] * I;
			double complex const fij = udag.comp[2] + udag.comp[1] * I;
			double complex const fji = -udag.comp[2] + udag.comp[1] * I;
			double complex const fjj = udag.comp[0] - udag.comp[3] * I;

			// link*=final
			for(int k = 0; k < NCOLOR; k++)
				{
				temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
				temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
				link->comp[m(k, i)] = temp0;
				link->comp[m(k, j)] = temp1;
				}

			// prod*=final
			for(int k = 0; k < NCOLOR; k++)
				{
				temp0 = prod.comp[m(k, i)] * fii + prod.comp[m(k, j)] * fji;
				temp1 = prod.comp[m(k, i)] * fij + prod.comp[m(k, j)] * fjj;
				prod.comp[m(k, i)] = temp0;
				prod.comp[m(k, j)] = temp1;
				}
			}
		}
	}


// unitarize A
ACC_ROUTINE_SEQ
static inline void unitarize_SuN(SuN *restrict A)
	{
	SuN F;                   // F = A^{dag}, force to unitarize A by cooling
	SuN G, G_old;            // current and previous guess for unitarized A
	SuN H, H_copy, H_square; // helpers to check convergence of unitarization

	// check if A needs re-unitarization: check_SuN(A) passes (=0) if
	// |A * A^{dag} - 1| < MIN_VALUE and |det(A) - 1| < MIN_VALUE
	if(scheck_SuN(A) == 1)
		{
		// use A^{dag} as force
		equal_dag_SuN(&F, A);

		// guess initialized to identity
		one_SuN(&G);
		double check = 1.0;
		while(check > MIN_VALUE)
			{
			// store old guess
			equal_SuN(&G_old, &G);

			// get new guess by cooling
			cool_SuN(&G, &F);

			// calculate the distance between old guess G_old and new guess G:
			// check = sqrt(|ReTr[(G-G_old)^2]|/N^2)
			equal_SuN(&H, &G);
			minus_equal_SuN(&H, &G_old);
			equal_SuN(&H_copy, &H);
			times_SuN(&H_square, &H, &H_copy);
			check = sqrt(fabs(retr_SuN(&H_square)) / (double) NCOLOR);
			}

		// replace A with G (U if C was applied)
		equal_SuN(A, &G);
		}
	}


// takes the traceless antihermitian part
ACC_ROUTINE_SEQ
static inline void ta_SuN(SuN *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	SuN aux, aux1;

	equal_SuN(&aux, A);
	equal_dag_SuN(&aux1, A);
	minus_equal_SuN(&aux, &aux1);
	times_equal_real_SuN(&aux, 0.5); // now aux=(A-A^{dag})/2

	double complex trace = aux.comp[m(0, 0)];
	for(int i = 1; i < NCOLOR; i++)
		{
		trace += aux.comp[m(i, i)];
		}
	trace /= (double) NCOLOR;

	for(int i = 0; i < NCOLOR; i++)
		{
		aux.comp[m(i, i)] -= trace;
		}

	equal_SuN(A, &aux);
	}


// exponential of the traceless antihermitian part
ACC_ROUTINE_SEQ
static inline void taexp_SuN(SuN *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	SuN aux, uno, res;

	equal_SuN(&aux, A);
	ta_SuN(&aux);

	one_SuN(&uno);

	// now aux is the traceless antihermitian part of the initial matrix
	// and we use
	// exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

	equal_SuN(&res, &aux);
	times_equal_real_SuN(&res, 1.0 / 5.0);
	plus_equal_SuN(&res, &uno);

	times_equal_SuN(&res, &aux);
	times_equal_real_SuN(&res, 1.0 / 4.0);
	plus_equal_SuN(&res, &uno);

	times_equal_SuN(&res, &aux);
	times_equal_real_SuN(&res, 1.0 / 3.0);
	plus_equal_SuN(&res, &uno);

	times_equal_SuN(&res, &aux);
	times_equal_real_SuN(&res, 1.0 / 2.0);
	plus_equal_SuN(&res, &uno);

	times_equal_SuN(&res, &aux);
	plus_equal_SuN(&res, &uno);

	unitarize_SuN(&res);
	equal_SuN(A, &res);
	}


// return 0 if matrix is traceless antihermitian, 1 otherwise
ACC_ROUTINE_SEQ
static inline int ta_check_SuN(SuN const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	int res = 0;

	double complex aux = 0.0;
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			aux += A->comp[m(i, j)] + conj(A->comp[m(j, i)]);
			}
		}
	if(cabs(aux) > MIN_VALUE) res = 1;

	aux = 0.0;
	for(int i = 0; i < NCOLOR; i++)
		{
		aux += A->comp[m(i, i)];
		}
	if(cabs(aux) > MIN_VALUE) res = 1;

	return res;
	}


// exponential of a TA matrix
ACC_ROUTINE_SEQ
static inline void exp_of_ta_SuN(SuN *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	// we use
	// exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

	#ifdef DEBUG
	ASSERT_DEV(ta_check_SuN(A) == 0, "exponential of a non-TA matrix");
	#endif

	SuN aux, uno;

	one_SuN(&uno);
	equal_SuN(&aux, A); // in aux the initial matrix is stored

	times_equal_real_SuN(A, 1.0 / 5.0);
	plus_equal_SuN(A, &uno);

	times_equal_SuN(A, &aux);
	times_equal_real_SuN(A, 1.0 / 4.0);
	plus_equal_SuN(A, &uno);

	times_equal_SuN(A, &aux);
	times_equal_real_SuN(A, 1.0 / 3.0);
	plus_equal_SuN(A, &uno);

	times_equal_SuN(A, &aux);
	times_equal_real_SuN(A, 1.0 / 2.0);
	plus_equal_SuN(A, &uno);

	times_equal_SuN(A, &aux);
	plus_equal_SuN(A, &uno);

	unitarize_SuN(A);
	}


// I/O operations

// print on screen
static inline void print_on_screen_SuN(SuN const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			fprintf(stdout, "(% 5.3f % 5.3f) ", creal(A->comp[m(i, j)]), cimag(A->comp[m(i, j)]));
			}
		fprintf(stdout, "\n");
		}
	fprintf(stdout, "\n");
	}


// print on file
static inline void print_on_file_SuN(FILE *fp, SuN const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			int const err = fprintf(fp, "% 18.12e % 18.12e ", creal(A->comp[m(i, j)]), cimag(A->comp[m(i, j)]));
			REQUIRE(err >= 0, "failed to write an SU(N) matrix on a file");
			}
		}
	fprintf(fp, "\n");
	}


// print on binary file without changing endiannes
static inline void print_on_binary_file_noswap_SuN(FILE *fp, SuN const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double re = creal(A->comp[m(i, j)]);
			double im = cimag(A->comp[m(i, j)]);

			size_t err = fwrite(&re, sizeof(double), 1, fp);
			REQUIRE(err == 1, "failed to write an SU(N) matrix on a file in binary mode");
			err = fwrite(&im, sizeof(double), 1, fp);
			REQUIRE(err == 1, "failed to write an SU(N) matrix on a file in binary mode");
			}
		}
	}


// print on binary file changing endiannes
static inline void print_on_binary_file_swap_SuN(FILE *fp, SuN const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double re = creal(A->comp[m(i, j)]);
			double im = cimag(A->comp[m(i, j)]);

			SwapBytesDouble(&re);
			SwapBytesDouble(&im);

			size_t err = fwrite(&re, sizeof(double), 1, fp);
			REQUIRE(err == 1, "failed to write an SU(N) matrix on a file in binary mode");
			err = fwrite(&im, sizeof(double), 1, fp);
			REQUIRE(err == 1, "failed to write an SU(N) matrix on a file in binary mode");
			}
		}
	}


// print on binary file in bigendian
static inline void print_on_binary_file_bigen_SuN(FILE *fp, SuN const *const A)
	{
	if(endian() == 0) // little endian machine
		{
		print_on_binary_file_swap_SuN(fp, A);
		}
	else
		{
		print_on_binary_file_noswap_SuN(fp, A);
		}
	}


// read from file
static inline void read_from_file_SuN(FILE *fp, SuN *A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double re, im;
			int const err = fscanf(fp, "%lg %lg", &re, &im);
			REQUIRE(err == 2, "failed to read the (%d, %d) component of an SU(N) matrix from file", i, j);
			A->comp[m(i, j)] = re + im * I;
			}
		}
	}


// read from binary file without changing endiannes
static inline void read_from_binary_file_noswap_SuN(FILE *fp, SuN *A)
	{
	double re, im;
	double aux[2];

	size_t err = 0;
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			err += fread(&re, sizeof(double), 1, fp);
			err += fread(&im, sizeof(double), 1, fp);
			aux[0] = re;
			aux[1] = im;

			memcpy((void *) &A->comp[m(i, j)], (void *) aux, sizeof(aux));
			//equivalent to A->comp[m(i,j)]=re+im*I;
			}
		}
	REQUIRE(err == 2 * NCOLOR * NCOLOR, "failed to read an SU(N) matrix from a file in binary mode");
	}


// read from binary file changing endianness
static inline void read_from_binary_file_swap_SuN(FILE *fp, SuN *A)
	{
	double re, im;
	double aux[2];

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			size_t err = 0;
			err += fread(&re, sizeof(double), 1, fp);
			err += fread(&im, sizeof(double), 1, fp);
			REQUIRE(err == 2, "failed to read the (%d, %d) component of an SU(N) matrix from a file", i, j);

			SwapBytesDouble(&re);
			SwapBytesDouble(&im);
			aux[0] = re;
			aux[1] = im;

			memcpy((void *) &A->comp[m(i, j)], (void *) aux, sizeof(aux));
			// equivalent to A->comp[m(i,j)]=re+im*I;
			}
		}
	}


// read from binary file written in bigendian
static inline void read_from_binary_file_bigen_SuN(FILE *fp, SuN *A)
	{
	if(endian() == 0) // little endian machine
		{
		read_from_binary_file_swap_SuN(fp, A);
		}
	else
		{
		read_from_binary_file_noswap_SuN(fp, A);
		}
	}


#endif
