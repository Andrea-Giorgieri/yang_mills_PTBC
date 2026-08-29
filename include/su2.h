#ifndef SU2_H
#define SU2_H

#include "macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "endianness.h"
#include "random.h"


// An SU(2) matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
typedef struct Su2
	{
	double comp[4] __attribute__((aligned(DOUBLE_ALIGN)));
	} Su2;

// basic operations

// A from vec
static inline void init_Su2(Su2 *restrict A, double vec[4])
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = vec[0];
	A->comp[1] = vec[1];
	A->comp[2] = vec[2];
	A->comp[3] = vec[3];
	}


// A=1
static inline void one_Su2(Su2 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = 1.0;
	A->comp[1] = 0.0;
	A->comp[2] = 0.0;
	A->comp[3] = 0.0;
	}


// A=0
static inline void zero_Su2(Su2 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = 0.0;
	A->comp[1] = 0.0;
	A->comp[2] = 0.0;
	A->comp[3] = 0.0;
	}


// A=B
static inline void equal_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0];
	A->comp[1] = B->comp[1];
	A->comp[2] = B->comp[2];
	A->comp[3] = B->comp[3];
	}


// A=B^{dag}
static inline void equal_dag_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0];
	A->comp[1] = -B->comp[1];
	A->comp[2] = -B->comp[2];
	A->comp[3] = -B->comp[3];
	}


// additions and subtractions

// A+=B
static inline void plus_equal_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] += B->comp[0];
	A->comp[1] += B->comp[1];
	A->comp[2] += B->comp[2];
	A->comp[3] += B->comp[3];
	}


// A+=B^{dag}
static inline void plus_equal_dag_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] += B->comp[0];
	A->comp[1] -= B->comp[1];
	A->comp[2] -= B->comp[2];
	A->comp[3] -= B->comp[3];
	}


// A-=B
static inline void minus_equal_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] -= B->comp[0];
	A->comp[1] -= B->comp[1];
	A->comp[2] -= B->comp[2];
	A->comp[3] -= B->comp[3];
	}


// A-=(r*B)
static inline void minus_equal_times_real_Su2(Su2 *restrict A, Su2 const *const restrict B, double const r)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] -= r * B->comp[0];
	A->comp[1] -= r * B->comp[1];
	A->comp[2] -= r * B->comp[2];
	A->comp[3] -= r * B->comp[3];
	}


// A-=B^{dag}
static inline void minus_equal_dag_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] -= B->comp[0];
	A->comp[1] += B->comp[1];
	A->comp[2] += B->comp[2];
	A->comp[3] += B->comp[3];
	}


// multiplications

// A*=r
static inline void times_equal_real_Su2(Su2 *restrict A, double const r)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] *= r;
	A->comp[1] *= r;
	A->comp[2] *= r;
	A->comp[3] *= r;
	}


// A*=r
static inline void times_equal_complex_Su2(Su2 *restrict A, double complex const r)
	{
	#ifdef DEBUG
	ASSERT(fabs(cimag(r)) < MIN_VALUE, "trying to multiply an SU(2) matrix by a non-real number");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] *= creal(r);
	A->comp[1] *= creal(r);
	A->comp[2] *= creal(r);
	A->comp[3] *= creal(r);
	}


// A*=B
static inline void times_equal_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	register double const a0 = A->comp[0];
	register double const a1 = A->comp[1];
	register double const a2 = A->comp[2];
	register double const a3 = A->comp[3];

	A->comp[0] = a0 * B->comp[0] - a1 * B->comp[1] - a2 * B->comp[2] - a3 * B->comp[3];
	A->comp[1] = a0 * B->comp[1] + a1 * B->comp[0] - a2 * B->comp[3] + a3 * B->comp[2];
	A->comp[2] = a0 * B->comp[2] + a2 * B->comp[0] + a1 * B->comp[3] - a3 * B->comp[1];
	A->comp[3] = a0 * B->comp[3] + a3 * B->comp[0] - a1 * B->comp[2] + a2 * B->comp[1];
	}


// A*=B^{dag}
static inline void times_equal_dag_Su2(Su2 *restrict A, Su2 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	register double const a0 = A->comp[0];
	register double const a1 = A->comp[1];
	register double const a2 = A->comp[2];
	register double const a3 = A->comp[3];

	A->comp[0] = a0 * B->comp[0] + a1 * B->comp[1] + a2 * B->comp[2] + a3 * B->comp[3];
	A->comp[1] = -a0 * B->comp[1] + a1 * B->comp[0] + a2 * B->comp[3] - a3 * B->comp[2];
	A->comp[2] = -a0 * B->comp[2] + a2 * B->comp[0] - a1 * B->comp[3] + a3 * B->comp[1];
	A->comp[3] = -a0 * B->comp[3] + a3 * B->comp[0] + a1 * B->comp[2] - a2 * B->comp[1];
	}


// A=B*C
static inline void times_Su2(Su2 *restrict A, Su2 const *const restrict B, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0] * C->comp[0] - B->comp[1] * C->comp[1] - B->comp[2] * C->comp[2] - B->comp[3] * C->comp[3];
	A->comp[1] = B->comp[0] * C->comp[1] + B->comp[1] * C->comp[0] - B->comp[2] * C->comp[3] + B->comp[3] * C->comp[2];
	A->comp[2] = B->comp[0] * C->comp[2] + B->comp[2] * C->comp[0] + B->comp[1] * C->comp[3] - B->comp[3] * C->comp[1];
	A->comp[3] = B->comp[0] * C->comp[3] + B->comp[3] * C->comp[0] - B->comp[1] * C->comp[2] + B->comp[2] * C->comp[1];
	}


// A=B^{dag}*C
static inline void times_dag1_Su2(Su2 *restrict A, Su2 const *const restrict B, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0] * C->comp[0] + B->comp[1] * C->comp[1] + B->comp[2] * C->comp[2] + B->comp[3] * C->comp[3];
	A->comp[1] = B->comp[0] * C->comp[1] - B->comp[1] * C->comp[0] + B->comp[2] * C->comp[3] - B->comp[3] * C->comp[2];
	A->comp[2] = B->comp[0] * C->comp[2] - B->comp[2] * C->comp[0] - B->comp[1] * C->comp[3] + B->comp[3] * C->comp[1];
	A->comp[3] = B->comp[0] * C->comp[3] - B->comp[3] * C->comp[0] + B->comp[1] * C->comp[2] - B->comp[2] * C->comp[1];
	}


// A=B*C^{dag}
static inline void times_dag2_Su2(Su2 *restrict A, Su2 const *const restrict B, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0] * C->comp[0] + B->comp[1] * C->comp[1] + B->comp[2] * C->comp[2] + B->comp[3] * C->comp[3];
	A->comp[1] = -B->comp[0] * C->comp[1] + B->comp[1] * C->comp[0] + B->comp[2] * C->comp[3] - B->comp[3] * C->comp[2];
	A->comp[2] = -B->comp[0] * C->comp[2] + B->comp[2] * C->comp[0] - B->comp[1] * C->comp[3] + B->comp[3] * C->comp[1];
	A->comp[3] = -B->comp[0] * C->comp[3] + B->comp[3] * C->comp[0] + B->comp[1] * C->comp[2] - B->comp[2] * C->comp[1];
	}


// A=B^{dag}*C^{dag}
static inline void times_dag12_Su2(Su2 *restrict A, Su2 const *const restrict B, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = B->comp[0] * C->comp[0] - B->comp[1] * C->comp[1] - B->comp[2] * C->comp[2] - B->comp[3] * C->comp[3];
	A->comp[1] = -B->comp[0] * C->comp[1] - B->comp[1] * C->comp[0] - B->comp[2] * C->comp[3] + B->comp[3] * C->comp[2];
	A->comp[2] = -B->comp[0] * C->comp[2] - B->comp[2] * C->comp[0] + B->comp[1] * C->comp[3] - B->comp[3] * C->comp[1];
	A->comp[3] = -B->comp[0] * C->comp[3] - B->comp[3] * C->comp[0] - B->comp[1] * C->comp[2] + B->comp[2] * C->comp[1];
	}


// linear combinations

// A=b*B+c*C
static inline void lin_comb_Su2(Su2 *restrict A, double const b, Su2 const *const restrict B, double const c, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = b * B->comp[0] + c * C->comp[0];
	A->comp[1] = b * B->comp[1] + c * C->comp[1];
	A->comp[2] = b * B->comp[2] + c * C->comp[2];
	A->comp[3] = b * B->comp[3] + c * C->comp[3];
	}


// A=b*B^{dag}+c*C
static inline void lin_comb_dag1_Su2(Su2 *restrict A, double const b, Su2 const *const restrict B, double const c, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = b * B->comp[0] + c * C->comp[0];
	A->comp[1] = -b * B->comp[1] + c * C->comp[1];
	A->comp[2] = -b * B->comp[2] + c * C->comp[2];
	A->comp[3] = -b * B->comp[3] + c * C->comp[3];
	}


// A=b*B+c*C^{dag}
static inline void lin_comb_dag2_Su2(Su2 *restrict A, double const b, Su2 const *const restrict B, double const c, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = b * B->comp[0] + c * C->comp[0];
	A->comp[1] = b * B->comp[1] - c * C->comp[1];
	A->comp[2] = b * B->comp[2] - c * C->comp[2];
	A->comp[3] = b * B->comp[3] - c * C->comp[3];
	}


// A=b*B^{dag}+c*C^{dag}
static inline void lin_comb_dag12_Su2(Su2 *restrict A, double const b, Su2 const *const restrict B, double const c, Su2 const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	__assume_aligned(&(C->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = b * B->comp[0] + c * C->comp[0];
	A->comp[1] = -b * B->comp[1] - c * C->comp[1];
	A->comp[2] = -b * B->comp[2] - c * C->comp[2];
	A->comp[3] = -b * B->comp[3] - c * C->comp[3];
	}


// random generation

// A=random
static inline void rand_matrix_Su2(Su2 *restrict A)
	{
	register double p0, p1, p2, p3, r2;

	do
		{
		p0 = 1.0 - 2.0 * casuale();
		p1 = 1.0 - 2.0 * casuale();
		p2 = 1.0 - 2.0 * casuale();
		p3 = 1.0 - 2.0 * casuale();

		r2 = p0 * p0 + p1 * p1 + p2 * p2 + p3 * p3;
		} while(r2 > 1.0);

	double const invr = 1.0 / sqrt(r2);

	A->comp[0] = p0 * invr;
	A->comp[1] = p1 * invr;
	A->comp[2] = p2 * invr;
	A->comp[3] = p3 * invr;
	}


// A=random with given p0 (used in the update)
static inline void rand_matrix_p0_Su2(double const p0, Su2 *restrict A)
	{
	register double p1, p2, p3, r2;

	do
		{
		p1 = 1.0 - 2.0 * casuale();
		p2 = 1.0 - 2.0 * casuale();
		p3 = 1.0 - 2.0 * casuale();

		r2 = p1 * p1 + p2 * p2 + p3 * p3;
		} while(r2 > 1.0);

	double const scale = sqrt((1.0 - p0 * p0) / r2);

	A->comp[0] = p0;
	A->comp[1] = p1 * scale;
	A->comp[2] = p2 * scale;
	A->comp[3] = p3 * scale;
	}


// norms and traces

// \sqrt{Det[A]}
static inline double sqrtdet_Su2(Su2 const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	return sqrt(A->comp[0] * A->comp[0] + A->comp[1] * A->comp[1] + A->comp[2] * A->comp[2] + A->comp[3] * A->comp[3]);
	}


// l2 norm
static inline double norm_Su2(Su2 const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	return sqrtdet_Su2(A);
	}


// ReTr[A]/N
static inline double retr_Su2(Su2 const *const restrict A)
	{
	return A->comp[0];
	}


// ImTr[A]/N
static inline double imtr_Su2(Su2 const *const restrict A)
	{
	(void) A; // to suppress compilation warning of unused variable
	return 0.0;
	}


// norm(A - B) / (1/2 * \sqrt{norm(A - 1)**2 + norm(B - 1)**2})
static inline double relative_dist_Su2(Su2 const *const restrict A, Su2 const *const restrict B)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	double const aux_ApB = 1.0 - 0.5 * (A->comp[0] + B->comp[0]);
	double const aux_AxB = 1.0 - (A->comp[0] * B->comp[0] + A->comp[1] * B->comp[1] + A->comp[2] * B->comp[2] + A->comp[3] * B->comp[3]);
	double const check = 2.83 * MIN_VALUE;
	if(aux_ApB <= check || aux_AxB <= check)
		{
		return 0.0;
		}
	return sqrt(aux_AxB / aux_ApB);
	}


// unitarization and exponentiation

// unitarize A
static inline void unitarize_Su2(Su2 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	double p = A->comp[0] * A->comp[0] + A->comp[1] * A->comp[1] + A->comp[2] * A->comp[2] + A->comp[3] * A->comp[3];
	p = 1.0 / sqrt(p);

	A->comp[0] *= p;
	A->comp[1] *= p;
	A->comp[2] *= p;
	A->comp[3] *= p;
	}


// traceless antihermitian part
static inline void ta_Su2(Su2 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp[0] = 0;
	}


// exponential of the traceless antihermitian part
static inline void taexp_Su2(Su2 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	// comp[0] is neglected since we consider the ta part
	double norm = A->comp[1] * A->comp[1];
	norm += A->comp[2] * A->comp[2];
	norm += A->comp[3] * A->comp[3];
	norm = sqrt(norm);

	double const v1 = A->comp[1] / norm;
	double const v2 = A->comp[2] / norm;
	double const v3 = A->comp[3] / norm;

	double const s = sin(norm);

	A->comp[0] = cos(norm);
	A->comp[1] = v1 * s;
	A->comp[2] = v2 * s;
	A->comp[3] = v3 * s;
	}


// I/O operations

// print on screen
static inline void print_on_screen_Su2(Su2 const *const restrict A)
	{
	//fprintf(stdout, "% 10.4e % 10.4e % 10.4e % 10.4e\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
	double complex f[2][2];
	f[0][0] = A->comp[0] + A->comp[3] * I;
	f[0][1] = A->comp[2] + A->comp[1] * I;
	f[1][0] = -A->comp[2] + A->comp[1] * I;
	f[1][1] = A->comp[0] - A->comp[3] * I;
	for(int i = 0; i < 2; i++)
		{
		for(int j = 0; j < 2; j++)
			{
			fprintf(stdout, "(% 5.3f % 5.3f) ", creal(f[i][j]), cimag(f[i][j]));
			}
		fprintf(stdout, "\n");
		}
	fprintf(stdout, "\n");
	}


// print on file
static inline void print_on_file_Su2(FILE *fp, Su2 const *const restrict A)
	{
	int const err = fprintf(fp, "% 18.12e % 18.12e % 18.12e % 18.12e\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
	REQUIRE(err >= 0, "failed to write an SU(2) matrix on a file");
	}


// print on binary file without changing endiannes
static inline void print_on_binary_file_noswap_Su2(FILE *fp, Su2 const *const restrict A)
	{
	size_t err = 0;
	err += fwrite(&A->comp[0], sizeof(double), 1, fp);
	err += fwrite(&A->comp[1], sizeof(double), 1, fp);
	err += fwrite(&A->comp[2], sizeof(double), 1, fp);
	err += fwrite(&A->comp[3], sizeof(double), 1, fp);
	REQUIRE(err == 4, "failed to write an SU(2) matrix on a file in binary mode");
	}


// print on binary file changing endiannes
static inline void print_on_binary_file_swap_Su2(FILE *fp, Su2 const *const restrict A)
	{
	double tmp;
	size_t err = 0;

	tmp = A->comp[0];
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);

	tmp = A->comp[1];
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);

	tmp = A->comp[2];
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);

	tmp = A->comp[3];
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);
	REQUIRE(err == 4, "failed to write an SU(2) matrix on a file in binary mode with swapped endianness");
	}


// print on binary file in big endian format
static inline void print_on_binary_file_bigen_Su2(FILE *fp, Su2 const *const restrict A)
	{
	if(endian() == 0) // little endian machine
		{
		print_on_binary_file_swap_Su2(fp, A);
		}
	else
		{
		print_on_binary_file_noswap_Su2(fp, A);
		}
	}


// read from file
static inline void read_from_file_Su2(FILE *fp, Su2 *restrict A)
	{
	int const err = fscanf(fp, "%lg %lg %lg %lg", &A->comp[0], &A->comp[1], &A->comp[2], &A->comp[3]);
	REQUIRE(err == 4, "failed to read an SU(2) matrix from a file");
	}


// read from binary file without changing endiannes
static inline void read_from_binary_file_noswap_Su2(FILE *fp, Su2 *restrict A)
	{
	size_t err = 0;
	err += fread(&A->comp[0], sizeof(double), 1, fp);
	err += fread(&A->comp[1], sizeof(double), 1, fp);
	err += fread(&A->comp[2], sizeof(double), 1, fp);
	err += fread(&A->comp[3], sizeof(double), 1, fp);
	REQUIRE(err == 4, "failed to read an SU(2) matrix from a file in binary mode");
	}


// read from binary file changing endiannes
static inline void read_from_binary_file_swap_Su2(FILE *fp, Su2 *restrict A)
	{
	size_t err = 0;
	err += fread(&A->comp[0], sizeof(double), 1, fp);
	err += fread(&A->comp[1], sizeof(double), 1, fp);
	err += fread(&A->comp[2], sizeof(double), 1, fp);
	err += fread(&A->comp[3], sizeof(double), 1, fp);
	REQUIRE(err == 4, "failed to read an SU(2) matrix from a file in binary mode with swapped endianness");

	SwapBytesDouble((void *) &A->comp[0]);
	SwapBytesDouble((void *) &A->comp[1]);
	SwapBytesDouble((void *) &A->comp[2]);
	SwapBytesDouble((void *) &A->comp[3]);
	}


// read from binary file written in big endian
static inline void read_from_binary_file_bigen_Su2(FILE *fp, Su2 *restrict A)
	{
	if(endian() == 0) // little endian machine
		{
		read_from_binary_file_swap_Su2(fp, A);
		}
	else
		{
		read_from_binary_file_noswap_Su2(fp, A);
		}
	}


#endif
