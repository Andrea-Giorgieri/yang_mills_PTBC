#ifndef U1_H
#define U1_H

#include "macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "endianness.h"
#include "random.h"


typedef struct U1
	{
	double complex comp __attribute__((aligned(DOUBLE_ALIGN)));
	} U1;


// basic operations

// A=vec
static inline void init_U1(U1 *restrict A, double complex const vec)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp = vec;
	}


// A=1
static inline void one_U1(U1 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp = 1.0;
	}


// A=0
static inline void zero_U1(U1 *restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp = 0.0;
	}


// A=B
static inline void equal_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp = B->comp;
	}


// A=B^{dag}
static inline void equal_dag_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp = conj(B->comp);
	}


// additions and subtractions

// A+=B
static inline void plus_equal_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp += B->comp;
	}


// A+=B^{dag}
static inline void plus_equal_dag_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp += conj(B->comp);
	}


// A-=B
static inline void minus_equal_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp -= B->comp;
	}


// A-=(r*B)
static inline void minus_equal_times_real_U1(U1 *restrict A, U1 const *const restrict B, double const r)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp -= r * B->comp;
	}


// A-=B^{dag}
static inline void minus_equal_dag_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp -= conj(B->comp);
	}


// multiplications

// A*=r
static inline void times_equal_real_U1(U1 *restrict A, double const r)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp *= r;
	}


// A*=r
static inline void times_equal_complex_U1(U1 *restrict A, double complex const r)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	A->comp *= r;
	}


// A*=B
static inline void times_equal_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp *= B->comp;
	}


// A*=B^{dag}
static inline void times_equal_dag_U1(U1 *restrict A, U1 const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	A->comp *= conj(B->comp);
	}


// A=B*C
static inline void times_U1(U1 *restrict A,
                            U1 const *const restrict B,
                            U1 const *const restrict C)
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

	A->comp = B->comp * C->comp;
	}


// A=B^{dag}*C
static inline void times_dag1_U1(U1 *restrict A,
                                 U1 const *const restrict B,
                                 U1 const *const restrict C)
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

	A->comp = conj(B->comp) * C->comp;
	}


// A=B*C^{dag}
static inline void times_dag2_U1(U1 *restrict A,
                                 U1 const *const restrict B,
                                 U1 const *const restrict C)
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

	A->comp = B->comp * conj(C->comp);
	}


// A=B^{dag}*C^{dag}
static inline void times_dag12_U1(U1 *restrict A,
                                  U1 const *const restrict B,
                                  U1 const *const restrict C)
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

	A->comp = conj(B->comp) * conj(C->comp);
	}


// linear combinations

// A=b*B+c*C
static inline void lin_comb_U1(U1 *restrict A,
                               double const b, U1 const *const restrict B,
                               double const c, U1 const *const restrict C)
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

	A->comp = b * B->comp + c * C->comp;
	}


// A=b*B^{dag}+c*C
static inline void lin_comb_dag1_U1(U1 *restrict A,
                                    double const b, U1 const *const restrict B,
                                    double const c, U1 const *const restrict C)
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

	A->comp = b * conj(B->comp) + c * C->comp;
	}


// A=b*B+c*C^{dag}
static inline void lin_comb_dag2_U1(U1 *restrict A,
                                    double const b, U1 const *const restrict B,
                                    double const c, U1 const *const restrict C)
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

	A->comp = b * B->comp + c * conj(C->comp);
	}


// A=b*B^{dag}+c*C^{dag}
static inline void lin_comb_dag12_U1(U1 *restrict A,
                                     double const b, U1 const *const restrict B,
                                     double const c, U1 const *const restrict C)
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

	A->comp = b * conj(B->comp) + c * conj(C->comp);
	}


// random generation

// A=random
static inline void rand_matrix_U1(U1 *restrict A)
	{
	double const theta = 2.0 * PI * casuale();
	A->comp = cexp(I * theta);
	}


// norms and traces

// l2 norm
static inline double norm_U1(U1 const *const restrict A)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	#endif

	return sqrt(creal(A->comp) * creal(A->comp) + cimag(A->comp) * cimag(A->comp));
	}


// ReTr[A]/N
static inline double retr_U1(U1 const *const restrict A)
	{
	return creal(A->comp);
	}


// ReTr[A]/N
static inline double imtr_U1(U1 const *const restrict A)
	{
	return cimag(A->comp);
	}


// norm(A - B) / (1/2 * \sqrt{norm(A - 1)**2 + norm(B - 1)**2})
static inline double relative_dist_U1(U1 const *const restrict A, U1 const *const restrict B)
	{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A->comp), DOUBLE_ALIGN);
	__assume_aligned(&(B->comp), DOUBLE_ALIGN);
	#endif

	double complex const tr_AxB = A->comp * conj(B->comp);
	double complex const tr_ApB = A->comp + B->comp;
	double const aux_ApB = 1.0 - 0.5 * creal(tr_ApB);
	double const aux_AxB = 1.0 - creal(tr_AxB);
	if(aux_ApB <= 2 * MIN_VALUE || aux_AxB <= 2 * MIN_VALUE)
		{
		return 0.0;
		}
	return sqrt(aux_AxB / aux_ApB);
	}


// unitarization and exponentiation

// unitarize A
static inline void unitarize_U1(U1 *restrict A)
	{
	double const p = norm_U1(A);
	A->comp /= p;
	}


// antihermitian part (NOT TRACELESS!)
static inline void ta_U1(U1 *restrict A)
	{
	A->comp = cimag(A->comp);
	}


// exponential of the antihermitian part (NOT TRACELESS!)
static inline void taexp_U1(U1 *restrict A)
	{
	double const angle = cimag(A->comp);
	A->comp = cexp(I * angle);
	}


// I/O operations

// print on screen
static inline void print_on_screen_U1(U1 const *const A)
	{
	fprintf(stdout, "% 10.4e % 10.4e\n", creal(A->comp), cimag(A->comp));
	}


// print on file
static inline void print_on_file_U1(FILE *fp, U1 const *const A)
	{
	int const err = fprintf(fp, "% 18.12e % 18.12e\n", creal(A->comp), cimag(A->comp));
	REQUIRE(err >= 0, "failed to write a U1 element on a file");
	}


// print on binary file without changing endiannes
static inline void print_on_binary_file_noswap_U1(FILE *fp, U1 const *const A)
	{
	size_t err = 0;
	double re, im;

	re = creal(A->comp);
	im = cimag(A->comp);

	err = fwrite(&re, sizeof(double), 1, fp);
	REQUIRE(err == 1, "failed to write a U1 element on a file in binary mode");
	err = fwrite(&im, sizeof(double), 1, fp);
	REQUIRE(err == 1, "failed to write a U1 element on a file in binary mode");
	}


// print on binary file changing endiannes
static inline void print_on_binary_file_swap_U1(FILE *fp, U1 const *const A)
	{
	double tmp;
	size_t err = 0;

	tmp = creal(A->comp);
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);

	tmp = cimag(A->comp);
	SwapBytesDouble(&tmp);
	err += fwrite(&tmp, sizeof(double), 1, fp);

	REQUIRE(err == 2, "failed to write a U1 element on a file in binary mode");
	}


// print on binary file in big endian format
static inline void print_on_binary_file_bigen_U1(FILE *fp, const U1 *const A)
	{
	if(endian() == 0) // little endian machine
		{
		print_on_binary_file_swap_U1(fp, A);
		}
	else
		{
		print_on_binary_file_noswap_U1(fp, A);
		}
	}


// read from file
static inline void read_from_file_U1(FILE *fp, U1 *A)
	{
	double re, im;

	int const err = fscanf(fp, "%lg %lg", &re, &im);
	REQUIRE(err == 2, "failed to read a U1 element from a file");

	A->comp = re + im * I;
	}


// read from binary file without changing endiannes
static inline void read_from_binary_file_noswap_U1(FILE *fp, U1 *A)
	{
	double re, im;

	size_t err = 0;
	err += fread(&re, sizeof(double), 1, fp);
	err += fread(&im, sizeof(double), 1, fp);
	REQUIRE(err == 2, "failed to read a U1 element from a file in binary mode");

	A->comp = re + I * im;
	}


// read from binary file changing endiannes
static inline void read_from_binary_file_swap_U1(FILE *fp, U1 *A)
	{
	double re, im;
	size_t err = 0;

	err += fread(&re, sizeof(double), 1, fp);
	err += fread(&im, sizeof(double), 1, fp);
	REQUIRE(err == 2, "failed to read a U1 element from a file in binary mode");

	SwapBytesDouble((void *) &re);
	SwapBytesDouble((void *) &im);

	A->comp = re + I * im;
	}


// read from binary file written in big endian
static inline void read_from_binary_file_bigen_U1(FILE *fp, U1 *A)
	{
	if(endian() == 0) // little endian machine
		{
		read_from_binary_file_swap_U1(fp, A);
		}
	else
		{
		read_from_binary_file_noswap_U1(fp, A);
		}
	}


#endif
