#ifndef TENS_PROD_H
#define TENS_PROD_H

#include "macro.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "endianness.h"
#include "gauge_group.h"


// see Luscher Weisz JHEP 0109 p. 010 (2001)   (hep-lat/0108014)
typedef struct TensProd
	{
	double complex comp[NCOLOR][NCOLOR][NCOLOR][NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
	} TensProd;


// initialize to zero
static inline void zero_TensProd(TensProd *restrict A)
	{
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] = 0.0 + 0.0 * I;
					}
				}
			}
		}
	}


// initialize to one
static inline void one_TensProd(TensProd *restrict A)
	{
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] = 0.0 + 0.0 * I;
					}
				}
			}
		}

	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			A->comp[i0][i0][i1][i1] = 1.0 + 0.0 * I;
			}
		}
	}


// initialize from two GAUGE_GROUPs
static inline void TensProd_init(TensProd *restrict TP, GAUGE_GROUP const *const restrict A1, GAUGE_GROUP const *const restrict A2)
	{
	#ifdef DEBUG
	ASSERT(A1 != A2, "the same pointer is used twice");
	#endif

	#ifdef __INTEL_COMPILER
	__assume_aligned(&(A1->comp), DOUBLE_ALIGN);
	__assume_aligned(&(A2->comp), DOUBLE_ALIGN);
	__assume_aligned(&(TP->comp), DOUBLE_ALIGN);
	#endif

	#if NCOLOR == 1
	TP->comp[0][0][0][0] = conj(A1->comp) * A2->comp;
	#else

	#if NCOLOR == 2
	// reconstruct the complex form of the matrices
	double complex aux1[4] __attribute__((aligned(DOUBLE_ALIGN)));
	double complex aux2[4] __attribute__((aligned(DOUBLE_ALIGN)));

	aux1[m(0, 0)] = A1->comp[0] + I * A1->comp[3];
	aux1[m(0, 1)] = A1->comp[2] + I * A1->comp[1];
	aux1[m(1, 0)] = -A1->comp[2] + I * A1->comp[1];
	aux1[m(1, 1)] = A1->comp[0] - I * A1->comp[3];

	aux2[m(0, 0)] = A2->comp[0] + I * A2->comp[3];
	aux2[m(0, 1)] = A2->comp[2] + I * A2->comp[1];
	aux2[m(1, 0)] = -A2->comp[2] + I * A2->comp[1];
	aux2[m(1, 1)] = A2->comp[0] - I * A2->comp[3];
	#endif

	for(int i = 0; i < 2; i++)
		{
		for(int j = 0; j < 2; j++)
			{
			for(int k = 0; k < 2; k++)
				{
				for(int l = 0; l < 2; l++)
					{
					#if NCOLOR == 2
					TP->comp[i][j][k][l] = conj(aux1[m(i, j)]) * aux2[m(k, l)];
					#else
					TP->comp[i][j][k][l] = conj(A1->comp[m(i, j)]) * A2->comp[m(k, l)];
					#endif
					}
				}
			}
		}
	#endif
	}


// A=B
static inline void equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] = B->comp[i0][i1][i2][i3];
					}
				}
			}
		}
	}


// A*=r
static inline void times_equal_real_TensProd(TensProd *restrict A, double const r)
	{
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] *= r;
					}
				}
			}
		}
	}


// A*=r
static inline void times_equal_complex_TensProd(TensProd *restrict A, double complex const r)
	{
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] *= r;
					}
				}
			}
		}
	}


// A+=B
static inline void plus_equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR; i3++)
					{
					A->comp[i0][i1][i2][i3] += B->comp[i0][i1][i2][i3];
					}
				}
			}
		}
	}


// A=B*C
static inline void times_TensProd(TensProd *restrict A,
                           TensProd const *const restrict B,
                           TensProd const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double complex sum = 0.0 + 0.0 * I;
					for(int a = 0; a < NCOLOR; a++)
						{
						for(int b = 0; b < NCOLOR; b++)
							{
							sum += B->comp[i][a][k][b] * C->comp[a][j][b][l];
							}
						}
					A->comp[i][j][k][l] = sum;
					}
				}
			}
		}
	}


// A*=B
static inline void times_equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	TensProd tmp __attribute__((aligned(DOUBLE_ALIGN)));

	equal_TensProd(&tmp, A);
	times_TensProd(A, &tmp, B);
	}


// ReTr[A]/N^2
static inline double retr_TensProd(TensProd const *const restrict A)
	{
	double complex tr = 0.0;
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			tr += A->comp[i0][i0][i1][i1];
			}
		}
	return creal(tr) / (NCOLOR * NCOLOR);
	}


// ImTr[A]/N^2
static inline double imtr_TensProd(TensProd const *const restrict A)
	{
	double complex tr = 0.0;
	for(int i0 = 0; i0 < NCOLOR; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR; i1++)
			{
			tr += A->comp[i0][i0][i1][i1];
			}
		}
	return cimag(tr) / (NCOLOR * NCOLOR);
	}


// I/O operations

static inline void print_on_screen_TensProd(TensProd const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					printf("%.16lf %.16lf ", creal(A->comp[i][j][k][l]), cimag(A->comp[i][j][k][l]));
					}
				}
			}
		}
	printf("\n");
	}


static inline void print_on_file_TensProd(FILE *fp, TensProd const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					int err = fprintf(fp, "%.16lf %.16lf", creal(A->comp[i][j][k][l]), cimag(A->comp[i][j][k][l]));
					REQUIRE(err == 2, "failed to write a TensProd on a file");
					}
				}
			}
		}
	fprintf(fp, "\n");
	}


static inline void print_on_binary_file_noswap_TensProd(FILE *fp, TensProd const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double re = creal(A->comp[i][j][k][l]);
					double im = cimag(A->comp[i][j][k][l]);

					size_t err = fwrite(&re, sizeof(double), 1, fp);
					REQUIRE(err == 1, "failed to write a TensProd on a file in binary mode");
					err = fwrite(&im, sizeof(double), 1, fp);
					REQUIRE(err == 1, "failed to write a TensProd on a file in binary mode");
					}
				}
			}
		}
	}


static inline void print_on_binary_file_swap_TensProd(FILE *fp, TensProd const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double re = creal(A->comp[i][j][k][l]);
					double im = cimag(A->comp[i][j][k][l]);

					SwapBytesDouble((void *) &re);
					SwapBytesDouble((void *) &im);

					size_t err = fwrite(&re, sizeof(double), 1, fp);
					REQUIRE(err == 1, "failed to write a TensProd on a file in binary mode");
					err = fwrite(&im, sizeof(double), 1, fp);
					REQUIRE(err == 1, "failed to write a TensProd on a file in binary mode");
					}
				}
			}
		}
	}


static inline void print_on_binary_file_bigen_TensProd(FILE *fp, TensProd const *const A)
	{
	if(endian() == 0) // little endian machine
		{
		print_on_binary_file_swap_TensProd(fp, A);
		}
	else
		{
		print_on_binary_file_noswap_TensProd(fp, A);
		}
	}


static inline void read_from_file_TensProd(FILE *fp, TensProd *A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double re, im;
					int err = fscanf(fp, "%lg %lg", &re, &im);
					REQUIRE(err == 2, "failed to read a TensProd from a file");
					A->comp[i][j][k][l] = re + I * im;
					}
				}
			}
		}
	}


static inline void read_from_binary_file_noswap_TensProd(FILE *fp, TensProd *A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double re, im;
					double aux[2];

					size_t err = 0;
					err += fread(&re, sizeof(double), 1, fp);
					err += fread(&im, sizeof(double), 1, fp);
					REQUIRE(err == 2, "failed to read a TensProd from a file in binary mode");

					aux[0] = re;
					aux[1] = im;

					memcpy((void *) &(A->comp[i][j][k][l]), (void *) aux, sizeof(aux));
					//equivalent to A->comp[i][j][k][l]=re+im*I;
					}
				}
			}
		}
	}


static inline void read_from_binary_file_swap_TensProd(FILE *fp, TensProd *A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			for(int k = 0; k < NCOLOR; k++)
				{
				for(int l = 0; l < NCOLOR; l++)
					{
					double re, im;
					double aux[2];

					size_t err = 0;
					err += fread(&re, sizeof(double), 1, fp);
					err += fread(&im, sizeof(double), 1, fp);
					REQUIRE(err == 2, "failed to read a TensProd from a file in binary mode");

					SwapBytesDouble(&re);
					SwapBytesDouble(&im);
					aux[0] = re;
					aux[1] = im;

					memcpy((void *) &(A->comp[i][j][k][l]), (void *) aux, sizeof(aux));
					//equivalent to A->comp[i][j][k][l]=re+im*I;
					}
				}
			}
		}
	}


static inline void read_from_binary_file_bigen_TensProd(FILE *fp, TensProd *A)
	{
	if(endian() == 0) // little endian machine
		{
		read_from_binary_file_swap_TensProd(fp, A);
		}
	else
		{
		read_from_binary_file_noswap_TensProd(fp, A);
		}
	}

#endif

