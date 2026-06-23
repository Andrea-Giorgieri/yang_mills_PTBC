#ifndef TENS_PROD_H
#define TENS_PROD_H

#include"macro.h"

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>

// see Luscher Weisz JHEP 0109 p. 010 (2001)   (hep-lat/0108014)
typedef struct TensProd
	{
	double complex comp[NCOLOR][NCOLOR][NCOLOR][NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
	} TensProd;


// initialize to zero
inline void zero_TensProd(TensProd *restrict A)
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
inline void one_TensProd(TensProd *restrict A)
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


// A=B
inline void equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
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
inline void times_equal_real_TensProd(TensProd *restrict A, double r)
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
inline void times_equal_complex_TensProd(TensProd *restrict A, double complex r)
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
inline void plus_equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
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
inline void times_TensProd(TensProd *restrict A,
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
inline void times_equal_TensProd(TensProd *restrict A, TensProd const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	TensProd tmp __attribute__((aligned(DOUBLE_ALIGN)));

	equal_TensProd(&tmp, A);
	times_TensProd(A, &tmp, B);
	}


inline double retr_TensProd(TensProd const *const restrict A)
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


inline double imtr_TensProd(TensProd const *const restrict A)
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


void print_on_screen_TensProd(TensProd const *const A);

void print_on_file_TensProd(FILE *fp, TensProd const *const A);

void print_on_binary_file_noswap_TensProd(FILE *fp, TensProd const *const A);

void print_on_binary_file_swap_TensProd(FILE *fp, TensProd const *const A);

void print_on_binary_file_bigen_TensProd(FILE *fp, TensProd const *const A);

void read_from_file_TensProd(FILE *fp, TensProd *A);

void read_from_binary_file_noswap_TensProd(FILE *fp, TensProd *A);

void read_from_binary_file_swap_TensProd(FILE *fp, TensProd *A);

void read_from_binary_file_bigen_TensProd(FILE *fp, TensProd *A);

#endif // TENS_PROD_H

