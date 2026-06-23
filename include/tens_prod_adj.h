#ifndef TENS_PROD_ADJ_H
#define TENS_PROD_ADJ_H

#include"macro.h"

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>

// see Luscher Weisz JHEP 0109 p. 010 (2001) (hep-lat/0108014)
typedef struct TensProdAdj
	{
	#if NCOLOR != 1
	double comp[NCOLOR * NCOLOR - 1][NCOLOR * NCOLOR - 1][NCOLOR * NCOLOR - 1][NCOLOR * NCOLOR - 1] __attribute__((aligned(DOUBLE_ALIGN)));
	#else // this will never be used, it is defined just to avoid warnings
	double comp[1][1][1][1] __attribute__((aligned(DOUBLE_ALIGN)));
	#endif
	} TensProdAdj;


// initialize to zero
inline void zero_TensProdAdj(TensProdAdj *restrict A)
	{
	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR * NCOLOR - 1; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR * NCOLOR - 1; i3++)
					{
					A->comp[i0][i1][i2][i3] = 0.0;
					}
				}
			}
		}
	}


// initialize to one
inline void one_TensProdAdj(TensProdAdj *restrict A)
	{
	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR * NCOLOR - 1; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR * NCOLOR - 1; i3++)
					{
					A->comp[i0][i1][i2][i3] = 0.0;
					}
				}
			}
		}

	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			A->comp[i0][i0][i1][i1] = 1.0;
			}
		}
	}


// A=B
inline void equal_TensProdAdj(TensProdAdj *restrict A, TensProdAdj const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR * NCOLOR - 1; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR * NCOLOR - 1; i3++)
					{
					A->comp[i0][i1][i2][i3] = B->comp[i0][i1][i2][i3];
					}
				}
			}
		}
	}


// A*=r
inline void times_equal_real_TensProdAdj(TensProdAdj *restrict A, double r)
	{
	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR * NCOLOR - 1; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR * NCOLOR - 1; i3++)
					{
					A->comp[i0][i1][i2][i3] *= r;
					}
				}
			}
		}
	}


// A+=B
inline void plus_equal_TensProdAdj(TensProdAdj *restrict A, TensProdAdj const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			for(int i2 = 0; i2 < NCOLOR * NCOLOR - 1; i2++)
				{
				for(int i3 = 0; i3 < NCOLOR * NCOLOR - 1; i3++)
					{
					A->comp[i0][i1][i2][i3] += B->comp[i0][i1][i2][i3];
					}
				}
			}
		}
	}


// A=B*C
inline void times_TensProdAdj(TensProdAdj *restrict A,
                              TensProdAdj const *const restrict B,
                              TensProdAdj const *const restrict C)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	ASSERT(A != C, "the same pointer is used twice");
	ASSERT(B != C, "the same pointer is used twice");
	#endif

	for(int i = 0; i < NCOLOR * NCOLOR - 1; i++)
		{
		for(int j = 0; j < NCOLOR * NCOLOR - 1; j++)
			{
			for(int k = 0; k < NCOLOR * NCOLOR - 1; k++)
				{
				for(int l = 0; l < NCOLOR * NCOLOR - 1; l++)
					{
					double sum = 0.0;
					for(int a = 0; a < NCOLOR * NCOLOR - 1; a++)
						{
						for(int b = 0; b < NCOLOR * NCOLOR - 1; b++)
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
inline void times_equal_TensProdAdj(TensProdAdj *restrict A, TensProdAdj const *const restrict B)
	{
	#ifdef DEBUG
	ASSERT(A != B, "the same pointer is used twice");
	#endif

	TensProdAdj tmp __attribute__((aligned(DOUBLE_ALIGN)));

	equal_TensProdAdj(&tmp, A);
	times_TensProdAdj(A, &tmp, B);
	}


inline double retr_TensProdAdj(TensProdAdj const *const restrict A)
	{
	double tr = 0.0;
	for(int i0 = 0; i0 < NCOLOR * NCOLOR - 1; i0++)
		{
		for(int i1 = 0; i1 < NCOLOR * NCOLOR - 1; i1++)
			{
			tr += A->comp[i0][i0][i1][i1];
			}
		}

	#if NCOLOR != 1
	tr /= (NCOLOR * NCOLOR - 1) * (NCOLOR * NCOLOR - 1);
	#endif

	return tr;
	}


inline double imtr_TensProdAdj(TensProdAdj const *const restrict A)
	{
	(void) A;

	return 0;
	}


#endif



