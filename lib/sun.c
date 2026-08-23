#ifndef SUN_C
#define SUN_C

#include "../include/macro.h"
#include "../include/aligncheck.h"
#include "../include/endianness.h"
#include "../include/random.h"
#include "../include/sun.h"
#include "../include/sun_upd.h"
#include "../include/tens_prod.h"
#include "../include/tens_prod_adj.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// A=1
void one_SuN(SuN *A);


// A=0
void zero_SuN(SuN *A);


// A=B
void equal_SuN(SuN *A, SuN const *const B);


// A=B^{dag}
void equal_dag_SuN(SuN *A, SuN const *const B);


// A+=B
void plus_equal_SuN(SuN *A, SuN const *const B);


// A+=B^{dag}
void plus_equal_dag_SuN(SuN *A, SuN const *const B);


// A-=B
void minus_equal_SuN(SuN *A, SuN const *const B);


// A-=(r*B)
void minus_equal_times_real_SuN(SuN *A, SuN const *const B, double r);


// A-=B^{dag}
void minus_equal_dag_SuN(SuN *A, SuN const *const B);


// A=b*B+c*C
void lin_comb_SuN(SuN *A, double b, SuN const *const B, double c, SuN const *const C);


// A=b*B^{dag}+c*C
void lin_comb_dag1_SuN(SuN *A, double b, SuN const *const B, double c, SuN const *const C);


// A=b*B+c*C^{dag}
void lin_comb_dag2_SuN(SuN *A, double b, SuN const *const B, double c, SuN const *const C);


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_SuN(SuN *A, double b, SuN const *const B, double c, SuN const *const C);


// A*=r
void times_equal_real_SuN(SuN *A, double r);


// A*=r
void times_equal_complex_SuN(SuN *A, double complex r);


// A*=B
void times_equal_SuN(SuN *A, SuN const *const B);


// A*=B^{dag}
void times_equal_dag_SuN(SuN *A, SuN const *const B);


// A=B*C
void times_SuN(SuN *A, SuN const *const B, SuN const *const C);


// A=B^{dag}*C
void times_dag1_SuN(SuN *A, SuN const *const B, SuN const *const C);


// A=B*C^{dag}
void times_dag2_SuN(SuN *A, SuN const *const B, SuN const *const C);


// A=B^{dag}*C^{dag}
void times_dag12_SuN(SuN *A, SuN const *const B, SuN const *const C);


// SU(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
void rand_matrix_SuN(SuN *A)
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
				p0 = 1.0 - 2.0 * casuale();
				p1 = 1.0 - 2.0 * casuale();
				p2 = 1.0 - 2.0 * casuale();
				p3 = 1.0 - 2.0 * casuale();

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


// generate a matrix in the algebra of SuN with gaussian
// random components in the base T_i such that Tr(T_iT_j)=delta_{ij}
void rand_algebra_gauss_matrix_SuN(SuN *A)
	{
	#if NCOLOR == 1
	(void) A; // just to avoid warnings
	#else
	double d1, d2, dd[NCOLOR - 1];
	const double factor = sqrt(2.0 / (double) (NCOLOR * NCOLOR - NCOLOR));

	zero_SuN(A);

	// out of diagonal elements
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			gauss2(&d1, &d2);
			A->comp[m(i, j)] = d1 - d2 * I;
			A->comp[m(j, i)] = d1 + d2 * I;
			}
		}

	// random numbers to be used in the diagonal
	for(int i = 0; i < NCOLOR - 1; i++)
		{
		dd[i] = gauss1();
		}

	// diagonal
	#if NCOLOR == 2
	A->comp[m(0, 0)] = dd[0];
	A->comp[m(1, 1)] = -dd[0];
	#else
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


// l2 norm of the matrix
double norm_SuN(SuN const *const A);


// real part of the trace /N
double retr_SuN(SuN const *const A);


// imaginary part of the trace /N
double imtr_SuN(SuN const *const A);


// carg() of the trace
inline double argtr_SuN(SuN const *const A);


// trace of A * B^{dag} / N
inline double complex tr_times_dag_SuN(SuN const *const A, SuN const *const B);


// relative distance between matrices
double relative_dist_SuN(SuN const *const A, SuN const *const B);


// LU decomposition with partial pivoting
// from Numerical Recipes in C, pag 46
void LU_SuN(SuN const *const restrict A, SuN *restrict res, int *restrict sign)
	{
	int i, j, k;
	double big, temp;
	double complex sum, dum;
	double vv[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

	int imax = 0;
	equal_SuN(res, A);

	(*sign) = 1;
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
				sum -= (res->comp[m(i, k)]) * (res->comp[m(k, j)]);
				}
			res->comp[m(i, j)] = sum;
			}

		big = 0.0;
		for(i = j; i < NCOLOR; i++)
			{
			sum = res->comp[m(i, j)];
			for(k = 0; k < j; k++)
				{
				sum -= (res->comp[m(i, k)]) * (res->comp[m(k, j)]);
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
			(*sign) *= (-1);
			vv[imax] = vv[j];
			}

		if(j != NCOLOR - 1)
			{
			dum = (1.0 + 0.0 * I) / (res->comp[m(j, j)]);
			for(i = j + 1; i < NCOLOR; i++)
				{
				(res->comp[m(i, j)]) *= dum;
				}
			}
		}
	}


// determinant
complex double det_SuN(SuN const *const A);


// gives 0 if the matrix is in SU(N) and 1 otherwise
int scheck_SuN(SuN const *const restrict A)
	{
	int res = 0;

	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double complex aux = 0.0 + 0.0 * I;
			for(int k = 0; k < NCOLOR; k++)
				{
				aux += (A->comp[m(i, k)]) * conj(A->comp[m(j, k)]);
				}
			if(i == j) aux -= (1.0 + 0.0 * I);
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


// sunitarize
void unitarize_SuN(SuN *restrict A)
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


// TODO: bad unitarize for testing, remove
// unitarize with heatbath at beta O(1e20) saving bad links to file fp if print_flag != 0,
// returns phase gained during unitarization (bad links if != 0)
double bad_unitarize_SuN(SuN *restrict A, double const beta, FILE *fp, int const print_flag)
	{
	double center_element = 0.0; // k such that gained phase is C = exp(-i 2pi/N k)
	SuN F;                       // F = A^{dag}, force to unitarize A by cooling
	SuN G, G_old;                // current and previous guess for unitarized A
	SuN H, H_copy, H_square;     // helpers to check convergence of unitarization
	SuN prod;

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

			// get new guess by heatbath
			single_heatbath_aux_SuN(&G, &F, beta); // maximize Tr(G*F) in large-beta limit

			// calculate the distance between old guess G_old and new guess G:
			// check = sqrt(|ReTr[(G-G_old)^2]|/N^2)
			equal_SuN(&H, &G);
			minus_equal_SuN(&H, &G_old);
			equal_SuN(&H_copy, &H);
			times_SuN(&H_square, &H, &H_copy);
			check = sqrt(fabs(retr_SuN(&H_square)) / (double) NCOLOR);
			}

		// Maximize ReTr[staple * C *link] for C \in Z(SU(N)) and update link *= C
		// \phi \equiv carg(Tr[staple * link]) => C = \exp{-i * 2\pi/N * round(\phi / (2*\pi/N))}
		equal_SuN(&prod, &F);                                   // prod=staple
		times_equal_SuN(&prod, &G);                             // prod=staple*link
		center_element = argtr_SuN(&prod);                      // center_element = phi
		center_element = round(center_element / PI2_N) * PI2_N; // round center_element to nearest center phase (PI2_N = 2*pi/N in marco.h)
		if(print_flag != 0 && fabs(center_element) > MIN_VALUE) // bad link: center_element != 0
			print_on_file_SuN(fp, A);

		// replace A with G
		equal_SuN(A, &G);
		}
	return center_element;
	}


// takes the traceless antihermitian part
void ta_SuN(SuN *A);


// exponential of the traceless antihermitian part
void taexp_SuN(SuN *A);


// return 0 if matrix is traceless antihermitian, 1 otherwise
int ta_check_SuN(const SuN *const A);


// exponential of a TA matrix
void exp_of_ta_SuN(SuN *A);


// print on screen
void print_on_screen_SuN(SuN const *const A)
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
void print_on_file_SuN(FILE *fp, SuN const *const A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			int err = fprintf(fp, "% 18.12e % 18.12e ", creal(A->comp[m(i, j)]), cimag(A->comp[m(i, j)]));
			REQUIRE(err >= 0, "failed to write an SU(N) matrix on a file");
			}
		}
	fprintf(fp, "\n");
	}


// print on binary file without changing endiannes
void print_on_binary_file_noswap_SuN(FILE *fp, SuN const *const A)
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
void print_on_binary_file_swap_SuN(FILE *fp, SuN const *const A)
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
void print_on_binary_file_bigen_SuN(FILE *fp, SuN const *const A)
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
void read_from_file_SuN(FILE *fp, SuN *A)
	{
	for(int i = 0; i < NCOLOR; i++)
		{
		for(int j = 0; j < NCOLOR; j++)
			{
			double re, im;
			int err = fscanf(fp, "%lg %lg", &re, &im);
			REQUIRE(err == 2, "failed to read the (%d, %d) component of an SU(N) matrix from file", i, j);
			A->comp[m(i, j)] = re + im * I;
			}
		}
	}


// read from binary file without changing endiannes
void read_from_binary_file_noswap_SuN(FILE *fp, SuN *A)
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

			memcpy((void *) &(A->comp[m(i, j)]), (void *) aux, sizeof(aux));
			//equivalent to A->comp[m(i,j)]=re+im*I;
			}
		}
	REQUIRE(err == 2 * NCOLOR * NCOLOR, "failed to read an SU(N) matrix from a file in binary mode");
	}


// read from binary file changing endianness
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A)
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

			memcpy((void *) &(A->comp[m(i, j)]), (void *) aux, sizeof(aux));
			// equivalent to A->comp[m(i,j)]=re+im*I;
			}
		}
	}


// read from binary file written in bigendian
void read_from_binary_file_bigen_SuN(FILE *fp, SuN *A)
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


// initialize tensor product
void TensProd_init_SuN(TensProd *TP, SuN const *const A1, SuN const *const A2);


// convert the fundamental representation matrix B to the adjoint representation matrix A
void fund_to_adj_SuN(SuNAdj *restrict A, SuN const *const restrict B);


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
void TensProdAdj_init_SuN(TensProdAdj *restrict TP, SuN const *const restrict A1, SuN const *const restrict A2);


#endif
