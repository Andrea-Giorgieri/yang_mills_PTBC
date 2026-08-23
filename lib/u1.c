#ifndef U1_C
#define U1_C

#include"../include/macro.h"
#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/random.h"
#include"../include/tens_prod.h"
#include"../include/u1.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>


// initialize
void init_U1(U1 *A, double complex vec);


// A=1
void one_U1(U1 *A);


// A=0
void zero_U1(U1 *A);


// A=B
void equal_U1(U1 *A, U1 const *const B);


// A=B^{dag}
void equal_dag_U1(U1 *A, U1 const *const B);


// A+=B
void plus_equal_U1(U1 *A, U1 const *const B);


// A+=B^{dag}
void plus_equal_dag_U1(U1 *A, U1 const *const B);


// A-=B
void minus_equal_U1(U1 *A, U1 const *const B);


// A-=(r*B)
void minus_equal_times_real_U1(U1 *A, U1 const *const B, double r);


// A-=B^{dag}
void minus_equal_dag_U1(U1 *A, U1 const *const B);


// A=b*B+c*C
void lin_comb_U1(U1 *A, double b, U1 const *const B, double c, U1 const *const C);


// A=b*B^{dag}+c*C
void lin_comb_dag1_U1(U1 *A, double b, U1 const *const B, double c, U1 const *const C);


// A=b*B+c*C^{dag}
void lin_comb_dag2_U1(U1 *A, double b, U1 const *const B, double c, U1 const *const C);


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_U1(U1 *A, double b, U1 const *const B, double c, U1 const *const C);


// A*=r
void times_equal_real_U1(U1 *A, double r);


// A*=r
void times_equal_complex_U1(U1 *A, double complex r);


// A*=B
void times_equal_U1(U1 *A, U1 const *const B);


// A*=B^{dag}
void times_equal_dag_U1(U1 *A, U1 const *const B);


// A=B*C
void times_U1(U1 *A, U1 const *const B, U1 const *const C);


// A=B^{dag}*C
void times_dag1_U1(U1 *A, U1 const *const B, U1 const *const C);


// A=B*C^{dag}
void times_dag2_U1(U1 *A, U1 const *const B, U1 const *const C);


// A=B^{dag}*C^{dag}
void times_dag12_U1(U1 *A, U1 const *const B, U1 const *const C);


// random matrix
void rand_matrix_U1(U1 *restrict A)
	{
	double const theta = 2.0 * PI * casuale();
	A->comp = cexp(I * theta);
	}


// l2 norm of the matrix
double norm_U1(U1 const *const A);


// real part of the trace
double retr_U1(U1 const *const A);


// imaginary part of the trace
double imtr_U1(U1 const *const A);


// relative distance between matrices
double relative_dist_U1(U1 const *const A, U1 const *const B);


// unitarize the matrix
void unitarize_U1(U1 *A);


// antihermitian part (NO TRACELESS!)
void ta_U1(U1 *A);


// exponential of the antihermitian part (NO TRACELESS!)
void taexp_U1(U1 *A);


// print on screen
void print_on_screen_U1(U1 const *const A)
	{
	fprintf(stdout, "% 10.4e % 10.4e\n", creal(A->comp), cimag(A->comp));
	}


// print on file
void print_on_file_U1(FILE *fp, U1 const *const A)
	{
	int err = fprintf(fp, "% 18.12e % 18.12e\n", creal(A->comp), cimag(A->comp));
	REQUIRE(err >= 0, "failed to write a U1 element on a file");
	}


// print on binary file without changing endiannes
void print_on_binary_file_noswap_U1(FILE *fp, U1 const *const A)
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
void print_on_binary_file_swap_U1(FILE *fp, U1 const *const A)
	{
	double tmp;
	size_t err = 0;

	tmp = creal(A->comp);
	SwapBytesDouble(&tmp);
	err += fwrite(&(tmp), sizeof(double), 1, fp);

	tmp = cimag(A->comp);
	SwapBytesDouble(&tmp);
	err += fwrite(&(tmp), sizeof(double), 1, fp);

	REQUIRE(err == 2, "failed to write a U1 element on a file in binary mode");
	}


// print on binary file in big endian format
void print_on_binary_file_bigen_U1(FILE *fp, const U1 *const A)
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
void read_from_file_U1(FILE *fp, U1 *A)
	{
	double re, im;

	int err = fscanf(fp, "%lg %lg", &re, &im);
	REQUIRE(err == 2, "failed to read a U1 element from a file");

	A->comp = re + im * I;
	}


// read from binary file without changing endiannes
void read_from_binary_file_noswap_U1(FILE *fp, U1 *A)
	{
	double re, im;

	size_t err = 0;
	err += fread(&re, sizeof(double), 1, fp);
	err += fread(&im, sizeof(double), 1, fp);
	REQUIRE(err == 2, "failed to read a U1 element from a file in binary mode");

	A->comp = re + I * im;
	}


// read from binary file changing endiannes
void read_from_binary_file_swap_U1(FILE *fp, U1 *A)
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
void read_from_binary_file_bigen_U1(FILE *fp, U1 *A)
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


// initialize tensor product
void TensProd_init_U1(TensProd *TP, U1 const *const A1, U1 const *const A2);


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
void TensProdAdj_init_U1(TensProdAdj *restrict TP, U1 const *const restrict A1, U1 const *const restrict A2);


#endif
