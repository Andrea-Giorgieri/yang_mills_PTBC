#ifndef MEMALIGN_C
#define MEMALIGN_C

#include "../include/macro.h"
#include "../include/memalign.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/gauge_conf.h"


// aligned allocation of size elements of elem_size bytes each, exit with an error message on failure
static void aligned_alloc_or_exit(void **const array,
                                  size_t alignment,
                                  long size,
                                  size_t elem_size,
                                  char const *const type_name,
                                  char const *const file,
                                  int line)
	{
	int const err = posix_memalign(array, alignment, (size_t) size * elem_size);
	if(err != 0)
		{
		fprintf(stderr, "Problems allocating an array of %s! (%s, %d):\n", type_name, file, line);
		fprintf(stderr, "    requested %ld elements of %zu bytes\n    error %d: %s\n", size, elem_size, err, strerror(err));
		exit(EXIT_FAILURE);
		}
	}


void allocate_array_int(int **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(int), "ints", file, line);
	}

void allocate_array_int_pointer(int ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(int *), "ptrs to int", file, line);
	}

void allocate_array_long(long **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(long), "longs", file, line);
	}

void allocate_array_long_pointer(long ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(long *), "ptrs to long", file, line);
	}

void allocate_array_long_pointer_pointer(long ****const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(long **), "ptrs to ptrs to long", file, line);
	}

void allocate_array_Rectangle(Rectangle **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(Rectangle), "Rectangles", file, line);
	}

void allocate_array_Rectangle_pointer(Rectangle ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) INT_ALIGN, size, sizeof(Rectangle *), "ptrs to Rectangle", file, line);
	}

void allocate_array_double(double **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(double), "doubles", file, line);
	}

void allocate_array_double_pointer(double ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(double *), "ptrs to doubles", file, line);
	}

void allocate_array_double_complex(double complex **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(double complex), "complex doubles", file, line);
	}

void allocate_array_double_complex_pointer(double complex ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(double complex *), "ptrs to complex double", file, line);
	}

void allocate_array_GAUGE_GROUP(GAUGE_GROUP **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(GAUGE_GROUP), "GAUGE_GROUPs", file, line);
	}

void allocate_array_GAUGE_GROUP_pointer(GAUGE_GROUP ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(GAUGE_GROUP *), "ptrs to GAUGE_GROUP", file, line);
	}

void allocate_array_GAUGE_GROUP_pointer_pointer(GAUGE_GROUP ****const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(GAUGE_GROUP **), "ptrs to ptrs to GAUGE_GROUP", file, line);
	}

void allocate_array_TensProd(TensProd **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(TensProd), "TensProds", file, line);
	}

void allocate_array_TensProd_pointer(TensProd ***const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(TensProd *), "ptrs to TensProd", file, line);
	}

void allocate_array_TensProd_pointer_pointer(TensProd ****const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(TensProd **), "ptrs to ptrs to TensProd", file, line);
	}

void allocate_array_Gauge_Conf(Gauge_Conf **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(Gauge_Conf), "Gauge_Confs", file, line);
	}

void allocate_array_Meas_Utils(Meas_Utils **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(Meas_Utils), "Meas_Utils", file, line);
	}

void allocate_array_rng(rng_t **const array, long size, char const *const file, int line)
	{
	aligned_alloc_or_exit((void **) array, (size_t) DOUBLE_ALIGN, size, sizeof(rng_t), "random number generators", file, line);
	}


// contiguous multidimensional arrays: one block for the elements and tables of pointers into the block,
// so that the syntax array[i][j] (array[i][j][k]) is unchanged; the block is array[0] (array[0][0])
// and is released by the corresponding free_array_*_2d (free_array_*_3d)

// [n1][n2]: array[i] = block + i * n2
void allocate_array_GAUGE_GROUP_2d(GAUGE_GROUP ***const array, long n1, long n2, char const *const file, int line)
	{
	GAUGE_GROUP *block;
	allocate_array_GAUGE_GROUP(&block, n1 * n2, file, line);
	allocate_array_GAUGE_GROUP_pointer(array, n1, file, line);
	for(long i = 0; i < n1; i++)
		(*array)[i] = block + i * n2;
	}

void allocate_array_double_complex_2d(double complex ***const array, long n1, long n2, char const *const file, int line)
	{
	double complex *block;
	allocate_array_double_complex(&block, n1 * n2, file, line);
	allocate_array_double_complex_pointer(array, n1, file, line);
	for(long i = 0; i < n1; i++)
		(*array)[i] = block + i * n2;
	}

void allocate_array_long_2d(long ***const array, long n1, long n2, char const *const file, int line)
	{
	long *block;
	allocate_array_long(&block, n1 * n2, file, line);
	allocate_array_long_pointer(array, n1, file, line);
	for(long i = 0; i < n1; i++)
		(*array)[i] = block + i * n2;
	}

void allocate_array_int_2d(int ***const array, long n1, long n2, char const *const file, int line)
	{
	int *block;
	allocate_array_int(&block, n1 * n2, file, line);
	allocate_array_int_pointer(array, n1, file, line);
	for(long i = 0; i < n1; i++)
		(*array)[i] = block + i * n2;
	}

// [n1][n2][n3]: array[i] = rows + i * n2, array[i][j] = block + (i * n2 + j) * n3
void allocate_array_GAUGE_GROUP_3d(GAUGE_GROUP ****const array, long n1, long n2, long n3, char const *const file, int line)
	{
	GAUGE_GROUP *block;
	GAUGE_GROUP **rows;
	allocate_array_GAUGE_GROUP(&block, n1 * n2 * n3, file, line);
	allocate_array_GAUGE_GROUP_pointer(&rows, n1 * n2, file, line);
	allocate_array_GAUGE_GROUP_pointer_pointer(array, n1, file, line);
	for(long i = 0; i < n1; i++)
		{
		(*array)[i] = rows + i * n2;
		for(long j = 0; j < n2; j++)
			(*array)[i][j] = block + (i * n2 + j) * n3;
		}
	}


void free_array_GAUGE_GROUP_2d(GAUGE_GROUP **array)
	{
	free(array[0]);
	free(array);
	}

void free_array_double_complex_2d(double complex **array)
	{
	free(array[0]);
	free(array);
	}

void free_array_long_2d(long **array)
	{
	free(array[0]);
	free(array);
	}

void free_array_int_2d(int **array)
	{
	free(array[0]);
	free(array);
	}

void free_array_GAUGE_GROUP_3d(GAUGE_GROUP ***array)
	{
	free(array[0][0]);
	free(array[0]);
	free(array);
	}

#endif
