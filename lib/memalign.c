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

#endif
