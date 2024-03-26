#ifndef MEMALIGN_C
#define MEMALIGN_C	

#include"../include/macro.h"
#include"../include/memalign.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"

#include<math.h>
#include<complex.h>
#include<stdio.h>
#include<stdlib.h>


void allocate_array_int(int **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)INT_ALIGN, (size_t)size*sizeof(int)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ints! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_long(long **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)INT_ALIGN, (size_t)size*sizeof(long)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of longs! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_long_pointer(long ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)INT_ALIGN, (size_t)size*sizeof(long*)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to long! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_Rectangle(Rectangle **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)INT_ALIGN, (size_t)size*sizeof(Rectangle)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of Rectangles! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_double(double **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(double)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_double_pointer(double ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(double*)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to doubles! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_double_complex(double complex **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(double complex)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of complex doubles! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_double_complex_pointer(double complex ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void **)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(double complex *)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to complex double! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_GAUGE_GROUP(GAUGE_GROUP **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(GAUGE_GROUP)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of GAUGE_GROUPs! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_GAUGE_GROUP_pointer(GAUGE_GROUP ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(GAUGE_GROUP*)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to GAUGE_GROUP! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_GAUGE_GROUP_pointer_pointer(GAUGE_GROUP ****const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(GAUGE_GROUP**)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to ptrs to GAUGE_GROUP! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProd(TensProd **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProd)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of TensProds! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProd_pointer(TensProd ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProd*)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to TensProd! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProd_pointer_pointer(TensProd ****const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProd**)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to ptrs to TensProd! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProdAdj(TensProdAdj **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProdAdj)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of TensProdAdjs! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProdAdj_pointer(TensProdAdj ***const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProdAdj*)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to TensProdAdj! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_TensProdAdj_pointer_pointer(TensProdAdj ****const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(TensProdAdj**)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of ptrs to ptrs to TensProdAdj! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_Gauge_Conf(Gauge_Conf **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(Gauge_Conf)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of Gauge_Confs! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}

void allocate_array_Meas_Utils(Meas_Utils **const array, long size, char const * const file, int line)
	{
	if(posix_memalign((void**)array, (size_t)DOUBLE_ALIGN, (size_t)size*sizeof(Meas_Utils)) != 0)
		{
		fprintf(stderr, "Problems allocating an array of Meas_Utils! (%s, %d)\n", file, line);
		exit(EXIT_FAILURE);
		}
	}
#endif

