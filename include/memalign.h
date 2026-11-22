#ifndef MEMALIGN_H
#define MEMALIGN_H

#include "macro.h"

#include <complex.h>

#include "gauge_conf.h"
#include "geometry.h"
#include "tens_prod.h"

void allocate_array_int(int **const array, long size, char const *const file, int const line);

void allocate_array_int_pointer(int ***const array, long size, char const *const file, int const line);

void allocate_array_long(long **const array, long size, char const *const file, int const line);

void allocate_array_long_pointer(long ***const array, long size, char const *const file, int const line);

void allocate_array_long_pointer_pointer(long ****const array, long size, char const *const file, int line);

void allocate_array_Rectangle(Rectangle **const array, long size, char const *const file, int line);

void allocate_array_Rectangle_pointer(Rectangle ***const array, long size, char const *const file, int line);

void allocate_array_double(double **const array, long size, char const *const file, int const line);

void allocate_array_double_pointer(double ***const array, long size, char const *const file, int const line);

void allocate_array_double_complex(double complex **const array, long size, char const *const file, int const line);

void allocate_array_double_complex_pointer(double complex ***const array, long size, char const *const file, int line);

void allocate_array_GAUGE_GROUP(GAUGE_GROUP **const array, long size, char const *const file, int line);

void allocate_array_GAUGE_GROUP_pointer(GAUGE_GROUP ***const array, long size, char const *const file, int line);

void allocate_array_GAUGE_GROUP_pointer_pointer(GAUGE_GROUP ****const array, long size, char const *const file, int line);

void allocate_array_TensProd(TensProd **const array, long size, char const *const file, int line);

void allocate_array_TensProd_pointer(TensProd ***const array, long size, char const *const file, int line);

void allocate_array_TensProd_pointer_pointer(TensProd ****const array, long size, char const *const file, int line);

void allocate_array_Gauge_Conf(Gauge_Conf **const array, long size, char const *const file, int line);

void allocate_array_Meas_Utils(Meas_Utils **const array, long size, char const *const file, int line);

void allocate_array_rng(rng_t **const array, long size, char const *const file, int line);


// contiguous multidimensional arrays: one block for the elements and tables of pointers into the block
// (syntax array[i][j] unchanged); the block is array[0] (array[0][0] for 3d)

// [n1][n2]: array[i] = block + i * n2
void allocate_array_GAUGE_GROUP_2d(GAUGE_GROUP ***const array, long n1, long n2, char const *const file, int line);

void allocate_array_double_complex_2d(double complex ***const array, long n1, long n2, char const *const file, int line);

void allocate_array_long_2d(long ***const array, long n1, long n2, char const *const file, int line);

void allocate_array_int_2d(int ***const array, long n1, long n2, char const *const file, int line);

// [n1][n2][n3]: array[i][j] = block + (i * n2 + j) * n3
void allocate_array_GAUGE_GROUP_3d(GAUGE_GROUP ****const array, long n1, long n2, long n3, char const *const file, int line);

// free the block and the tables of the arrays allocated with the functions above
void free_array_GAUGE_GROUP_2d(GAUGE_GROUP **array);

void free_array_double_complex_2d(double complex **array);

void free_array_long_2d(long **array);

void free_array_int_2d(int **array);

void free_array_GAUGE_GROUP_3d(GAUGE_GROUP ***array);

#endif
