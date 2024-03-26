#ifndef MEMALIGN_H
#define MEMALIGN_H

#include"macro.h"
#include<complex.h>

#include"u1.h"
#include"su2.h"
#include"sun.h"
#include"gparam.h"
#include"geometry.h"
#include"gauge_conf.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

void allocate_array_int(int **const array, long size, char const * const file, int const line);
void allocate_array_long(long **const array, long size, char const * const file, int const line);
void allocate_array_long_pointer(long ***const array, long size, char const * const file, int const line);
void allocate_array_Rectangle(Rectangle **const array, long size, char const * const file, int line);
void allocate_array_double(double **const array, long size, char const * const file, int const line);
void allocate_array_double_pointer(double ***const array, long size, char const * const file, int const line);
void allocate_array_double_complex(double complex **const array, long size, char const * const file, int const line);
void allocate_array_double_complex_pointer(double complex ***const array, long size, char const * const file, int line);
void allocate_array_GAUGE_GROUP(GAUGE_GROUP **const array, long size, char const * const file, int line);
void allocate_array_GAUGE_GROUP_pointer(GAUGE_GROUP ***const array, long size, char const * const file, int line);
void allocate_array_GAUGE_GROUP_pointer_pointer(GAUGE_GROUP ****const array, long size, char const * const file, int line);
void allocate_array_TensProd(TensProd **const array, long size, char const * const file, int line);
void allocate_array_TensProd_pointer(TensProd ***const array, long size, char const * const file, int line);
void allocate_array_TensProd_pointer_pointer(TensProd ****const array, long size, char const * const file, int line);
void allocate_array_TensProdAdj(TensProdAdj **const array, long size, char const * const file, int line);
void allocate_array_TensProdAdj_pointer(TensProdAdj ***const array, long size, char const * const file, int line);
void allocate_array_TensProdAdj_pointer_pointer(TensProdAdj ****const array, long size, char const * const file, int line);
void allocate_array_Gauge_Conf(Gauge_Conf **const array, long size, char const * const file, int line);
void allocate_array_Meas_Utils(Meas_Utils **const array, long size, char const * const file, int line);

#endif
