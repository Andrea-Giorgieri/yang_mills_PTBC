#ifndef ACCEL_H
#define ACCEL_H

#include "../config.h"

// macros to compile the same source for CPU (OpenMP) and GPU (OpenACC)

// stringify the argument of _Pragma
#define ACC_PRAGMA(x) _Pragma(#x)

// mark a function as callable from OpenACC compute regions (sequential inside each thread)
// usage: put ACC_ROUTINE_SEQ on the line before the definition and before every prototype
#ifdef OPENACC_MODE
#define ACC_ROUTINE_SEQ ACC_PRAGMA(acc routine seq)
#else
#define ACC_ROUTINE_SEQ
#endif

#endif
