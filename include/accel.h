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

// vector_length (threads per gang) of the OpenACC compute kernels, by family; override at build time with
// -DACC_VL_UPD=<n> etc. The simple copy kernels keep the compiler default.
//   ACC_VL_UPD   update, cooling, clovers, accept/restore with unitarization (no meaningful effect with no reductions)
//   ACC_VL_FLOW  steps of the gradient flow
//   ACC_VL_MEAS  measures and the reductions that drive decisions (plaquette sums of the swaps, multicanonic delta Q)
// Defaults (from benchmarks on A100): 32 for the update (small rectangles of the hierarchical update spread over more
// SMs; full-lattice kernels unchanged), 128 for flow and measures (32 brings nothing there).
#ifndef ACC_VL_UPD
#define ACC_VL_UPD 32
#endif
#ifndef ACC_VL_FLOW
#define ACC_VL_FLOW 128
#endif
#ifndef ACC_VL_MEAS
#define ACC_VL_MEAS 128
#endif

#endif
