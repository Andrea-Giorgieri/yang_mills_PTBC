#!/bin/bash

set -euo pipefail

# ========================================================================
# Build configuration
# ========================================================================

NCOLORS=3     # Number of colors
STDIM=4       # Spacetime dimension

NTHREADS=112  # Number of OpenMP threads (if OpenMP is enabled)
NLEVELS=1     # Number of levels in multilevel algorithm

# 1=yes, 0=no
ENABLE_THETA=0         # Add theta term to the action
ENABLE_MULTICAN=0      # Add multicanonical potential to the action
ENABLE_MEAS_REPLICA=0  # Perform measurements on all PTBC replica
ENABLE_OPENMP=1        # Multithreading with OpenMP
ENABLE_DEBUG=0         # Perform extra sanity checks and debugging tests

# Available configurations: gcc, leonardo, galileo100, marconi
COMPILE_CONFIG=leonardo
COMPILE_TARGETS=(all)

BUILD_JOBS=112  # Number of parallel compilation jobs
BUILD_CLEAN=1   # Completely clean build (1=yes, 0=no)

# ========================================================================
# Compiler configuration
# ========================================================================

case "${COMPILE_CONFIG}" in

    gcc)
        CC=gcc
        CFLAGS='-O3'
        ;;

    leonardo)
        module load intel-oneapi-compilers/2023.2.1
        CC=icc
        CFLAGS='-O3 -axCORE-AVX512 -mtune=icelake -ip -ipo'
        ;;

    galileo100)
        module load intel
        CC=icc
        CFLAGS='-O3 -axCORE-AVX512 -mtune=cascadelake -ip -ipo'
        ;;

    marconi)
        module load intel
        CC=icc
        CFLAGS='-O3 -axCORE-AVX512 -mtune=skylake -ip -ipo'
        ;;

    *)
        echo "Error: unknown compiler configuration '${COMPILE_CONFIG}'." >&2
        echo "Available configurations: gcc, leonardo, galileo100, marconi" >&2
        exit 1
        ;;

esac

# Libraries
LIBS="-ldl -lz -lc"


# ========================================================================
# Display configuration
# ========================================================================

echo "========================================"
echo "Build configuration"
echo "========================================"
echo "N_c:            ${NCOLORS}"
echo "ST_dim:         ${STDIM}"
echo "Num_threads:    ${NTHREADS}"
echo "Num_levels:     ${NLEVELS}"
echo "OpenMP:         ${ENABLE_OPENMP}"
echo "Theta:          ${ENABLE_THETA}"
echo "Multicanonical: ${ENABLE_MULTICAN}"
echo "Replica meas:   ${ENABLE_MEAS_REPLICA}"
echo "Debug:          ${ENABLE_DEBUG}"
echo "========================================"
echo ""
echo "========================================"
echo "Compiler configuration"
echo "========================================"
echo "Compiler:       ${CC}"
echo "CFLAGS:         ${CFLAGS}"
echo "LIBS:           ${LIBS}"
echo "========================================"
echo ""


# ========================================================================
# Generate configure script
# ========================================================================

if (( BUILD_CLEAN )); then
    make distclean 2>/dev/null || true
fi
autoreconf -fi

# ========================================================================
# Configure
# ========================================================================

CONFIGURE_ARGS=(
    "N_c=${NCOLORS}"
    "ST_dim=${STDIM}"
    "Num_levels=${NLEVELS}"
    "Num_threads=${NTHREADS}"
    "CC=${CC}"
    "CFLAGS=${CFLAGS}"
    "LIBS=${LIBS}"
)

(( ENABLE_THETA ))        && CONFIGURE_ARGS+=(--enable-use-theta)
(( ENABLE_OPENMP ))       && CONFIGURE_ARGS+=(--enable-use-openmp)
(( ENABLE_MULTICAN ))     && CONFIGURE_ARGS+=(--enable-use-multicanonical)
(( ENABLE_MEAS_REPLICA )) && CONFIGURE_ARGS+=(--enable-meas-replica)
(( ENABLE_DEBUG ))        && CONFIGURE_ARGS+=(--enable-debug)

./configure "${CONFIGURE_ARGS[@]}"


# ========================================================================
# Build
# ========================================================================

make "${COMPILE_TARGETS[@]}" -j "${BUILD_JOBS}"
