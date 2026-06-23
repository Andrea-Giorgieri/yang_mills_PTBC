#ifndef MACRO_H
#define MACRO_H

#include"../config.h"

#if NCOLOR == 1
    #define GAUGE_GROUP U1
#elif NCOLOR == 2
    #define GAUGE_GROUP Su2
#else
    #define GAUGE_GROUP SuN
#endif

#define MIN_VALUE 1.0e-13
#define INT_ALIGN 16
#define DOUBLE_ALIGN 32
#define STD_STRING_LENGTH 150

static const double PI = 3.141592653589793238462643383279502884197169399375105820974944;
static const double PI2 = 6.283185307179586476925286766559005768394338798750211641949889;
static const double HALF_PI = 1.570796326794896619231321691639751442098584699687552910487472;
static const double PI2_N = PI2 / (double)NCOLOR;
static const int MAX_POLY_PWR = (int)(NCOLOR / 2 + 1);


// function to access matrix elements
#define m(X,Y) ((X) * NCOLOR + (Y))
#define madj(X,Y) ((X) * (NCOLOR * NCOLOR - 1) + (Y))

// function to validate conditions
#define REQUIRE(cond, ...)                                    \
    do {                                                      \
        if (!(cond)) {                                        \
            fprintf(stderr,                                   \
                    "ERROR (%s:%d:%s): ",                     \
                    __FILE__, __LINE__, __func__);            \
            fprintf(stderr, __VA_ARGS__);                     \
            fprintf(stderr, "\n");                            \
            exit(EXIT_FAILURE);                               \
        }                                                     \
    } while (0)

// function for debugging checks
#ifdef DEBUG
#define ASSERT(cond, ...)                                     \
    do {                                                      \
        if (!(cond)) {                                        \
            fprintf(stderr,                                   \
                    "ASSERTION FAILED (%s:%d:%s): ",          \
                    __FILE__, __LINE__, __func__);            \
            fprintf(stderr, __VA_ARGS__);                     \
            fprintf(stderr, "\n");                            \
            abort();                                          \
        }                                                     \
    } while (0)
#else
#define ASSERT(cond, fmt, ...) ((void)0)
#endif

// function for probing numeric variables (cast as double)
#define PROBE(x)                                              \
    fprintf(stderr, "[PROBE] %s:%d:%s: %s = %g\n",            \
            __FILE__, __LINE__, __func__, #x, (double)(x))


// function to print flags
#define HERE(msg)                                             \
    do {                                                      \
        fprintf(stderr, "[HERE] %s:%d:%s: %s\n",              \
                __FILE__, __LINE__, __func__, (msg));         \
    } while(0)


// way to print a macro: if
// #define val1 val2
// then QUOTEME(val1) give the string "val2"
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

// to activate posix_memalign in stdlib.h
#define _POSIX_C_SOURCE 200809L

#endif
