#ifndef RANDOM_H
#define RANDOM_H

#include "macro.h"

#include <math.h>
#include <stdint.h>

// Random number generators: the state of the generator (rng_t) is passed explicitly to every function
// that needs random numbers.
//
// default (dSFMT): one dSFMT generator per OpenMP thread, as in the original code
//   rng_thread() = generator of the calling thread; rng_host() = rng_thread()
//
// RNG_PHILOX_MODE (always with OpenACC, optional on CPU): counter-based Philox4x32-10 generator
//   (Salmon, Moraes, Dror, Shaw, "Parallel random numbers: as easy as 1, 2, 3", SC11)
//   - one stream per (replica, lattice site), stored in Gauge_Conf.rng_site and used by the link updates
//   - one stream for the serial code on the host, rng_host()
//   streams are distinct by construction: counter = (draw index [64 bit], site, stream), key = (seed, constant)

#ifdef RNG_PHILOX_MODE

typedef struct Rng_State
	{
	uint32_t ctr[4]; // ctr[0], ctr[1]: 64-bit draw counter; ctr[2]: site; ctr[3]: stream
	uint32_t key[2]; // key[0]: seed; key[1]: fixed constant
	} rng_t;

#define RNG_HOST_STREAM 0xFFFFFFFFu // stream index of the host generator
#define RNG_KEY_CONSTANT 0x6A09E667u // second word of the key (fractional part of sqrt(2))

extern rng_t g_rng_host;   // generator for the serial code on the host
extern uint32_t g_rng_seed; // seed set by initrand


// Philox4x32 with 10 rounds: out = bijection_key(ctr)
ACC_ROUTINE_SEQ
static inline void philox4x32_10(uint32_t const ctr[4], uint32_t const key[2], uint32_t out[4])
	{
	uint32_t c0 = ctr[0], c1 = ctr[1], c2 = ctr[2], c3 = ctr[3];
	uint32_t k0 = key[0], k1 = key[1];

	for(int round = 0; round < 10; round++)
		{
		uint64_t const p0 = (uint64_t) 0xD2511F53u * (uint64_t) c0;
		uint64_t const p1 = (uint64_t) 0xCD9E8D57u * (uint64_t) c2;

		c0 = (uint32_t) (p1 >> 32) ^ c1 ^ k0;
		c1 = (uint32_t) p1;
		c2 = (uint32_t) (p0 >> 32) ^ c3 ^ k1;
		c3 = (uint32_t) p0;

		k0 += 0x9E3779B9u; // key schedule (Weyl sequence)
		k1 += 0xBB67AE85u;
		}

	out[0] = c0;
	out[1] = c1;
	out[2] = c2;
	out[3] = c3;
	}


// initialize the stream (stream, site) with the given seed
ACC_ROUTINE_SEQ
static inline void rng_init_stream(rng_t *rng, uint32_t const seed, uint32_t const stream, uint32_t const site)
	{
	rng->ctr[0] = 0;
	rng->ctr[1] = 0;
	rng->ctr[2] = site;
	rng->ctr[3] = stream;
	rng->key[0] = seed;
	rng->key[1] = RNG_KEY_CONSTANT;
	}


// random number in (0,1): 52 random bits from one Philox block, then the counter is advanced
ACC_ROUTINE_SEQ
static inline double casuale(rng_t *rng)
	{
	uint32_t out[4];
	philox4x32_10(rng->ctr, rng->key, out);

	rng->ctr[0]++;
	if(rng->ctr[0] == 0) rng->ctr[1]++;

	uint64_t const x = ((uint64_t) out[0] << 32) | (uint64_t) out[1];
	return ((double) (x >> 12) + 0.5) * 2.220446049250313080847263336181640625e-16; // ((x >> 12) + 1/2) / 2^52
	}


// generator for the serial code on the host
static inline rng_t *rng_host(void)
	{
	return &g_rng_host;
	}

#else // dSFMT

#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "dSFMT.h"

typedef dsfmt_t rng_t;

extern dsfmt_t rng_status[NTHREADS]; // one generator per thread


// generator of the calling thread
static inline rng_t *rng_thread(void)
	{
	#ifdef OPENMP_MODE
	return &(rng_status[omp_get_thread_num()]);
	#else
	return &(rng_status[0]);
	#endif
	}


// random number in (0,1)
static inline double casuale(rng_t *rng)
	{
	return dsfmt_genrand_open_open(rng);
	}


// generator for the serial code on the host
static inline rng_t *rng_host(void)
	{
	return rng_thread();
	}

#endif


void initrand(unsigned int const s); // initialize random generator(s)


// normal gaussian random number generator (polar method, knuth vol 2, p. 117)
ACC_ROUTINE_SEQ
static inline double gauss1(rng_t *rng)
	{
	double v1, v2, s;

	do
		{
		v1 = 1.0 - 2.0 * casuale(rng);
		v2 = 1.0 - 2.0 * casuale(rng);
		s = v1 * v1 + v2 * v2;
		} while(s >= 1);

	return v1 * sqrt(-2 * log(s) / s);
	}


// normal gaussian random number generator (polar method, knuth vol 2, p. 117)
ACC_ROUTINE_SEQ
static inline void gauss2(double *res1, double *res2, rng_t *rng)
	{
	double v1, v2, s;

	do
		{
		v1 = -1.0 + 2.0 * casuale(rng);
		v2 = -1.0 + 2.0 * casuale(rng);
		s = v1 * v1 + v2 * v2;
		} while(s >= 1);

	*res1 = v1 * sqrt(-2 * log(s) / s);
	*res2 = v2 * sqrt(-2 * log(s) / s);
	}

#endif
