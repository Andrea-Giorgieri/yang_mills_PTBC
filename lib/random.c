#ifndef RANDOM_C
#define RANDOM_C

#include "../include/macro.h"

#include <limits.h>
#include <time.h>

#include "../include/random.h"


#ifdef RNG_PHILOX_MODE
rng_t g_rng_host;        // generator for the serial code on the host
uint32_t g_rng_seed = 1; // seed set by initrand
#else
dsfmt_t rng_status[NTHREADS]; // one generator per thread
#endif


// initialize random generator(s)
void initrand(unsigned int const s)
	{
	unsigned int seed = s;

	if(s == 0)
		{
		seed = ((unsigned int) time(NULL) + 10) % UINT_MAX;
		if(seed == 0)
			{
			seed = 1;
			}
		}

	#ifdef RNG_PHILOX_MODE
	// the streams of the lattice sites are initialized together with the gauge configurations (init_rng_gauge_conf)
	g_rng_seed = (uint32_t) seed;
	rng_init_stream(&g_rng_host, g_rng_seed, RNG_HOST_STREAM, 0);
	#else
	for(int thread_id = 0; thread_id < NTHREADS; thread_id++)
		{
		dsfmt_init_gen_rand(&(rng_status[thread_id]), seed);
		seed = (seed + 10) % UINT_MAX;
		if(seed == 0)
			{
			seed = 1;
			}
		}
	#endif
	}


#endif
