#ifndef DEBUG_RNG_C
#define DEBUG_RNG_C

#include "../include/macro.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/random.h"


// print a few random numbers of the host generator for the given seed
static void print_sequence(unsigned int seed)
	{
	initrand(seed);
	if(seed != 0) printf("seed = %u\n", seed);
	else printf("seed=time()\n");

	for(int i = 0; i < 5; i++)
		{
		printf("    random[%d] = %.16g\n", i, casuale(rng_host()));
		}
	printf("\n");
	}


// print table header for uniform distribution tests
static void print_uniform_header(void)
	{
	printf("\n");
	printf("%-28s %12s %12s %8s %12s %12s %8s %9s\n",
	       "TEST", "MEAN", "EXPECTED", "SIGMA",
	       "VARIANCE", "EXPECTED", "SIGMA", "STATUS");

	printf("%-28s %12s %12s %8s %12s %12s %8s %9s\n",
	       "----------------------------", "------------", "------------", "--------", "------------", "------------", "--------", "---------");
	}


// mean and variance of n uniform numbers in (0,1):
// expected mean = 1/2, variance = 1/12
static int check_uniform(rng_t *rng, long n, char const *const name)
	{
	double sum = 0.0, sum2 = 0.0, min = 1.0, max = 0.0;

	for(long i = 0; i < n; i++)
		{
		double const x = casuale(rng);
		sum += x;
		sum2 += x * x;
		if(x < min) min = x;
		if(x > max) max = x;
		}

	double const mean = sum / (double) n;
	double const var = sum2 / (double) n - mean * mean;

	// 5 sigma tolerances:
	// sigma(mean) = sqrt(1/12/n)
	// sigma(var)  = sqrt(1/180/n)
	double const dev_mean = fabs(mean - 0.5) / sqrt(1.0 / 12.0 / (double) n);
	double const dev_var = fabs(var - 1.0 / 12.0) / sqrt(1.0 / 180.0 / (double) n);

	int const ok = (dev_mean < 5.0 && dev_var < 5.0 && min > 0.0 && max < 1.0);
	printf("%-28s %12.6f %12.6f %8.2f %12.6f %12.6f %8.2f %9s\n",
	       name, mean, 0.5, dev_mean, var, 1.0 / 12.0, dev_var, ok ? "OK" : "FAILED");
	return ok;
	}


#ifdef RNG_PHILOX_MODE

// known-answer tests of Philox4x32-10 (KAT vectors of Random123)
static int check_philox_kat(void)
	{
	uint32_t const ctr[3][4] = {
			{0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u},
			{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu},
			{0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u}
			};

	uint32_t const key[3][2] = {
			{0x00000000u, 0x00000000u},
			{0xffffffffu, 0xffffffffu},
			{0xa4093822u, 0x299f31d0u}
			};

	uint32_t const expected[3][4] = {
			{0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u},
			{0x408f276du, 0x41c83b0eu, 0xa20bc7c6u, 0x6d5451fdu},
			{0xd16cfe09u, 0x94fdccebu, 0x5001e420u, 0x24126ea1u}
			};

	int ok = 1;
	printf("\n");
	printf("%-6s   %-35s   %-35s   %9s\n", "TEST", "ACTUAL", "EXPECTED", "STATUS");
	printf("%-6s   %-35s   %-35s   %9s\n", "------", "-----------------------------------", "-----------------------------------", "---------");

	for(int t = 0; t < 3; t++)
		{
		uint32_t out[4];
		philox4x32_10(ctr[t], key[t], out);
		int const this_ok = (out[0] == expected[t][0] && out[1] == expected[t][1] && out[2] == expected[t][2] && out[3] == expected[t][3]);
		printf("%-6d   %08x %08x %08x %08x   %08x %08x %08x %08x   %9s\n",
		       t,
		       out[0], out[1], out[2], out[3],
		       expected[t][0], expected[t][1],
		       expected[t][2], expected[t][3],
		       this_ok ? "OK" : "FAILED");
		ok = ok && this_ok;
		}

	return ok;
	}


// correlation between streams of two neighbouring sites:
// expected correlation = 0 within 5 sigma
static int check_site_correlation(long n)
	{
	rng_t a, b;

	rng_init_stream(&a, g_rng_seed, 0, 0);
	rng_init_stream(&b, g_rng_seed, 0, 1);
	double sum = 0.0;
	for(long i = 0; i < n; i++)
		{
		sum += (casuale(&a) - 0.5) * (casuale(&b) - 0.5);
		}
	double const corr = 12.0 * sum / (double) n;
	double const dev_corr = fabs(corr) * sqrt((double) n);
	int const ok = (dev_corr < 5.0);
	printf("%-28s %12.6e %12.6e %8.2f %9s\n", "sites (0,0)-(0,1)", corr, 0.0, dev_corr, ok ? "OK" : "FAILED");
	return ok;
	}

#endif


int main(void)
	{
	print_sequence(1);
	print_sequence(2);
	print_sequence(1);
	print_sequence(0);

	int ok = 1;
	long const n = 10000000;

	#ifdef RNG_PHILOX_MODE
	printf("\nRNG: Philox4x32-10\n");
	#else
	printf("\nRNG: dSFMT\n");
	#endif

	initrand(1);

	print_uniform_header();
	ok = check_uniform(rng_host(), n, "host generator") && ok;

	#ifdef RNG_PHILOX_MODE

	printf("\nPhilox4x32-10 known-answer tests\n");
	ok = check_philox_kat() && ok;

	rng_t site;
	rng_init_stream(&site, g_rng_seed, 3, 12345);

	print_uniform_header();
	ok = check_uniform(&site, n, "site generator (3, 12345)") && ok;

	printf("\nSite stream correlation\n");
	printf("%-28s %12s %12s %8s %9s\n", "TEST", "MEASURED", "EXPECTED", "SIGMA", "STATUS");
	printf("%-28s %12s %12s %8s %9s\n", "----------------------------", "------------", "------------", "--------", "---------");
	ok = check_site_correlation(n) && ok;

	#endif

	printf("\n%s\n", ok ? "ALL TESTS PASSED" : "SOME TESTS FAILED");

	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
	}

#endif
