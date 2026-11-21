#ifndef DEBUG_SUN_C
#define DEBUG_SUN_C

#include "../include/macro.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/sun.h"
#include "../include/sun_upd.h"


int main(void)
	{
	unsigned int const seed = 0;
	GParam param;

	SuN M, N, L, T;

	// initialize random seed
	initrand(seed);

	// fix a value for d_beta
	param.d_beta = 2.3;

	printf("\n***********************************\n");
	printf("PROGRAM FOR THE DEBUGGING OF SU(N)\n");
	printf("***********************************\n\n");

	printf("N = %s", QUOTEME(NCOLOR));

	printf("\n\n");
	printf("VERIFY THAT THE RANDOM MATRIX IS IN SU(N)\n\n");
	printf("    Random matrix: ");
	rand_matrix_SuN(&M, rng_host());
	REQUIRE(scheck_SuN(&M) == 0, "TEST FAILED");
	printf("TEST PASSED\n");

	printf("\n\n");
	printf("VERIFY THAT UPDATE SU(N)->SU(N)\n\n");

	// heatbath
	printf("    Heatbath: ");
	rand_matrix_SuN(&M, rng_host());
	rand_matrix_SuN(&N, rng_host());
	rand_matrix_SuN(&L, rng_host());
	plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)
	single_heatbath_SuN(&M, &N, &param, rng_host());
	REQUIRE(scheck_SuN(&M) == 0, "TEST FAILED");
	printf("TEST PASSED\n");

	// overrelaxation
	printf("    Overrelaxation: ");
	rand_matrix_SuN(&M, rng_host());
	rand_matrix_SuN(&N, rng_host());
	rand_matrix_SuN(&L, rng_host());
	plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple) */
	single_overrelaxation_SuN(&M, &N, rng_host());
	REQUIRE(scheck_SuN(&M) == 0, "TEST FAILED");
	printf("TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
	rand_matrix_SuN(&M, rng_host());
	rand_matrix_SuN(&N, rng_host());
	rand_matrix_SuN(&L, rng_host());
	plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

	times_SuN(&T, &M, &N);        // T=M*N
	double energy = retr_SuN(&T); // initial energy
	single_overrelaxation_SuN(&M, &N, rng_host());
	times_SuN(&T, &M, &N);  // T=M*N
	energy -= retr_SuN(&T); // -=final energy
	REQUIRE(fabs(energy) < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
	rand_matrix_SuN(&M, rng_host());
	rand_matrix_SuN(&N, rng_host());
	rand_matrix_SuN(&L, rng_host());
	plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

	times_SuN(&T, &M, &N); // T=M*N
	energy = retr_SuN(&T); // initial energy

	cool_SuN(&M, &N);

	times_SuN(&T, &M, &N);  // T=M*N
	energy -= retr_SuN(&T); // -=final energy
	REQUIRE(energy < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED ;)\n");

	printf("\nALL TESTS PASSED\n\n");

	return EXIT_SUCCESS;
	}

#endif
