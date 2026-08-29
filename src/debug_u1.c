#ifndef DEBUG_U1_C
#define DEBUG_U1_C

#include "../include/macro.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/u1.h"
#include "../include/u1_upd.h"

int main(void)
	{
	unsigned int const seed = 0;
	GParam param;

	U1 M, N, L, T, id_U1;

	// initialize random seed
	initrand(seed);

	// fix a value for d_beta
	param.d_beta = 0.9;

	printf("\n**********************************\n");
	printf("PROGRAM FOR THE DEBUGGING OF U(1)\n");
	printf("**********************************\n");

	printf("\n");
	printf("VERIFY THAT THE RANDOM MATRIX IS IN U(1)\n\n");

	printf("    Random matrix: ");
	rand_matrix_U1(&M);
	equal_U1(&N, &M);
	times_dag2_U1(&T, &M, &N);
	one_U1(&id_U1);
	minus_equal_U1(&T, &id_U1);

	REQUIRE(norm_U1(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT UPDATE U(1)->U(1)\n\n");

	// heatbath
	printf("    Heatbath: ");
	rand_matrix_U1(&M);
	rand_matrix_U1(&N);
	rand_matrix_U1(&L);
	plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)
	single_heatbath_U1(&M, &N, &param);
	equal_U1(&N, &M);
	times_dag2_U1(&T, &M, &N);
	minus_equal_U1(&T, &id_U1);
	REQUIRE(norm_U1(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");

	// overrelaxation
	printf("    Overrelaxation: ");
	rand_matrix_U1(&M);
	rand_matrix_U1(&N);
	rand_matrix_U1(&L);
	plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)
	single_overrelaxation_U1(&M, &N);
	equal_U1(&N, &M);
	times_dag2_U1(&T, &M, &N);
	minus_equal_U1(&T, &id_U1);
	REQUIRE(norm_U1(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
	rand_matrix_U1(&M);
	rand_matrix_U1(&N);
	rand_matrix_U1(&L);
	plus_equal_U1(&N, &L);       // N+=L,  M in U(1), N no   (M=link, N=staple)
	times_U1(&T, &M, &N);        // T=M*N
	double energy = retr_U1(&T); // initial energy
	single_overrelaxation_U1(&M, &N);
	times_U1(&T, &M, &N);  // T=M*N
	energy -= retr_U1(&T); // -=final energy
	REQUIRE(fabs(energy) < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
	rand_matrix_U1(&M);
	rand_matrix_U1(&N);
	rand_matrix_U1(&L);
	plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)

	times_U1(&T, &M, &N); // T=M*N
	energy = retr_U1(&T); // initial energy

	cool_U1(&M, &N);

	times_U1(&T, &M, &N);  // T=M*N
	energy -= retr_U1(&T); // -=final energy
	REQUIRE(energy < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED\n");

	printf("\nALL TESTS PASSED\n\n");

	return EXIT_SUCCESS;
	}


#endif
