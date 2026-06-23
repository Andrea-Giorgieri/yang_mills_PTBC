#ifndef DEBUG_SU2_C
#define DEBUG_SU2_C

#include"../include/macro.h"
#include"../include/gparam.h"
#include"../include/random.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main(void)
	{
	unsigned int seed = 0;
	double energy;
	GParam param;

	Su2 M, N, L, T, id_SU2;

	// initialize random seed
	initrand(seed);

	// fix a value for d_beta
	param.d_beta = 2.3;

	printf("\n***********************************\n");
	printf("PROGRAM FOR THE DEBUGGING OF SU(2)\n");
	printf("***********************************\n");

	printf("\n");
	printf("VERIFY THAT THE RANDOM MATRIX IS IN SU(2)\n\n");

	printf("    Random matrix: ");
	rand_matrix_Su2(&M);
	equal_Su2(&N, &M);
	one_Su2(&id_SU2);
	times_dag2_Su2(&T, &M, &N);
	minus_equal_Su2(&T, &id_SU2);
	REQUIRE(norm_Su2(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT UPDATE SU(2)->SU(2)\n\n");

	// heatbath
	printf("    Heatbath: ");
	rand_matrix_Su2(&M);
	rand_matrix_Su2(&N);
	rand_matrix_Su2(&L);
	plus_equal_Su2(&N, &L); // N+=L,  M in SU(2), N no   (M=link, N=staple)
	single_heatbath_Su2(&M, &N, &param);
	equal_Su2(&N, &M);
	times_dag2_Su2(&T, &M, &N);
	minus_equal_Su2(&T, &id_SU2);
	REQUIRE(norm_Su2(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");

	// overrelaxation
	printf("    Overrelaxation: ");
	rand_matrix_Su2(&M);
	rand_matrix_Su2(&N);
	rand_matrix_Su2(&L);
	plus_equal_Su2(&N, &L); // N+=L,  M in SU(2), N no   (M=link, N=staple)
	single_overrelaxation_Su2(&M, &N);
	equal_Su2(&N, &M);
	times_dag2_Su2(&T, &M, &N);
	minus_equal_Su2(&T, &id_SU2);
	REQUIRE(norm_Su2(&T) <= MIN_VALUE, "TEST FAILED");
	printf("TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
	rand_matrix_Su2(&M);
	rand_matrix_Su2(&N);
	rand_matrix_Su2(&L);
	plus_equal_Su2(&N, &L); // N+=L,  M in SU(2), N no   (M=link, N=staple)

	times_Su2(&T, &M, &N); // T=M*N
	energy = retr_Su2(&T); // initial energy
	single_overrelaxation_Su2(&M, &N);
	times_Su2(&T, &M, &N);  // T=M*N
	energy -= retr_Su2(&T); // -=final energy
	REQUIRE(fabs(energy) < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED\n");


	printf("\n\n");
	printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
	rand_matrix_Su2(&M);
	rand_matrix_Su2(&N);
	rand_matrix_Su2(&L);
	plus_equal_Su2(&N, &L); // N+=L,  M in SU(2), N no   (M=link, N=staple)

	times_Su2(&T, &M, &N); // T=M*N
	energy = retr_Su2(&T); // initial energy

	cool_Su2(&M, &N);

	times_Su2(&T, &M, &N);  // T=M*N
	energy -= retr_Su2(&T); // -=final energy
	REQUIRE(energy < MIN_VALUE, "    TEST FAILED: Delta E = %f", energy);
	printf("    TEST PASSED\n");

	printf("\nALL TESTS PASSED\n\n");

	return EXIT_SUCCESS;
	}

#endif
