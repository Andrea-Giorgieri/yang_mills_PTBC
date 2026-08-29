#ifndef DEBUG_RNG_C
#define DEBUG_RNG_C

#include <stdio.h>
#include <stdlib.h>

#include "../include/random.h"

int main(void)
	{
	unsigned int seed = 1;
	initrand(seed);
	printf("seed = %u\n", seed);
	for(int i = 0; i < 5; i++)
		{
		printf("    random[%d] = %.16g\n", i, casuale());
		}

	printf("\n");

	seed = 2;
	initrand(seed);
	printf("seed = %u\n", seed);
	for(int i = 0; i < 5; i++)
		{
		printf("    random[%d] = %.16g\n", i, casuale());
		}

	printf("\n");

	seed = 1;
	initrand(seed);
	printf("seed = %u\n", seed);
	for(int i = 0; i < 5; i++)
		{
		printf("    random[%d] = %.16g\n", i, casuale());
		}

	printf("\n");

	seed = 0;
	initrand(seed);
	printf("seed=time()\n");
	for(int i = 0; i < 5; i++)
		{
		printf("    random[%d] = %.16g\n", i, casuale());
		}

	printf("\n");

	return EXIT_SUCCESS;
	}

#endif
