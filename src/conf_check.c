#ifndef CONF_CHECK_C
#define CONF_CHECK_C

#include "../include/macro.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HASH_MODE
#include <openssl/md5.h>
#endif

#include "../include/gauge_group.h"
#include "../include/memalign.h"

// get the spacetime dimension
void getspacetimedim(char *infile, int *dim)
	{
	FILE *fp = fopen(infile, "r");
	REQUIRE(fp != NULL, "failed to open input file %s", infile);

	int err = fscanf(fp, "%d", dim);
	REQUIRE(err == 1, "failed to read number of dimensions from input file %s", infile);

	fclose(fp);
	}


// read the size and the hash of the configuration
void getsizeandhash(char *infile, int *sides, char *hash)
	{
	long update_index;
	int dim;

	FILE *fp = fopen(infile, "r");
	REQUIRE(fp != NULL, "failed to open input file %s", infile);
	int err = fscanf(fp, "%d", &dim);
	REQUIRE(err == 1, "failed to read number of dimensions from input file %s", infile);

	for(int i = 0; i < dim; i++)
		{
		err = fscanf(fp, "%d", &(sides[i]));
		REQUIRE(err == 1, "failed to read %d-th size from input file %s", i, infile);
		}

	err = fscanf(fp, "%ld %s", &update_index, hash);
	REQUIRE(err == 2, "failed to read update index or hash from input file %s", infile);

	fclose(fp);
	}


// compute the hash
void computehash(char *infile, int dim, long volume, char *hash)
	{
	#ifdef HASH_MODE

	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];

	// open the configuration file in binary
	FILE *fp = fopen(infile, "rb");
	REQUIRE(fp != NULL, "failed to open input file %s in binary mode", infile);

	// read again the header:
	int tmp = 0;
	while(tmp != '\n')
		{
		tmp = fgetc(fp);
		}

	// read the configuration & compute md5sum
	MD5_Init(&mdContext);
	for(long r = 0; r < volume; r++)
		{
		for(int i = 0; i < dim; i++)
			{
			GAUGE_GROUP link;
			read_from_binary_file_bigen(fp, &link);
			#if NCOLOR == 1
			MD5_Update(&mdContext, &(link.comp), sizeof(double complex));
			#elif NCOLOR == 2
			for(int j = 0; j < 4; j++)
				{
				MD5_Update(&mdContext, &(link.comp[j]), sizeof(double));
				}
			#else
			for(int j = 0; j < NCOLOR * NCOLOR; j++)
				{
				MD5_Update(&mdContext, &(link.comp[j]), sizeof(double complex));
				}
			#endif
			}
		}
	MD5_Final(c, &mdContext);
	for(long r = 0; r < MD5_DIGEST_LENGTH; r++)
		{
		sprintf(&(hash[2 * r]), "%02x", c[r]);
		}

	fclose(fp);

	#else

	// just to avoid compile time warnings
	(void) infile;
	(void) dim;
	(void) volume;
	*hash = '0';

	#endif
	}


// MAIN
int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		printf("Usage: %s conf_file\n\n", argv[0]);
		print_compilation_details();

		return EXIT_SUCCESS;
		}

	char infile[STD_STRING_LENGTH];
	int dim, *sides;

	#ifdef HASH_MODE
	long volume;
	char md5sum_old[2 * MD5_DIGEST_LENGTH + 1];
	char md5sum_new[2 * MD5_DIGEST_LENGTH + 1];
	#else
	char md5sum_old[2 * STD_STRING_LENGTH + 1] = {0};
	#endif

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");

	strcpy(infile, argv[1]);

	// get spacetime dim
	getspacetimedim(infile, &dim);

	allocate_array_int(&sides, dim, __FILE__, __LINE__);

	// get lattice size and initial hash
	getsizeandhash(infile, sides, md5sum_old);

	#ifdef HASH_MODE
	// total volume
	volume = 1;
	for(int i = 0; i < dim; i++)
		{
		volume *= sides[i];
		}

	// compute the hash
	computehash(infile, dim, volume, md5sum_new);

	// check md5sum computed and stored
	int err = strncmp(md5sum_old, md5sum_new, 2 * MD5_DIGEST_LENGTH + 1);
	REQUIRE(err == 0, "the configuration %s is corrupted!", infile);
	fprintf(stdout, "Check passed\n");

	#endif

	free(sides);

	return EXIT_SUCCESS;
	}

#endif
