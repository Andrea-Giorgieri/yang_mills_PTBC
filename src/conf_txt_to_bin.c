#ifndef CONF_TXT_TO_BIN_C
#define CONF_TXT_TO_BIN_C

#include"../include/macro.h"
#ifdef HASH_MODE
#include<openssl/md5.h>
#endif

#ifdef OPENMP_MODE
#include<omp.h>
#endif

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/memalign.h"
#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"


void get_conf_files(char const *const arg1, char const *const arg2, char const *const arg3,
                    char *conf_file_in, char *conf_file_out, char *twist_file_out)
	{
	REQUIRE(strlen(arg1) < STD_STRING_LENGTH, "conf_file_in name too long. Increase STD_STRING_LENGTH in include/macro.h");
	REQUIRE(strlen(arg2) < STD_STRING_LENGTH, "conf_file_out name too long. Increase STD_STRING_LENGTH in include/macro.h");
	REQUIRE(strlen(arg3) < STD_STRING_LENGTH, "twist_file_out name too long. Increase STD_STRING_LENGTH in include/macro.h");

	strcpy(conf_file_in, arg1);
	strcpy(conf_file_out, arg2);
	strcpy(twist_file_out, arg3);
	}


void get_twist_params(char const *const arg1, char const *const arg2, char const *const arg3,
                      int *twist_mu, int *twist_nu, int *twist_k)
	{
	REQUIRE(arg1 != NULL, "first twist plane direction is NULL");
	REQUIRE(arg2 != NULL, "second twist plane direction is NULL");
	REQUIRE(arg3 != NULL, "twist factor is NULL");

	REQUIRE(sscanf(arg1, "%d", twist_mu) == 1, "failed to read the first twist plane direction (arg=%s)", arg1);
	REQUIRE(sscanf(arg2, "%d", twist_nu) == 1, "failed to read the second twist plane direction (arg=%s)", arg2);
	REQUIRE(sscanf(arg3, "%d", twist_k) == 1, "failed to read the twist factor (arg=%s)", arg3);
	}


void get_conf_sizes(char **argv, int argc, int *sizes)
	{
	REQUIRE(argc == STDIM, "expected %d sizes, got %d", STDIM, argc);
	for(int i = 0; i < STDIM; i++)
		{
		REQUIRE(argv[i] != NULL, "NULL argument at position %d", i);
		REQUIRE(sscanf(argv[i], "%d", &sizes[i]) == 1, "failed to read the conf size %d (arg=%s)", i, argv[i]);
		}
	}


void set_param_for_writing_conf(int const *const sizes, int twist_mu, int twist_nu, int twist_k, char const *const conf_file_out,
                                char const *const twist_file_out, GParam *const param)
	{
	set_defaults(param);

	for(int i = 0; i < STDIM; i++) param->d_size[i] = sizes[i];

	param->d_k_twist[dirs_to_si(twist_mu, twist_nu)] = twist_k;
	strcpy(param->d_conf_file, conf_file_out);
	strcpy(param->d_twist_file, twist_file_out);

	param->d_volume = 1;
	for(int i = 0; i < STDIM; i++) param->d_space_vol[i] = 1;
	for(int i = 0; i < STDIM; i++)
		{
		(param->d_volume) *= (param->d_size[i]);
		for(int j = 0; j < STDIM; j++) if(j != i) (param->d_space_vol[j]) *= (param->d_size[i]);
		}
	param->d_n_planes = STDIM * (STDIM - 1);
	}


void parse_line(char *line, int *cartcoord, int *mu, int *i, int *j, double *re, double *im)
	{
	int step;
	char const *line_ptr = line;

	REQUIRE(sscanf(line_ptr, "%*d %*d%n", &step) == 0, "failed to read indices to be ignored");
	line_ptr += step;

	for(int k = 0; k < STDIM; k++)
		{
		REQUIRE(sscanf(line_ptr, "%d%n", &cartcoord[k], &step) == 1, "failed to read cartesian coordinate %d", k);
		line_ptr += step;
		}

	REQUIRE(sscanf(line_ptr, "%d%n", mu, &step) == 1, "failed to read link direction");
	line_ptr += step;

	REQUIRE(sscanf(line_ptr, "%d %d%n", i, j, &step) == 2, "failed to read color indices");
	line_ptr += step;

	REQUIRE(sscanf(line_ptr, " (%lf,%lf)", re, im) == 2, "failed to read complex link component");
	}


void check_valid_line(int const *const cartcoord, int mu, int i, int j, GParam const *const param)
	{
	for(int k = 0; k < STDIM; k++)
		{
		REQUIRE(cartcoord[k] >= 0 && cartcoord[k] < param->d_size[k], "invalid cartesian coordinate %d along direction %d", cartcoord[k], k);
		}
	REQUIRE(mu >= 0 && mu < STDIM, "invalid spacetime direction %d", mu);
	REQUIRE(i >= 1 && i <= NCOLOR && j >= 1 && j <= NCOLOR, "invalid color indices (%d,%d)", i, j);
	}


void init_gauge_conf_from_txt(char const *const conf_file_in, Gauge_Conf *const GC, GParam const *const param)
	{
	char line[STD_STRING_LENGTH];
	int mu, i, j, cartcoord[STDIM];
	double re, im;

	FILE *fp = fopen(conf_file_in, "r");
	REQUIRE(fp != NULL, "failed to open file %s", conf_file_in);

	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice), param->d_volume, __FILE__, __LINE__);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice[r]), STDIM, __FILE__, __LINE__);
		}

	long num_lines = 0;

	#if NCOLOR < 3
	REQUIRE(0, "not implemented");
	#else
	while(fscanf(fp, "%[^\n]\n", line) != EOF)
		{
		parse_line(line, cartcoord, &mu, &i, &j, &re, &im);
		check_valid_line(cartcoord, mu, i, j, param);
		i -= 1;
		j -= 1;
		GC->lattice[cart_to_si(cartcoord, param)][mu].comp[m(i, j)] = re + I * im;
		num_lines++;
		}
	#endif

	fclose(fp);

	long const num_lines_expected = STDIM * NCOLOR * NCOLOR * param->d_volume;
	REQUIRE(num_lines == num_lines_expected, "expected %ld lines, got %ld", num_lines_expected, num_lines);

	GC->conf_label = 0;
	GC->replica_index = 0;
	GC->update_index = 0;
	}


void init_twist_cond_from_txt(Gauge_Conf *const GC, GParam const *const param)
	{
	allocate_array_double_complex_pointer(&(GC->Z), param->d_volume, __FILE__, __LINE__);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_double_complex(&(GC->Z[r]), param->d_n_planes + 1, __FILE__, __LINE__);
		for(int j = 0; j < param->d_n_planes + 1; j++) GC->Z[r][j] = 1.0 + 0.0 * I;
		}

	int const x_mu = 0;
	int const x_nu = 0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		int cartcoord[STDIM];
		si_to_cart(cartcoord, r, param);
		for(int i = 0; i < STDIM; i++)
			for(int j = i + 1; j < STDIM; j++)
				if(cartcoord[i] == x_mu && cartcoord[j] == x_nu)
					{
					GC->Z[r][dirs_to_si(i, j)] = cexp(I * PI2_N * (param->d_k_twist[dirs_to_si(i, j)]));
					GC->Z[r][dirs_to_si(j, i)] = conj(GC->Z[r][dirs_to_si(i, j)]);
					}
		}
	}


void free_gauge_conf_from_txt(Gauge_Conf *GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		free(GC->lattice[r]);
		free(GC->Z[r]);
		}
	free(GC->lattice);
	free(GC->Z);
	}


int main(int argc, char **argv)
	{
	#if NCOLOR < 3

	printf("This program can be used only for NCOLOR > 2\n");

	#else

	char conf_file_in[STD_STRING_LENGTH], conf_file_out[STD_STRING_LENGTH], twist_file_out[STD_STRING_LENGTH];
	int twist_mu, twist_nu, twist_k, sizes[STDIM];
	Gauge_Conf GC;
	GParam param;

	#ifdef OPENMP_MODE
	omp_set_max_active_levels(1);
	#endif

	if(argc < 8)
		{
		printf("Usage: %s conf_file_in conf_file_out twist_file_out twist_mu twist_nu twist_k size_0 size_1 ...\n\n", argv[0]);
		print_compilation_details();
		return EXIT_SUCCESS;
		}

	get_conf_files(argv[1], argv[2], argv[3], conf_file_in, conf_file_out, twist_file_out);
	get_twist_params(argv[4], argv[5], argv[6], &twist_mu, &twist_nu, &twist_k);
	get_conf_sizes(&argv[7], argc - 7, sizes);

	set_param_for_writing_conf(sizes, twist_mu, twist_nu, twist_k, conf_file_out, twist_file_out, &param);

	init_indexing_lexeo();
	init_gauge_conf_from_txt(conf_file_in, &GC, &param);
	init_twist_cond_from_txt(&GC, &param);
	fprintf(stdout, "Done reading input conf\n");

	write_conf_on_file(&GC, &param);
	fprintf(stdout, "Done writing output conf\n");

	//print_on_file(stderr, &(GC.lattice[511][0]));

	free_gauge_conf_from_txt(&GC, &param);

	#endif

	return EXIT_SUCCESS;
	}

#endif
