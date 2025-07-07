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


void get_conf_files(char const * const arg1, char const * const arg2, char const * const arg3,
					char * conf_file_in, char * conf_file_out, char * twist_file_out)
	{
	if(strlen(arg1) >= STD_STRING_LENGTH)
		{
		fprintf(stderr, "conf_file_in name too long. Increase STD_STRING_LENGTH in include/macro.h\n");
		exit(EXIT_FAILURE);
		}
	strcpy(conf_file_in, arg1);
	if(strlen(arg2) >= STD_STRING_LENGTH)
		{
		fprintf(stderr, "conf_file_out name too long. Increase STD_STRING_LENGTH in include/macro.h\n");
		exit(EXIT_FAILURE);
		}
	strcpy(conf_file_out, arg2);
	if(strlen(arg3) >= STD_STRING_LENGTH)
		{
		fprintf(stderr, "twist_file_out name too long. Increase STD_STRING_LENGTH in include/macro.h\n");
		exit(EXIT_FAILURE);
		}
	strcpy(twist_file_out, arg3);
	}


void get_twist_params(char const * const arg1, char const * const arg2, char const * const arg3,
						int * twist_mu, int * twist_nu, int * twist_k)
	{
	if(sscanf(arg1, "%d", twist_mu) != 1)
		{
		fprintf(stderr, "Error reading the first twist plane direction (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(sscanf(arg2, "%d", twist_nu) != 1)
		{
		fprintf(stderr, "Error reading the second twist plane direction (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(sscanf(arg3, "%d", twist_k) != 1)
		{
		fprintf(stderr, "Error reading the twist factor (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}


void get_conf_sizes(char ** argv, int argc, int * sizes)
	{
	if(argc != STDIM)
		{
		fprintf(stderr, "Expected %d sizes, got %d (%s, %d)\n", STDIM, argc, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	for(int i = 0; i < STDIM; i++)
		{
		if(sscanf(argv[i], "%d", &sizes[i]) != 1)
			{
			fprintf(stderr, "Error reading the conf size %d (%s, %d)\n", i, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		}
	}


void set_param_for_writing_conf(int const * const sizes, int twist_mu, int twist_nu, int twist_k, char const * const conf_file_out,
								char const * const twist_file_out, GParam * const param)
	{
	int i, j;
	
	set_defaults(param);
	
	for(i=0; i<STDIM; i++) param->d_size[i] = sizes[i];
	
	param->d_k_twist[dirs_to_si(twist_mu, twist_nu)] = twist_k;
	strcpy(param->d_conf_file, conf_file_out);
	strcpy(param->d_twist_file, twist_file_out);
	
	param->d_volume=1;
	for(i=0; i<STDIM; i++) param->d_space_vol[i] = 1;
	for(i=0; i<STDIM; i++)
		{
		(param->d_volume)*=(param->d_size[i]);
		for(j=0; j<STDIM; j++) if (j != i) (param->d_space_vol[j]) *= (param->d_size[i]);
		}
	param->d_n_planes = STDIM * (STDIM-1);
	}


void parse_line(char * line, int *cartcoord, int *mu, int *i, int *j, double *re, double *im)
	{
	int k, step;
	char *line_ptr;
	
	line_ptr = line;
	if(sscanf(line_ptr, "%*d %*d%n", &step) != 0)
		{
		fprintf(stderr, "Error reading indeces to be ignored (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	line_ptr += step * (int)sizeof(char);

	for(k = 0; k < STDIM; k++)
		{
		if(sscanf(line_ptr, "%d%n", &cartcoord[k], &step) != 1)
			{
			fprintf(stderr, "Error reading a cartesian coordinate along direction %d (%s, %d)\n", k, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		line_ptr += step * (int)sizeof(char);
		}

	if(sscanf(line_ptr, "%d%n", mu, &step) != 1)
		{
		fprintf(stderr, "Error reading a link direction (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	line_ptr += step * (int)sizeof(char);

	if(sscanf(line_ptr, "%d %d%n", i, j, &step) != 2)
		{
		fprintf(stderr, "Error reading a color index (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	line_ptr += step * (int)sizeof(char);

	if(sscanf(line_ptr, " (%lf,%lf)", re, im) != 2)
		{
		fprintf(stderr, "Error reading a link component (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}


void check_valid_line(int const * const cartcoord, int mu, int i, int j, GParam const * const param)
	{
	for(int k = 0; k < STDIM; k++)
		if(cartcoord[k] < 0 || cartcoord[k] >= param->d_size[k])
			{
			fprintf(stderr, "Invalid cartesian coordinate %d along direction %d (%s, %d)\n", cartcoord[k], k, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	if(mu < 0 || mu >= STDIM)
		{
		fprintf(stderr, "Invalid spacetime direction %d (%s, %d)\n", mu, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i < 1 || i > NCOLOR || j < 1 || j > NCOLOR)
		{
		fprintf(stderr, "Invalid color index (%d,%d) (%s, %d)\n", i, j, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}


void init_gauge_conf_from_txt(char const * const conf_file_in, Gauge_Conf * const GC, GParam const * const param)
	{
	char line[STD_STRING_LENGTH];
	int mu, i, j, cartcoord[STDIM];
	long r, num_lines;
	double re, im;
	
	FILE *fp = fopen(conf_file_in, "r");
	if(fp == NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_in, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice), param->d_volume, __FILE__, __LINE__);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r)
	#endif
	for(r = 0; r < (param->d_volume); r++) 
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice[r]), STDIM, __FILE__, __LINE__);
		}
	
	num_lines = 0;
	while(fscanf(fp, "%[^\n]\n", line) != EOF)
		{
		parse_line(line, cartcoord, &mu, &i, &j, &re, &im);
		check_valid_line(cartcoord, mu, i, j, param);
		i -= 1;
		j -= 1;
		GC->lattice[cart_to_si(cartcoord, param)][mu].comp[m(i,j)] = re + I * im;
		num_lines++;
		}
	fclose(fp);
	
	r = STDIM * NCOLOR * NCOLOR * param->d_volume;
	if(num_lines != r)
		{
		fprintf(stderr, "Expected %ld lines, got %ld (%s, %d)\n", num_lines, r, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	GC->conf_label = 0;
	GC->replica_index = 0;
	GC->update_index = 0;
	}


void init_twist_cond_from_txt(Gauge_Conf * const GC, GParam const * const param)
	{
	long r;
	int i, j, x_mu, x_nu;
	int cartcoord[STDIM];

	allocate_array_double_complex_pointer(&(GC->Z), param->d_volume, __FILE__, __LINE__);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r, j)
	#endif 
	for(r=0; r<(param->d_volume); r++)
		{
		allocate_array_double_complex(&(GC->Z[r]), param->d_n_planes, __FILE__, __LINE__);
		for(j=0; j<param->d_n_planes; j++) GC->Z[r][j] = 1.0 + 0.0 * I;
		}

	x_mu = 0;
	x_nu = 0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r, i, j)
	#endif 
	for(r=0; r<param->d_volume; r++)
		{
		si_to_cart(cartcoord, r, param);
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
				if (cartcoord[i] == x_mu && cartcoord[j] == x_nu)
					{
					GC->Z[r][dirs_to_si(i,j)] = cexp(I*PI2_N*(param->d_k_twist[dirs_to_si(i,j)]));
					GC->Z[r][dirs_to_si(j,i)] = conj(GC->Z[r][dirs_to_si(i,j)]);
					}
		}
	}


void free_gauge_conf_from_txt(Gauge_Conf *GC, GParam const * const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif 
	for(long r=0; r<param->d_volume; r++)
		{
		free(GC->lattice[r]);
		free(GC->Z[r]);
		}
	free(GC->lattice);
	free(GC->Z);
	}


int main (int argc, char **argv)
	{
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
	get_conf_sizes(&argv[7], argc-7, sizes);

	set_param_for_writing_conf(sizes, twist_mu, twist_nu, twist_k, conf_file_out, twist_file_out, &param);

	init_indexing_lexeo();
	init_gauge_conf_from_txt(conf_file_in, &GC, &param);
	init_twist_cond_from_txt(&GC, &param);
	fprintf(stdout, "Done reading input conf\n");

	write_conf_on_file(&GC, &param);
	fprintf(stdout, "Done writing output conf\n");
	
	//print_on_file(stderr, &(GC.lattice[511][0]));
	
	free_gauge_conf_from_txt(&GC, &param);

	return EXIT_SUCCESS;
	}

#endif