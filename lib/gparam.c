#ifndef GPARAM_C
#define GPARAM_C

#include "../include/macro.h"
#include "../include/endianness.h"
#include "../include/gparam.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// remove from input file white/empty lines and comments
// comments start with the charachter #
void remove_white_line_and_comments(FILE *input)
	{
	int temp_i;

	temp_i = getc(input);
	if(temp_i == '\n' || temp_i == ' ' || temp_i == '\043') // scan for white lines and comments
		{
		ungetc(temp_i, input);

		temp_i = getc(input);
		if(temp_i == '\n' || temp_i == ' ') // white line
			{
			do {
				temp_i = getc(input);
				} while(temp_i == '\n' || temp_i == ' ');
			}
		ungetc(temp_i, input);

		temp_i = getc(input);
		if(temp_i == '\043') // comment, \043 = ascii oct for #
			{
			do {
				temp_i = getc(input);
				} while(temp_i != '\n');
			}
		else
			{
			ungetc(temp_i, input);
			}

		remove_white_line_and_comments(input);
		}
	else
		{
		ungetc(temp_i, input);
		}
	}

void readinput(char *in_file, GParam *param)
	{
	FILE *input;
	char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
	double temp_d;
	int temp_i, i;
	int err, end = 1;
	unsigned int temp_ui;

	// this is to avoid unnecessary checks in case the multilevel is not used
	for(i = 0; i < NLEVELS; i++) { param->d_ml_step[i] = 0; }

	// just to avoid possible mistakes with uninitialized stuff
	for(i = 0; i < NCOLOR; i++) { param->d_h[i] = 0.0; }
	param->d_theta = 0.0;

	// to avoid possible mistakes with uninitialized stuff
	for(i = 0; i < STDIM - 1; i++) { param->d_L_defect[i] = 0; }
	param->d_defect_dir = 0;
	param->d_N_replica_pt = 1;

	// to avoid possible mistakes with uninitialized stuff
	param->d_ngfsteps = 0;
	param->d_gf_meas_each = 1;
	param->d_gfstep = 0.01;

	param->d_agf_length = 0.0;
	param->d_agf_meas_each = 1.0;
	param->d_agf_step = 0.01;
	param->d_agf_delta = 0.001;
	param->d_agf_time_bin = 0;

	param->d_coolsteps = 0;
	param->d_coolrepeat = 0;

	// default = compute only plaquette and topological charge
	param->d_plaquette_meas = 1;
	param->d_clover_energy_meas = 0;
	param->d_charge_meas = 1;
	param->d_chi_prime_meas = 0;
	param->d_polyakov_meas = 0;
	param->d_topcharge_tprof_meas = 0;

	input = fopen(in_file, "r"); // open the input file
	if(input == NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		while(end == 1) // slide the file
			{
			remove_white_line_and_comments(input);

			err = fscanf(input, "%s", str);
			if(err != 1)
				{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				printf("err=%d\n", err);
				exit(EXIT_FAILURE);
				}

			if(strncmp(str, "size", 5) == 0)
				{
				for(i = 0; i < STDIM; i++)
					{
					err = fscanf(input, "%d", &temp_i);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_size[i] = temp_i;
					}
				}

			else if(strncmp(str, "beta", 4) == 0)
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_beta = temp_d;
				}

			else if(strncmp(str, "htracedef", 9) == 0)
				{
				for(i = 0; i < (int)floor(NCOLOR / 2.0); i++)
					{
					err = fscanf(input, "%lf", &temp_d);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_h[i] = temp_d;
					}
				}

			else if(strncmp(str, "theta", 5) == 0)
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_theta = temp_d;
				}

			else if(strncmp(str, "sample", 6) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_sample = temp_i;
				}

			else if(strncmp(str, "thermal", 7) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_thermal = temp_i;
				}

			else if(strncmp(str, "overrelax", 9) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_overrelax = temp_i;
				}

			else if(strncmp(str, "measevery", 9) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_measevery = temp_i;
				}

			else if(strncmp(str, "start", 5) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_start = temp_i;
				}

			else if(strncmp(str, "saveconf_back_every", 19) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_saveconf_back_every = temp_i;
				}

			else if(strncmp(str, "saveconf_analysis_every", 23) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_saveconf_analysis_every = temp_i;
				}

			else if(strncmp(str, "epsilon_metro", 13) == 0)
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_epsilon_metro = temp_d;
				}

			else if(strncmp(str, "coolsteps", 9) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_coolsteps = temp_i;
				}

			else if(strncmp(str, "coolrepeat", 10) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_coolrepeat = temp_i;
				}

			else if(strncmp(str, "plaquette_meas", 14) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_plaquette_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: plaquette_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "clover_energy_meas", 18) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_clover_energy_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: clover_energy_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "charge_meas", 11) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_charge_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: charge_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "polyakov_meas", 13) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_polyakov_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: polyakov_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "chi_prime_meas", 14) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_chi_prime_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: chi_prime_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "topcharge_tprof_meas", 20) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 1) || (temp_i == 0)) param->d_topcharge_tprof_meas = temp_i;
				else
					{
					fprintf(stderr, "Error: topcharge_tprof_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "gfstep", 6) == 0) // integration step
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d > 0) param->d_gfstep = temp_d;
				else
					{
					fprintf(stderr, "Error: gfstep must be positive in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "num_gfsteps", 11) == 0) // number of integration steps
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_i > -1) param->d_ngfsteps = temp_i;
				else
					fprintf(stderr, "Error: num_gfsteps must be non negative in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "gf_meas_each", 12) == 0) // number of integration steps
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_i > 0) param->d_gf_meas_each = temp_i;
				else
					fprintf(stderr, "Error: gf_meas_each must be positive in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "agf_length", 10) == 0) // length of adaptive gradflow evolution
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d >= 0) param->d_agf_length = temp_d;
				else
					fprintf(stderr, "Error: agf_length must be non-negative in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "agf_meas_each", 13) == 0) // time interval of adaptive gradflow measures
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d > 0) param->d_agf_meas_each = temp_d;
				else
					fprintf(stderr, "Error: agf_meas_each must be positive in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "agf_step", 8) == 0) // initial size of integration with adaptive step
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d > 0) param->d_agf_step = temp_d;
				else
					fprintf(stderr, "Error: agf_step must be positive in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "agf_delta", 9) == 0) // error threshold for adaptive gradflow integration
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d > 0) param->d_agf_delta = temp_d;
				else
					fprintf(stderr, "Error: adf_delta must be positive in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "agf_time_bin", 12) == 0) // error threshold for adaptive gradflow time of measure
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if(temp_d >= 0) param->d_agf_time_bin = temp_d;
				else
					fprintf(stderr, "Error: adf_time_bin must be non-negative in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				}

			else if(strncmp(str, "multihit", 8) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_multihit = temp_i;
				}

			else if(strncmp(str, "ml_step", 7) == 0)
				{
				for(i = 0; i < NLEVELS; i++)
					{
					err = fscanf(input, "%d", &temp_i);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_ml_step[i] = temp_i;
					}
				}

			else if(strncmp(str, "ml_upd", 6) == 0)
				{
				for(i = 0; i < NLEVELS; i++)
					{
					err = fscanf(input, "%d", &temp_i);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_ml_upd[i] = temp_i;
					}
				}

			else if(strncmp(str, "ml_level0_repeat", 16) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_ml_level0_repeat = temp_i;
				}

			else if(strncmp(str, "dist_poly", 9) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_dist_poly = temp_i;
				}

			else if(strncmp(str, "transv_dist", 11) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_trasv_dist = temp_i;
				}

			else if(strncmp(str, "plaq_dir", 8) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_plaq_dir[0] = temp_i;

				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_plaq_dir[1] = temp_i;
				}

			else if(strncmp(str, "conf_file", 9) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_conf_file, temp_str);
				}

			else if(strncmp(str, "data_file", 9) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_data_file, temp_str);
				}

			else if(strncmp(str, "chiprime_data_file", 18) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_chiprime_file, temp_str);
				}

			else if(strncmp(str, "topcharge_tprof_file", 20) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_topcharge_tprof_file, temp_str);
				}

			else if(strncmp(str, "log_file", 8) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_log_file, temp_str);
				}

			else if(strncmp(str, "ml_file", 7) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_ml_file, temp_str);
				}

			else if(strncmp(str, "randseed", 8) == 0)
				{
				err = fscanf(input, "%u", &temp_ui);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_randseed = temp_ui;
				}

			else if(strncmp(str, "defect_dir", 10) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				if((temp_i == 0) || (temp_i == 1) || (temp_i == 2) || (temp_i == 3)) { param->d_defect_dir = temp_i; }
				else
					{
					fprintf(stderr, "Error in reading the file %s, defect_dir must be either 0,1,2 or 3 (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}

			else if(strncmp(str, "defect_size", 11) == 0)
				{
				for(i = 0; i < STDIM - 1; i++)
					{
					err = fscanf(input, "%d", &temp_i);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_L_defect[i] = temp_i;
					}
				}

			else if(strncmp(str, "N_replica_pt", 12) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_N_replica_pt = temp_i;
				if(posix_memalign((void **)&(param->d_pt_bound_cond_coeff), (size_t)DOUBLE_ALIGN, (size_t)param->d_N_replica_pt * sizeof(double)) != 0)
					{
					fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				for(i = 0; i < param->d_N_replica_pt; i++)
					{
					err = fscanf(input, "%lf", &temp_d);
					if(err != 1)
						{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					param->d_pt_bound_cond_coeff[i] = temp_d;
					}
				}

			else if(strncmp(str, "swap_acc_file", 13) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_swap_acc_file, temp_str);
				}

			else if(strncmp(str, "swap_track_file", 15) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_swap_tracking_file, temp_str);
				}

			else if(strncmp(str, "hierarc_upd", 11) == 0)
				{
				err = fscanf(input, "%d", &temp_i);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_N_hierarc_levels = temp_i;
				if(param->d_N_hierarc_levels > 0)
					{
					if(posix_memalign((void **)&(param->d_L_rect), (size_t)INT_ALIGN, (size_t)param->d_N_hierarc_levels * sizeof(int)) != 0)
						{
						fprintf(stderr, "Problems allocating an array of ints! (%s, %d)\n", __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					for(i = 0; i < param->d_N_hierarc_levels; i++)
						{
						err = fscanf(input, "%d", &temp_i);
						if(err != 1)
							{
							fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
							exit(EXIT_FAILURE);
							}
						param->d_L_rect[i] = temp_i;
						}
					if(posix_memalign((void **)&(param->d_N_sweep_rect), (size_t)INT_ALIGN, (size_t)param->d_N_hierarc_levels * sizeof(int)) != 0)
						{
						fprintf(stderr, "Problems allocating an array of ints! (%s, %d)\n", __FILE__, __LINE__);
						exit(EXIT_FAILURE);
						}
					for(i = 0; i < param->d_N_hierarc_levels; i++)
						{
						err = fscanf(input, "%d", &temp_i);
						if(err != 1)
							{
							fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
							exit(EXIT_FAILURE);
							}
						param->d_N_sweep_rect[i] = temp_i;
						}
					} // closes if( num_hierarc_levels > 0 )
				}

			else if(strncmp(str, "multicanonic_acc_file", 21) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_multicanonic_acc_file, temp_str);
				}

			else if(strncmp(str, "topo_potential_file", 19) == 0)
				{
				err = fscanf(input, "%s", temp_str);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				strcpy(param->d_topo_potential_file, temp_str);
				}

			else if(strncmp(str, "grid_step ", 9) == 0)
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_grid_step = temp_d;
				}

			else if(strncmp(str, "grid_max", 8) == 0)
				{
				err = fscanf(input, "%lf", &temp_d);
				if(err != 1)
					{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				param->d_grid_max = temp_d;
				}

			else
				{
				fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, in_file, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}

			remove_white_line_and_comments(input);

			// check if the read line is the last one
			temp_i = getc(input);
			if(temp_i == EOF) { end = 0; }
			else
				{
				ungetc(temp_i, input);
				}
			}

		fclose(input);

		// VARIOUS CHECKS
		if(param->d_ml_step[0] != 0)
			{
			if(param->d_size[0] % param->d_ml_step[0] || param->d_size[0] < param->d_ml_step[0])
				{
				fprintf(stderr, "Error: size[0] has to be divisible by ml_step[0] and satisfy ml_step[0]<=size[0] (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			for(i = 1; i < NLEVELS; i++)
				{
				if(param->d_ml_step[i - 1] % param->d_ml_step[i] || param->d_ml_step[i - 1] <= param->d_ml_step[i])
					{
					fprintf(stderr, "Error: ml_step[%d] has to be divisible by ml_step[%d] and larger than it (%s, %d)\n", i - 1, i, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}
			if(param->d_ml_step[NLEVELS - 1] == 1)
				{
				fprintf(stderr, "Error: ml_step[%d] has to be larger than 1 (%s, %d)\n", NLEVELS - 1, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			}

		#ifdef OPENMP_MODE
		for(i = 0; i < STDIM; i++)
			{
			temp_i = param->d_size[i] % 2;
			if(temp_i != 0)
				{
				fprintf(stderr, "Error: size[%d] is not even.\n", i);
				fprintf(stderr, "When using OpenMP all the sides of the lattice have to be even! (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			}
		#endif

		err = 0;
		for(i = 0; i < STDIM; i++)
			{
			if(param->d_size[i] == 1) { err = 1; }
			}
		if(err == 1) { fprintf(stderr, "Error: all sizes has to be larger than 1: the totally reduced case is not implemented! (%s, %d)\n", __FILE__, __LINE__); }

		// various checks on parallel tempering parameters
		if(param->d_L_defect[0] > param->d_size[0])
			{
			fprintf(stderr, "Error: defect's t-length is greater than lattice's t-length (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(param->d_L_defect[1] > param->d_size[2])
			{
			fprintf(stderr, "Error: defect's y-length is greater than lattice's y-length (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(param->d_L_defect[2] > param->d_size[3])
			{
			fprintf(stderr, "Error: defect's z-length is greater than lattice's z-length (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(param->d_N_replica_pt < 1)
			{
			fprintf(stderr, "Error: number of replica of parallel tempering must be greater than 0 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}

		// check on gradflow parameters
		if(param->d_agf_meas_each <= param->d_agf_time_bin)
			{
			fprintf(stderr, "Error: agf_meas_each must be greater than agf_time_bin (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(param->d_agf_meas_each <= MIN_VALUE)
			{
			fprintf(stderr, "Error: agf_meas_each must be greater than MIN_VALUE in /include/macro.h (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}

		init_derived_constants(param);
		}
	}

void init_derived_constants(GParam *param)
	{
	int i;

	// derived constants
	param->d_volume = 1;
	for(i = 0; i < STDIM; i++) { (param->d_volume) *= (param->d_size[i]); }

	param->d_space_vol = 1;
	// direction 0 is time
	for(i = 1; i < STDIM; i++) { (param->d_space_vol) *= (param->d_size[i]); }

	param->d_inv_vol = 1.0 / ((double)param->d_volume);
	param->d_inv_space_vol = 1.0 / ((double)param->d_space_vol);

	// volume of the defect
	param->d_volume_defect = 1;
	for(i = 0; i < STDIM - 1; i++) { param->d_volume_defect *= param->d_L_defect[i]; }

	// number of grid points (multicanonic only)
	param->d_n_grid = (int)((2.0 * param->d_grid_max / param->d_grid_step) + 1.0);
	}

// initialize data file
void init_data_file(FILE **dataf, FILE **chiprimef, FILE **topchar_tprof_f, GParam const *const param)
	{
	int i;

	if(param->d_start == 2)
		{
		// open std data file (plaquette, polyakov, topological charge)
		*dataf = fopen(param->d_data_file, "r");
		if(*dataf != NULL) // file exists
			{
			fclose(*dataf);
			*dataf = fopen(param->d_data_file, "a");
			}
		else
			{
			*dataf = fopen(param->d_data_file, "w");
			fprintf(*dataf, "%d ", STDIM);
			for(i = 0; i < STDIM; i++) fprintf(*dataf, "%d ", param->d_size[i]);
			fprintf(*dataf, "\n");
			}
		// open chi prime data file
		if(param->d_chi_prime_meas == 1)
			{
			*chiprimef = fopen(param->d_chiprime_file, "r");
			if(*chiprimef != NULL) // file exists
				{
				fclose(*chiprimef);
				*chiprimef = fopen(param->d_chiprime_file, "a");
				}
			else
				{
				*chiprimef = fopen(param->d_chiprime_file, "w");
				fprintf(*chiprimef, "%d ", STDIM);
				for(i = 0; i < STDIM; i++) fprintf(*chiprimef, "%d ", param->d_size[i]);
				fprintf(*chiprimef, "\n");
				}
			}
		else
			{
			(void)chiprimef;
			}
		// open topocharge_tprof data file
		if(param->d_topcharge_tprof_meas == 1)
			{
			*topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "r");
			if(*topchar_tprof_f != NULL) // file exists
				{
				fclose(*topchar_tprof_f);
				*topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "a");
				}
			else
				{
				*topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "w");
				fprintf(*topchar_tprof_f, "%d ", STDIM);
				for(i = 0; i < STDIM; i++) fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
				fprintf(*topchar_tprof_f, "\n");
				}
			}
		else
			{
			(void)topchar_tprof_f;
			}
		}
	else
		{
		// open std data file
		*dataf = fopen(param->d_data_file, "w");
		fprintf(*dataf, "%d ", STDIM);
		for(i = 0; i < STDIM; i++) { fprintf(*dataf, "%d ", param->d_size[i]); }
		fprintf(*dataf, "\n");
		// open chi prime data file
		if(param->d_chi_prime_meas == 1)
			{
			*chiprimef = fopen(param->d_chiprime_file, "w");
			fprintf(*chiprimef, "%d ", STDIM);
			for(i = 0; i < STDIM; i++) fprintf(*chiprimef, "%d ", param->d_size[i]);
			fprintf(*chiprimef, "\n");
			}
		else
			{
			(void)chiprimef;
			}
		// open topocharge_tprof data file
		if(param->d_topcharge_tprof_meas == 1)
			{
			*topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "w");
			fprintf(*topchar_tprof_f, "%d ", STDIM);
			for(i = 0; i < STDIM; i++) fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
			fprintf(*topchar_tprof_f, "\n");
			}
		else
			{
			(void)topchar_tprof_f;
			}
		}
	fflush(*dataf);
	if(param->d_chi_prime_meas == 1) fflush(*chiprimef);
	else
		{
		(void)chiprimef;
		}
	if(param->d_topcharge_tprof_meas == 1) fflush(*topchar_tprof_f);
	else
		{
		(void)topchar_tprof_f;
		}
	}

// free allocated memory for hierarc update parameters
void free_hierarc_params(GParam *param)
	{
	if(param->d_N_hierarc_levels == 0)
		{
		(void)param; // to avoid compiler warning about unused variable
		}
	else
		{
		free(param->d_L_rect);
		free(param->d_N_sweep_rect);
		}
	}

// print simulation parameters aux

void print_configuration_parameters(FILE *fp)
	{
	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif
	
	if(endian()==0) fprintf(fp, "Little endian machine\n\n");
	else fprintf(fp, "Big endian machine\n\n");
	}

void print_pt_parameters(FILE *fp, GParam const * const param)
	{
	int i;
	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
	fprintf(fp, "defect: %d", param->d_L_defect[0]);
	for(i=1; i<STDIM-1; i++) fprintf(fp, "x%d", param->d_L_defect[i]);
	fprintf(fp, "\n\n");
	
	fprintf(fp,"number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
	fprintf(fp,"boundary condition constants: ");
	for(i=0;i<param->d_N_replica_pt;i++) fprintf(fp,"%lf ",param->d_pt_bound_cond_coeff[i]);
	fprintf(fp,"\n");
	
	fprintf(fp,"number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
	if(param->d_N_hierarc_levels>0)
		{
		fprintf(fp,"extention of rectangles: ");
		for(i=0;i<param->d_N_hierarc_levels;i++)
			{
			fprintf(fp,"%d ", param->d_L_rect[i]);
			}
		fprintf(fp,"\n");
		fprintf(fp,"number of sweeps per hierarchical level: ");
		for(i=0;i<param->d_N_hierarc_levels;i++)
			{
			fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
			}
		}
	fprintf(fp,"\n\n");
	}

void print_multicanonic_parameters(FILE *fp, GParam const * const param)
	{
	fprintf(fp,"Multicanonic topo-potential read from file %s\nPotential defined on a grid with step=%.10lf and max=%.10lf\n", param->d_topo_potential_file, param->d_grid_step, param->d_grid_max);
	fprintf(fp,"\n\n");
	}

void print_simul_parameters(FILE *fp, GParam const * const param)
	{
	int i;
	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);
	
	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i=1; i<STDIM; i++) fprintf(fp, "x%d", param->d_size[i]);
	fprintf(fp, "\n\n");
	
	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");
	
	fprintf(fp, "sample:	%d\n", param->d_sample);
	fprintf(fp, "thermal:	%d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");
	
	fprintf(fp, "plaquette_meas: %d\n", param->d_plaquette_meas);
	fprintf(fp, "clover_energy_meas: %d\n", param->d_clover_energy_meas);
	fprintf(fp, "charge_meas: %d\n", param->d_charge_meas);
	fprintf(fp, "polyakov_meas: %d\n", param->d_polyakov_meas);
	fprintf(fp, "chi_prime_meas: %d\n", param->d_chi_prime_meas);
	fprintf(fp, "topcharge_tprof_meas: %d\n", param->d_topcharge_tprof_meas);
	fprintf(fp, "\n");

	fprintf(fp, "start:                     %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:       %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every:   %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");
	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");
	}

void print_adaptive_gradflow_parameters(FILE *fp, GParam const * const param)
	{
	fprintf(fp, "agf_length     %lf\n", param->d_agf_length);
	fprintf(fp, "agf_step:      %lf\n", param->d_agf_step);
	fprintf(fp, "agf_meas_each  %lf\n", param->d_agf_meas_each);
	fprintf(fp, "agf_delta      %e\n",  param->d_agf_delta);
	fprintf(fp, "\n");
	}

void print_gradflow_parameters(FILE *fp, GParam const * const param)
	{
	fprintf(fp, "gfstep:        %lf\n", param->d_gfstep);
	fprintf(fp, "num_gfsteps    %d\n",  param->d_ngfsteps);
	fprintf(fp, "gf_meas_each   %d\n",    param->d_gf_meas_each);
	fprintf(fp, "\n");
	}

void print_cooling_parameters(FILE *fp, GParam const * const param)
	{
	fprintf(fp, "coolsteps:     %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:    %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");
	}

void print_multilevel_parameters(FILE *fp, GParam const * const param)
	{
	int i;
	fprintf(fp, "multihit:	%d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i=0; i<NLEVELS; i++)
		{
		fprintf(fp, "%d ", param->d_ml_step[i]);
		}
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i=0; i<NLEVELS; i++)
		{
		fprintf(fp, "%d ", param->d_ml_upd[i]);
		}
	fprintf(fp, "\n\n");
	}

void print_metro_parameters(FILE *fp, GParam const * const param, double acc)
	{
	fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
	fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
	fprintf(fp, "\n");
	}

// print simulation parameters

void print_parameters_local(GParam const *const param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+-----------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local |\n");
	fprintf(fp, "+-----------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");

	fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

void print_parameters_local_pt_multicanonic(GParam const *const param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+---------------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_pt_multicanonic |\n");
	fprintf(fp, "+---------------------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
	fprintf(fp, "defect: %d", param->d_L_defect[0]);
	for(i = 1; i < STDIM - 1; i++) { fprintf(fp, "x%d", param->d_L_defect[i]); }
	fprintf(fp, "\n\n");
	fprintf(fp, "number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
	fprintf(fp, "boundary condition constants: ");
	for(i = 0; i < param->d_N_replica_pt; i++) fprintf(fp, "%lf ", param->d_pt_bound_cond_coeff[i]);
	fprintf(fp, "\n");
	fprintf(fp, "number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
	if(param->d_N_hierarc_levels > 0)
		{
		fprintf(fp, "extention of rectangles: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) { fprintf(fp, "%d ", param->d_L_rect[i]); }
		fprintf(fp, "\n");
		fprintf(fp, "number of sweeps per hierarchical level: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) { fprintf(fp, "%d ", param->d_N_sweep_rect[i]); }
		}
	fprintf(fp, "\n\n");

	fprintf(fp, "Multicanonic topo-potential read from file %s\nPotential defined on a grid with step=%.10lf and max=%.10lf\n", param->d_topo_potential_file, param->d_grid_step, param->d_grid_max);

	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");

	fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

void print_parameters_local_pt(GParam const *const param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+--------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_pt |\n");
	fprintf(fp, "+--------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
	fprintf(fp, "defect: %d", param->d_L_defect[0]);
	for(i = 1; i < STDIM - 1; i++) { fprintf(fp, "x%d", param->d_L_defect[i]); }
	fprintf(fp, "\n\n");
	fprintf(fp, "number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
	fprintf(fp, "boundary condition constants: ");
	for(i = 0; i < param->d_N_replica_pt; i++) fprintf(fp, "%lf ", param->d_pt_bound_cond_coeff[i]);
	fprintf(fp, "\n");
	fprintf(fp, "number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
	if(param->d_N_hierarc_levels > 0)
		{
		fprintf(fp, "extention of rectangles: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) { fprintf(fp, "%d ", param->d_L_rect[i]); }
		fprintf(fp, "\n");
		fprintf(fp, "number of sweeps per hierarchical level: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) { fprintf(fp, "%d ", param->d_N_sweep_rect[i]); }
		}
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");

	fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_local_pt_agf(GParam const * const param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	double diff_sec;

	fp=fopen(param->d_log_file, "w");
	fprintf(fp, "+------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_pt_agf |\n");
	fprintf(fp, "+------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	print_simul_parameters(fp, param);
	print_adaptive_gradflow_parameters(fp, param);
	print_cooling_parameters(fp, param);
	
	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
	fprintf(fp, "\n");
	
	fclose(fp);
	}

// print simulation parameters
void print_parameters_polycorr_long(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+-------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_polycorr_long |\n");
	fprintf(fp, "+-------------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "\n");

	fprintf(fp, "multihit:   %d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_step[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_upd[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
	fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_polycorr(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+--------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_polycorr |\n");
	fprintf(fp, "+--------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif

	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "\n");

	fprintf(fp, "multihit:   %d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_step[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_upd[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "dist_poly:  %d\n", param->d_dist_poly);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_t0(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+--------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_t0 |\n");
	fprintf(fp, "+--------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "gfstep:    %lf\n", param->d_gfstep);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_gf(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+-------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_gradient_flow |\n");
	fprintf(fp, "+-------------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "gfstep:        %lf\n", param->d_gfstep);
	fprintf(fp, "num_gfsteps    %d\n", param->d_ngfsteps);
	fprintf(fp, "gf_meas_each   %d\n", param->d_gf_meas_each);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_agf(GParam * param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp=fopen(param->d_log_file, "w");
	fprintf(fp, "+----------------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_adaptive_gradient_flow |\n");
	fprintf(fp, "+----------------------------------------------------------+\n\n");
	
	print_configuration_parameters(fp);
	
	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);
	
	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i=1; i<STDIM; i++)
		{
		fprintf(fp, "x%d", param->d_size[i]);
		}
	fprintf(fp, "\n\n");
	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");
	
	print_adaptive_gradflow_parameters(fp, param);
	
	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
	fprintf(fp, "\n");
	
	fclose(fp);
	}

// print simulation parameters for the tracedef case
void print_parameters_tracedef(GParam const *const param, time_t time_start, time_t time_end, double acc)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+--------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tracedef |\n");
	fprintf(fp, "+--------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "h: %.10lf ", param->d_h[0]);
	for(i = 1; i < (int)floor(NCOLOR / 2.0); i++) { fprintf(fp, "%.10lf ", param->d_h[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");

	fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
	fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
	fprintf(fp, "\n");

	fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_tube_disc(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+---------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tube_disc |\n");
	fprintf(fp, "+---------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "\n");

	fprintf(fp, "multihit:   %d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_step[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_upd[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
	fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
	fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_tube_conn(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+---------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tube_conn |\n");
	fprintf(fp, "+---------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "\n");

	fprintf(fp, "multihit:   %d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_step[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_upd[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
	fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
	fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print simulation parameters
void print_parameters_tube_conn_long(GParam *param, time_t time_start, time_t time_end)
	{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+--------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tube_conn_long |\n");
	fprintf(fp, "+--------------------------------------------------+\n\n");

	#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) { fprintf(fp, "x%d", param->d_size[i]); }
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "\n");

	fprintf(fp, "multihit:   %d\n", param->d_multihit);
	fprintf(fp, "levels for multileves: %d\n", NLEVELS);
	fprintf(fp, "multilevel steps: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_step[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "updates for levels: ");
	for(i = 0; i < NLEVELS; i++) { fprintf(fp, "%d ", param->d_ml_upd[i]); }
	fprintf(fp, "\n");
	fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
	fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
	fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
	fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if(endian() == 0) { fprintf(fp, "Little endian machine\n\n"); }
	else
		{
		fprintf(fp, "Big endian machine\n\n");
		}

	fclose(fp);
	}

// print template input aux

void print_template_volume_parameters(FILE *fp)
	{
	fprintf(fp,"size 12 4 4 12  # Nt Nx Ny Nz\n");
	fprintf(fp,"\n");
	}

void print_template_simul_parameters(FILE *fp)
	{
	fprintf(fp,"# Simulations parameters\n");
	fprintf(fp, "beta  6.4881\n");
	fprintf(fp,"\n");
	fprintf(fp, "sample     10\n");
	fprintf(fp, "thermal    0\n");
	fprintf(fp, "overrelax  5\n");
	fprintf(fp, "measevery  1\n");
	fprintf(fp,"\n");
	
	fprintf(fp, "start                    0  # 0=all links to identity  1=random  2=from saved configuration\n");
	fprintf(fp, "saveconf_back_every      5  # if 0 does not save, else save backup configurations every ... updates\n");
	fprintf(fp, "saveconf_analysis_every  5  # if 0 does not save, else save configurations for analysis every ... updates\n");
	fprintf(fp, "\n");
	
	fprintf(fp, "randseed 0    # (0=time)\n");
	fprintf(fp, "\n");
	
	fprintf(fp, "# Observables to measure\n");
	fprintf(fp, "plaquette_meas        0  # 1=YES, 0=NO\n");
	fprintf(fp, "clover_energy_meas    1  # 1=YES, 0=NO\n");
	fprintf(fp, "charge_meas           1  # 1=YES, 0=NO\n");
	fprintf(fp, "polyakov_meas         0  # 1=YES, 0=NO\n");
	fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
	fprintf(fp, "topcharge_tprof_meas  0  # 1=YES, 0=NO\n");
	fprintf(fp,"\n");
	}

void print_template_pt_parameters(FILE *fp)
	{
	fprintf(fp,"# Parallel tempering parameters\n");
	fprintf(fp,"defect_dir    0             # choose direction of defect boundary: 0->t, 1->x, 2->y, 3->z\n");
	fprintf(fp,"defect_size   2 2 2         # size of the defect (order: y-size z-size t-size)\n");
	fprintf(fp,"N_replica_pt  2    1.0 0.0  # number of parallel tempering replica ____ boundary conditions coefficients\n");
	fprintf(fp,"\n");
	fprintf(fp,"# Hierarchical update parameters\n");
	fprintf(fp,"# Order: num of hierarc levels ____ extension of rectangles ____ num of sweeps per rectangle\n");
	fprintf(fp,"hierarc_upd 2    2 1    1 1\n");
	fprintf(fp,"\n");
	}

void print_template_adaptive_gradflow_parameters(FILE *fp)
	{
	fprintf(fp, "# For adaptive gradient flow evolution\n");
	fprintf(fp, "agf_length       10    # total integration time for adaptive gradient flow\n");
	fprintf(fp, "agf_step       0.01    # initial integration step for adaptive gradient flow\n");
	fprintf(fp, "agf_meas_each     1    # time interval between measures during adaptive gradient flow\n");
	fprintf(fp, "agf_delta     0.001    # error threshold on gauge links for adaptive gradient flow\n");
	fprintf(fp, "agf_time_bin  0.001    # error threshold on time of measures for adaptive gradient flow\n");
	fprintf(fp, "\n");
	}

void print_template_gradflow_parameters(FILE *fp)
	{
	fprintf(fp, "# For gradient flow evolution\n");
	fprintf(fp, "gfstep      0.01    # integration step for gradient flow\n");
	fprintf(fp, "num_gfsteps 100     # number of integration steps for gradient flow\n");
	fprintf(fp, "gf_meas_each 5      # compute observables every <gfstep_each> integration steps during the gradient flow\n");
	fprintf(fp, "\n");
	}

void print_template_cooling_parameters(FILE *fp)
	{
	fprintf(fp, "# For cooling\n");
	fprintf(fp, "coolsteps             3  # number of cooling steps to be used\n");
	fprintf(fp, "coolrepeat            5  # number of times 'coolsteps' are repeated\n");
	fprintf(fp, "\n");
	}

void print_template_metro_parameters(FILE *fp)
	{
	fprintf(fp, "epsilon_metro    0.25  #distance from the identity of the random matrix for metropolis\n");
	fprintf(fp, "\n");
	}

void print_template_multicanonic_parameters(FILE *fp)
	{
	fprintf(fp,"# Multicanonic parameters\n");
	fprintf(fp,"grid_step             0.05\n");
	fprintf(fp,"grid_max              3.0\n");
	fprintf(fp,"topo_potential_file   topo_potential\n");
	fprintf(fp, "multicanonic_acc_file   multicanonic_acc.dat\n");
	fprintf(fp, "\n");
	}

void print_template_multilevel_parameters(FILE *fp)
	{
    fprintf(fp, "# For multilevel\n");
    fprintf(fp, "multihit         10  # number of multihit step\n");
    fprintf(fp, "ml_step           2  # timeslices for multilevel (from largest to smallest)\n");
    fprintf(fp, "ml_upd           10  # number of updates for various levels\n");
	fprintf(fp, "ml_file      ml.dat  # multilevel output file\n");
	fprintf(fp, "\n");
	}

void print_template_output_parameters(FILE *fp)
	{
	fprintf(fp, "# Output files\n");
	fprintf(fp, "conf_file             conf.dat\n");
	fprintf(fp, "data_file             dati.dat\n");
	fprintf(fp, "chiprime_data_file    chi_prime_cool.dat\n");
	fprintf(fp, "topcharge_tprof_file  topo_tprof_cool.dat\n");
	fprintf(fp, "log_file              log.dat\n");
	fprintf(fp, "swap_acc_file         swap_acc.dat\n");
	fprintf(fp, "swap_track_file       swap_track.dat\n");
	fprintf(fp, "\n");
	}

// print program details

void print_authors(int parallel_tempering)
	{
	printf("\n");
	if(parallel_tempering==1)
		{
		printf("SU(N) Hasenbusch Parallel Tempering implemented by Claudio Bonanno (claudiobonanno93@gmail.com)\n");
		printf("within yang-mills package\n\n");
		}
	
	printf("Details about yang-mills package:\n");
	printf("\tPackage %s version: %s\n", PACKAGE_NAME, PACKAGE_VERSION);
	printf("\tAuthor: Claudio Bonati %s\n\n", PACKAGE_BUGREPORT);
	}

void print_compilation_details()
	{
	printf("Compilation details:\n");
	printf("\tN_c (number of colors): %d\n", NCOLOR);
	printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
	printf("\tNum_levels (number of levels): %d\n", NLEVELS);
	printf("\n");
	printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
	printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));
	
	#ifdef DEBUG
	printf("\n\tDEBUG mode\n");
	#endif
	
	#ifdef OPENMP_MODE
	printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
	#endif
	
	#ifdef THETA_MODE
	printf("\n\tusing imaginary theta\n");
	#endif
	
	#ifdef OPT_MULTIHIT
	printf("\tcompiled for multihit optimization\n");
	#endif
	
	#ifdef OPT_MULTILEVEL
	printf("\tcompiled for multilevel optimization\n");
	#endif
	
	printf("\n");
	
	#ifdef __INTEL_COMPILER
	printf("\tcompiled with icc\n");
	#elif defined(__clang__)
	printf("\tcompiled with clang\n");
	#elif defined( __GNUC__ )
	printf("\tcompiled with gcc version: %d.%d.%d\n",
		__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	#endif
	printf("\n");
	}

#endif
