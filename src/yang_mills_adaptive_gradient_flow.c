#ifndef YM_AGF_C
#define YM_AGF_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#ifdef OPENMP_MODE
#include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"
#include"../include/timing.h"

void real_main(char *in_file, long step, long stop_index)
	{
	Gauge_Conf GC;
	Geometry geo;
	GParam param;
	Meas_Utils meas_aux;
	Time_Utils timers;
	
	// to disable nested parallelism
	#ifdef OPENMP_MODE
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// read input file
	readinput(in_file, &param);

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// this code has to start from saved conf.
	param.d_start = 2;
	
	// not to overwrite files of runs with online gradient flow 
	strcat(param.d_data_file, "_agf");
	strcat(param.d_energydensity_file, "_agf");
	strcat(param.d_chiprime_file, "_agf");
	strcat(param.d_topcharge_tcorr_file, "_agf");
	strcat(param.d_log_file, "_agf");

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);
	
	// find and init first conf
	while(init_gauge_conf_step(&GC, &param, step++) == 0 && step <= stop_index);
	if (step > stop_index)
		{
		fprintf(stderr, "No configuration found up to update index %ld (%s, %d)\n", stop_index, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	stop_timer(&(timers.init_timer));

	// perform measures and load next conf
	while(step <= stop_index)
		{
		start_timer(&(timers.step_timer));
		
		start_timer(&(timers.meas_timer));
		perform_measures_localobs_with_adaptive_gradflow(&GC, &geo, &param, &meas_aux);
		stop_timer(&(timers.meas_timer));
		
		while(read_gauge_conf_step(&GC, &param, step++) == 0 && step <= stop_index);
		
		stop_timer(&(timers.step_timer));
		if (wall_time_check(&timers) == 1) break;
		}
	
	stop_timer(&(timers.prog_timer));
	
	// free gauge conf
	free_gauge_conf(&GC, &param);
	
	// free meas utils
	free_meas_utils(meas_aux, &param, 0);
	
	// free geometry
	free_geometry(&geo, &param);

	// print simulation details
	print_parameters_agf(&param, &timers);
	}

void print_template_input(void)
	{
	FILE *fp;

	fp = fopen("template_input.example", "w");

	if(fp == NULL)
		{
		fprintf(stderr, "Error in opening the file template_input.example (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		print_template_volume_parameters(fp);
		print_template_pt_parameters(fp);
		print_template_twist_parameters(fp);
		#ifdef MULTICANONICAL_MODE
		print_template_multicanonic_parameters(fp);
		#endif
		print_template_simul_parameters(fp);
		print_template_adaptive_gradflow_parameters(fp);
		print_template_output_parameters(fp);
		fclose(fp);
		}
	}

int main(int argc, char **argv)
	{
	char in_file[500];

	if(argc != 4)
		{
		int parallel_tempering = 0;
		int twisted_bc = 1;
		print_authors(parallel_tempering, twisted_bc);
		
		printf("Usage: %s input_file start_index stop_index\n\n", argv[0]);
		
		print_compilation_details();
		print_template_input();
		
		return EXIT_SUCCESS;
		}
	else
		{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
			{
			fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
			}
		else
			{
			strcpy(in_file, argv[1]);
			}
		}
	
	long start_index = strtol(argv[2], NULL, 10);
	long stop_index = strtol(argv[3], NULL, 10);
	if(start_index < 0) fprintf(stderr, "argument 'start_index' must be a non-negative integer\n");
	if(stop_index < 0) fprintf(stderr, "argument 'stop_index' must be a non-negative integer\n");
	if(stop_index < start_index) fprintf(stderr, "argument 'stop_index' must not be smaller than argument 'start_index' \n");
	real_main(in_file, start_index, stop_index);

	return EXIT_SUCCESS;
	}

#endif