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
	strcat(param.d_energy_slices_file, "_agf");
	strcat(param.d_charge_slices_file, "_agf");
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
	REQUIRE(step <= stop_index, "no configuration found up to update index %ld", stop_index);
	stop_timer(&(timers.init_timer));

	// perform measures and load next conf
	while(step <= stop_index)
		{
		start_timer(&(timers.step_timer));

		start_timer(&(timers.meas_timer));
		perform_measures_localobs(&GC, &geo, &param, &meas_aux);
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
	FILE *fp = fopen("template_input.example", "w");
	REQUIRE(fp != NULL, "failed to open template_input.example");

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

int main(int argc, char **argv)
	{
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
	
	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");
	
	char *end;

	long start_index = strtol(argv[2], &end, 10);
	REQUIRE(*end == '\0', "start_index is not a valid integer");

	long stop_index = strtol(argv[3], &end, 10);
	REQUIRE(*end == '\0', "stop_index is not a valid integer");

	REQUIRE(start_index >= 0, "start_index must be non-negative");
	REQUIRE(stop_index >= 0, "stop_index must be non-negative");
	REQUIRE(stop_index >= start_index, "stop_index must be >= start_index");

	real_main(argv[1], start_index, stop_index);

	return EXIT_SUCCESS;
	}

#endif