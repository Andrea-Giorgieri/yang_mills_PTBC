#ifndef YM_MULTILEVEL_C
#define YM_MULTILEVEL_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/timing.h"

void real_main(char const *in_file)
	{
	Gauge_Conf GC;
	Geometry geo;
	GParam param;
	Acc_Utils acc_counters;
	Meas_Utils meas_aux;
	Time_Utils timers;

	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// read input file
	readinput(in_file, &param);

	int size_1 = param.d_size[1];
	for(int i = 2; i < STDIM; i++)
		REQUIRE(param.d_size[i] == size_1, "all the spatial sizes must be equal");

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_geometry(&geo, &param);

	// initialize gauge configuration
	init_gauge_conf(&GC, &geo, &param);

	// initialize multilevel arrays
	allocate_multilevel_arrays(&GC, &param, param.d_ml_obs);

	// initialize meas utils
	init_meas_utils(&meas_aux, &param, 0);

	stop_timer(&(timers.init_timer));

	// if number of MC updates is set to 0, perform measures on initialized configuration
	if(param.d_sample == 0)
		{
		start_timer(&(timers.step_timer));
		start_timer(&(timers.meas_timer));
		perform_multilevel_update_and_measures(&GC, &geo, &param, &meas_aux, param.d_ml_obs);
		start_timer(&(timers.meas_timer));
		stop_timer(&(timers.step_timer));
		}

	// Monte Carlo begin
	for(int count = 0; count < param.d_sample; count++)
		{
		start_timer(&(timers.step_timer));

		// update configuration
		start_timer(&(timers.update_timer));
		update(&GC, &geo, &param, &acc_counters);
		stop_timer(&(timers.update_timer));

		// perform measures
		if(GC.update_index % param.d_measevery == 0 && GC.update_index >= param.d_thermal)
			{
			start_timer(&(timers.meas_timer));
			perform_multilevel_update_and_measures(&GC, &geo, &param, &meas_aux, param.d_ml_obs);
			stop_timer(&(timers.meas_timer));
			}

		// save configuration for backup
		if(param.d_saveconf_back_every != 0)
			{
			if(GC.update_index % param.d_saveconf_back_every == 0)
				{
				// simple
				write_conf_on_file(&GC, &param);
				// backup copy
				write_conf_on_file_back(&GC, &param);
				}
			}

		// save configuration for offline analysis
		if(param.d_saveconf_analysis_every != 0)
			{
			if(GC.update_index % param.d_saveconf_analysis_every == 0)
				{
				strcpy(name, param.d_conf_file);
				strcat(name, "_step_");
				sprintf(aux, "%ld", GC.update_index);
				strcat(name, aux);
				write_conf_on_file_with_name(&GC, &param, name);

				strcpy(name, param.d_twist_file);
				strcat(name, "_step_");
				strcat(name, aux);
				write_twist_on_file_with_name(&GC, &param, name);
				}
			}

		stop_timer(&(timers.step_timer));
		if(wall_time_check(&timers) == 1) break;
		}

	// Monte Carlo end
	stop_timer(&(timers.prog_timer));

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);

	// save configuration
	if(param.d_saveconf_back_every != 0)
		{
		write_conf_on_file(&GC, &param);
		}

	// print simulation details
	print_parameters_multilevel(&param, &timers);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free multilevel arrays
	free_multilevel_arrays(&GC, &param, param.d_ml_obs);

	// free geometry
	free_geometry(&geo, &param);
	}


void print_template_input(void)
	{
	FILE *fp = fopen("template_input.example", "w");
	REQUIRE(fp != NULL, "failed to open template_input.example");

	print_template_volume_parameters(fp);
	print_template_twist_parameters(fp);
	print_template_simul_parameters(fp);
	print_template_multilevel_parameters(fp);
	print_template_output_parameters(fp);

	fclose(fp);
	}


int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		int parallel_tempering = 0;
		int twisted_bc = 1;
		print_authors(parallel_tempering, twisted_bc);

		printf("Usage: %s input_file\n\n", argv[0]);

		print_compilation_details();
		print_template_input();

		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");

	real_main(argv[1]);

	return EXIT_SUCCESS;
	}

#endif
