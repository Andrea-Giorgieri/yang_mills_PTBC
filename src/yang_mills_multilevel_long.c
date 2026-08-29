#ifndef YM_MULTILEVEL_LONG_C
#define YM_MULTILEVEL_LONG_C

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

	// switch to long ml_obs for calling the appropriate functions
	switch(param.d_ml_obs)
		{
		case NONE:
			break;
		case POLYCORR:
		case POLYCORR_LONG:
			param.d_ml_obs = POLYCORR_LONG;
			break;
		case TUBE_CONN:
		case TUBE_CONN_LONG:
			param.d_ml_obs = TUBE_CONN_LONG;
			break;
		default:
			REQUIRE(0, "unsupported multilevel observable %s for long simulations\n", param.d_ml_obs_str);
		}

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

	// initialize multilevel state
	int iteration = -1;
	if(param.d_start == 2)
		{
		read_multilevel_status_from_file(&GC, &param, &iteration, param.d_ml_obs);
		}

	stop_timer(&(timers.init_timer));

	// Monte Carlo begin
	start_timer(&(timers.step_timer));
	if(iteration == -1) // only standard update, no multilevel
		{
		for(int count = 0; count < param.d_measevery; count++)
			{
			start_timer(&(timers.update_timer));
			update(&GC, &geo, &param, &acc_counters);
			stop_timer(&(timers.update_timer));
			}
		// next time, start multilevel
		iteration = 0;
		}
	else
		{
		start_timer(&(timers.update_timer));
		perform_multilevel_long_update_zero(&GC, &geo, &param, iteration, param.d_ml_obs);
		iteration += 1;
		stop_timer(&(timers.update_timer));
		if(iteration == param.d_ml_level0_repeat)
			{
			// end of multilevel, only standard updates next time
			start_timer(&(timers.meas_timer));
			perform_multilevel_update_and_measures(&GC, &geo, &param, &meas_aux, param.d_ml_obs);
			stop_timer(&(timers.meas_timer));
			iteration = -1;
			}
		}
	stop_timer(&(timers.step_timer));

	// Monte Carlo end
	stop_timer(&(timers.prog_timer));

	// save multilevel stuff
	write_multilevel_status_on_file(&GC, &param, iteration, param.d_ml_obs);

	// save configuration
	write_conf_on_file(&GC, &param);
	if(param.d_saveconf_back_every != 0)
		{
		// backup copy
		write_conf_on_file_back(&GC, &param);
		}

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);

	// print simulation details
	print_parameters_multilevel_long(&param, &timers);

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
	fprintf(fp, "ml_level0_repeat  1 # number of times level0 is repeated in long simulations\n");
	fprintf(fp, "\n");
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
