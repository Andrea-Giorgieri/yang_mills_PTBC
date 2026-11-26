#ifndef YM_CONF_ANALYSIS_C
#define YM_CONF_ANALYSIS_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/device.h"
#include "../include/device_data.h"
#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/timing.h"

void real_main(char *in_file, long step, long stop_index)
	{
	Gauge_Conf GC;
	Geometry geo;
	GParam param;
	Meas_Utils meas_aux;
	Time_Utils timers;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	omp_set_max_active_levels(1);
	#endif

	// read input file
	readinput(in_file, &param);
	edit_params_meas_only(&param);

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_geometry(&geo, &param);

	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);

	// find and init first conf
	while(step < stop_index && init_gauge_conf_step(&GC, &param, step) == 0)
		step++;

	REQUIRE(step < stop_index, "no configuration found before update index %ld", stop_index);

	// copy the shared data and the configuration to the device (one configuration, no rectangles)
	init_device_utils(&param, 1);
	init_device_data(&GC, &geo, &param, NULL, &meas_aux, 1);
	#ifdef MULTICANONICAL_MODE
	GC.stored_topcharge = multicanonic_charge_conf(&GC, &geo, &param); // charge of the multicanonic weight
	#endif
	stop_timer(&(timers.init_timer));

	// perform measures and load next conf
	while(step < stop_index)
		{
		start_timer(&(timers.step_timer));

		start_timer(&(timers.meas_timer));
		perform_measures_localobs(&GC, &geo, &param, &meas_aux);
		stop_timer(&(timers.meas_timer));

		do
			step++;
		while(step < stop_index && read_gauge_conf_step(&GC, &param, step) == 0);

		// new configuration (read on the host with its copies and twist factors) to the device
		if(step < stop_index)
			{
			copy_gauge_conf_to_device(&GC, &param, GC_LATTICE | GC_COPY | GC_Z);
			#ifdef MULTICANONICAL_MODE
			GC.stored_topcharge = multicanonic_charge_conf(&GC, &geo, &param); // charge of the multicanonic weight
			#endif
			}

		stop_timer(&(timers.step_timer));
		if(wall_time_check(&timers) == 1) break;
		}

	stop_timer(&(timers.prog_timer));

	// free device copies and gauge conf
	free_device_data(&GC, &geo, &param, NULL, &meas_aux, 1);
	free_gauge_conf(&GC, &param);

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);

	// free geometry
	free_geometry(&geo, &param);

	// print simulation details
	print_parameters_conf_analysis(&param, &timers);

	// free input parameters
	free_gparam(&param);
	}

void print_template_input(void)
	{
	char const *const template_name = "template_input.example";
	FILE *fp = fopen(template_name, "w");
	REQUIRE(fp != NULL, "failed to open %s", template_name);

	print_template_volume_parameters(fp);
	print_template_pt_parameters(fp);
	print_template_twist_parameters(fp);
	#ifdef MULTICANONICAL_MODE
	print_template_multicanonic_parameters(fp);
	#endif
	print_template_simul_parameters(fp);
	print_template_adaptive_gradflow_parameters(fp);
	print_template_cooling_parameters(fp);
	print_template_output_parameters(fp);

	fclose(fp);

	printf("Generated template input: %s\n\n", template_name);
	}

int main(int argc, char **argv)
	{
	if(argc != 4)
		{
		print_program_details("measurements on configurations saved by other Yang-Mills programs",
		              "determined by the input file of the Yang-Mills program that generated the configurations",
		              "no updates; gradient flow (fixed or adaptive step) and cooling",
		              "determined by the input file of the Yang-Mills program that generated the configurations",
		              argv[0], "input_file start_index stop_index");

		print_template_input();
		print_compilation_details();

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
	REQUIRE(stop_index > start_index, "stop_index must be larger than start_index");

	real_main(argv[1], start_index, stop_index);

	return EXIT_SUCCESS;
	}

#endif
