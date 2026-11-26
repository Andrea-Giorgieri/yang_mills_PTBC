#ifndef YM_LOCAL_AGF_C
#define YM_LOCAL_AGF_C

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
	omp_set_max_active_levels(1);
	#endif

	// read input file
	readinput(in_file, &param);

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

	// initialize multicanonic Metropolis acceptance arrays, open multicanonic acceptance file
	init_acc_utils(&acc_counters, &param);

	// initialize meas utils
	init_meas_utils(&meas_aux, &param, 0);

	// copy the shared data and the configuration to the device (one configuration, no rectangles)
	init_device_utils(&param, 1);
	init_device_data(&GC, &geo, &param, NULL, &meas_aux, 1);

	stop_timer(&(timers.init_timer));

	// if number of MC updates is set to 0, perform measures on initialized configuration
	if(param.d_sample == 0)
		{
		start_timer(&(timers.step_timer));
		start_timer(&(timers.meas_timer));
		perform_measures_localobs(&GC, &geo, &param, &meas_aux);
		stop_timer(&(timers.meas_timer));
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

		// print mean multicanonic acceptance
		#ifdef MULTICANONICAL_MODE
		print_multicanonic_acceptance(&GC, &param, &acc_counters);
		#endif

		// perform measures
		if(GC.update_index % param.d_measevery == 0 && GC.update_index >= param.d_thermal)
			{
			start_timer(&(timers.meas_timer));
			perform_measures_localobs(&GC, &geo, &param, &meas_aux);
			stop_timer(&(timers.meas_timer));
			}

		// save configuration for backup
		if(param.d_saveconf_back_every != 0)
			{
			if(GC.update_index % param.d_saveconf_back_every == 0)
				{
				copy_gauge_conf_to_host(&GC, &param, GC_LATTICE | GC_Z);
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
				copy_gauge_conf_to_host(&GC, &param, GC_LATTICE | GC_Z);
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

	// configuration back to the host, free device copies
	copy_gauge_conf_to_host(&GC, &param, GC_LATTICE | GC_Z);
	free_device_data(&GC, &geo, &param, NULL, &meas_aux, 1);

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);

	// save configuration
	if(param.d_saveconf_back_every != 0)
		{
		write_conf_on_file(&GC, &param);
		}

	// print simulation details
	print_parameters_local_agf(&param, &timers);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free acceptance counters, close multicanonic acceptance file
	free_acc_utils(&acc_counters, &param);

	// free geometry
	free_geometry(&geo, &param);

	// free input parameters
	free_gparam(&param);
	}


void print_template_input(void)
	{
	char const *const template_name = "template_input.example";
	FILE *fp = fopen(template_name, "w");
	REQUIRE(fp != NULL, "failed to open %s", template_name);

	print_template_volume_parameters(fp);
	print_template_twist_parameters(fp);
	#ifdef MULTICANONICAL_MODE
	print_template_multicanonic_parameters(fp);
	#endif
	print_template_simul_parameters(fp);
	print_template_adaptive_gradflow_parameters(fp);
	print_template_gradflow_parameters(fp);
	print_template_cooling_parameters(fp);
	print_template_output_parameters(fp);

	fclose(fp);

	printf("Generated template input: %s\n\n", template_name);
	}


int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		print_program_details("simulation of U(1)/SU(N) Yang-Mills theory in any dimension",
		                      "Wilson action, periodic/twisted/open BC; options (if Q defined): imaginary theta, V(Q)",
		                      "heatbath, overrelaxation, Metropolis for V(Q)",
		                      "plaquette, energy, charge, Polyakov loops; gradient flow, cooling",
		                      argv[0], "input_file");

		print_template_input();
		print_compilation_details();

		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");

	real_main(argv[1]);

	return EXIT_SUCCESS;
	}

#endif
