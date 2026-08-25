#ifndef YM_POLYCORR_LONG_C
#define YM_POLYCORR_LONG_C

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

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);

	// initialize gauge configuration
	init_gauge_conf(&GC, &geo, &param);

	// initialize ml_polycorr arrays
	alloc_polycorr_stuff(&GC, &param);

	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);

	stop_timer(&(timers.init_timer));

	// Monte Carlo begin
	if(param.d_start != 2) // NEW SIMULATION
		{
		start_timer(&(timers.step_timer));

		for(int count = 0; count < param.d_measevery; count++)
			{
			update(&GC, &geo, &param, &acc_counters);
			}

		// save configuration
		write_conf_on_file(&GC, &param);
		// backup copy
		write_conf_on_file_back(&GC, &param);

		// save ml polycorr arrays
		write_polycorr_on_file(&GC, &param, 0);

		stop_timer(&(timers.step_timer));
		}
	else // CONTINUATION OF PREVIOUS SIMULATION
		{
		start_timer(&(timers.step_timer));

		int iteration;

		// read multilevel stuff
		read_polycorr_from_file(&GC, &param, &iteration);

		if(iteration < 0) // update the conf, no multilevel
			{
			for(int count = 0; count < param.d_measevery; count++)
				{
				update(&GC, &geo, &param, &acc_counters);
				}

			// save configuration
			write_conf_on_file(&GC, &param);
			// backup copy
			write_conf_on_file_back(&GC, &param);

			// save multilevel stuff
			write_polycorr_on_file(&GC, &param, 0);
			}
		else // iteration >=0, perform multilevel
			{
			multilevel_polycorr_long_zero(&GC, &geo, &param, iteration);
			iteration += 1;
			if(iteration == param.d_ml_level0_repeat)
				{
				// print the measure
				perform_measures_polycorr_long(&GC, &param, &meas_aux);

				iteration = -1; // next time the conf will be updated, no multilevel
				}

			// save multilevel stuff
			write_polycorr_on_file(&GC, &param, iteration);
			}
		stop_timer(&(timers.step_timer));
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
	print_parameters_polycorr_long(&param, &timers);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free ml_polycorr
	free_polycorr_stuff(&GC, &param);

	// free geometry
	free_geometry(&geo, &param);
	}

void print_template_input(void)
	{
	FILE *fp = fopen("template_input.example", "w");
	REQUIRE(fp != NULL, "failed to open template_input.example");

	print_template_volume_parameters(fp);
	print_template_simul_parameters(fp);
	print_template_multilevel_parameters(fp);
	fprintf(fp, "ml_level0_repeat 1 # number of times level0 is repeated in long sim.\n");
	fprintf(fp, "dist_poly        2 # distance between the polyakov loop\n");
	fprintf(fp, "\n");
	print_template_output_parameters(fp);
	fclose(fp);
	}


int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		int parallel_tempering = 0;
		int twisted_bc = 0;
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
