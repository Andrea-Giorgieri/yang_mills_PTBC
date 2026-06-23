#ifndef YM_TRACEDEF_C
#define YM_TRACEDEF_C

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

void real_main(char *in_file)
	{
	Gauge_Conf GC;
	Geometry geo;
	GParam param;
	Meas_Utils meas_aux;
	Time_Utils timers;

	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	double acc, acc_local;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
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
	init_indexing_lexeo();
	init_geometry(&geo, &param);

	// initialize gauge configuration
	init_gauge_conf(&GC, &geo, &param);

	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);

	// acceptance of the metropolis update
	acc = 0.0;

	stop_timer(&(timers.init_timer));

	// Monte Carlo begin (count starts from 1 to avoid problems using %)
	for(int count = 1; count < param.d_sample + 1; count++)
		{
		start_timer(&(timers.step_timer));

		update_with_trace_def(&GC, &geo, &param, &acc_local);
		acc += acc_local;

		if(count % param.d_measevery == 0 && count >= param.d_thermal)
			{
			perform_measures_localobs(&GC, &geo, &param, &meas_aux);
			}

		// save configuration for backup
		if(param.d_saveconf_back_every != 0)
			{
			if(count % param.d_saveconf_back_every == 0)
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
			if(count % param.d_saveconf_analysis_every == 0)
				{
				strcpy(name, param.d_conf_file);
				sprintf(aux, "%ld", GC.update_index);
				strcat(name, aux);
				write_conf_on_file_with_name(&GC, &param, name);
				}
			}
		stop_timer(&(timers.step_timer));
		if(wall_time_check(&timers) == 1) break;
		}

	// Monte Carlo end
	stop_timer(&(timers.prog_timer));

	acc /= (double) param.d_sample;

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);

	// save configuration
	if(param.d_saveconf_back_every != 0)
		{
		write_conf_on_file(&GC, &param);
		}

	// print simulation details
	print_parameters_tracedef(&param, &timers, acc);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free geometry
	free_geometry(&geo, &param);
	}


void print_template_input(void)
	{
	FILE *fp = fopen("template_input.example", "w");
	REQUIRE(fp != NULL, "failed to open template_input.example");

	print_template_volume_parameters(fp);
	fprintf(fp, "htracedef  1.1\n");
	fprintf(fp, "\n");
	print_template_simul_parameters(fp);
	print_template_metro_parameters(fp);
	print_template_cooling_parameters(fp);
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
