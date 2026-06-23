#ifndef YM_TUBE_CONN_C
#define YM_TUBE_CONN_C

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

	// allocate ml_polycorr, ml_polyplaq and ml_polyplaqconn arrays
	alloc_tube_conn_stuff(&GC, &param);

	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);

	stop_timer(&(timers.init_timer));

	// Monte Carlo begin (count starts from 1 to avoid problems using %)
	for(int count = 1; count < param.d_sample + 1; count++)
		{
		start_timer(&(timers.step_timer));

		update(&GC, &geo, &param, &acc_counters);

		if(count % param.d_measevery == 0 && count >= param.d_thermal)
			{
			perform_measures_tube_conn(&GC, &geo, &param, &meas_aux);
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
		stop_timer(&(timers.step_timer));
		if(wall_time_check(&timers) == 1) break;
		}

	//Monte Carlo end
	stop_timer(&(timers.prog_timer));

	// init meas utils
	free_meas_utils(meas_aux, &param, 0);

	// save configuration
	if(param.d_saveconf_back_every != 0)
		{
		write_conf_on_file(&GC, &param);
		}

	// print simulation details
	print_parameters_tube_conn(&param, &timers);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free ml_polycorr and ml_polyplaq
	free_tube_conn_stuff(&GC, &param);

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
	fprintf(fp, "dist_poly      2 # distance between the polyakov loop\n");
	fprintf(fp, "transv_dist    2 # transverse distance from the polyakov correlator\n");
	fprintf(fp, "plaq_dir     1 0 # plaquette orientation for flux tube\n");
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

