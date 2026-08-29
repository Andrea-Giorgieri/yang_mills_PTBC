#ifndef KP_TEST_C
#define KP_TEST_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/gauge_group.h"
#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/timing.h"

void real_main(char const *in_file)
	{
	GAUGE_GROUP link, staple, matrix, link_aux;
	GParam param;
	Time_Utils timers;

	// read input file
	readinput(in_file, &param);

	// initialize staple
	FILE *filep = fopen(param.d_conf_file, "r");
	read_from_file(filep, &staple);
	fclose(filep);

	// initialize link
	one(&link);

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// initialize random generator
	initrand(param.d_randseed);

	// init file for measurements
	FILE *datafilep = fopen(param.d_data_file, "w+");
	fprintf(datafilep, "# upd_index reE imE reEU imEU reEUb imEUb\n");

	stop_timer(&(timers.init_timer));

	// Monte Carlo begin
	for(long count = 0; count < param.d_sample; count++)
		{
		start_timer(&(timers.step_timer));

		// perform a heatbath update
		start_timer(&(timers.update_timer));
		single_heatbath(&link, &staple, &param);
		stop_timer(&(timers.update_timer));

		// perform measures

		start_timer(&(timers.meas_timer));

		// plaquette with link
		times(&matrix, &link, &staple);
		double const reE = retr(&matrix);
		double const imE = imtr(&matrix);

		// plaquette with unitarized link
		equal(&link_aux, &link);
		unitarize(&link_aux);
		times(&matrix, &link_aux, &staple);
		double const reEU = retr(&matrix);
		double const imEU = imtr(&matrix);

		// plaquette with badly-unitarized link
		equal(&link_aux, &link);
		#if NCOLOR > 2
		bad_unitarize_SuN(&link_aux, 1e10, stdout, 1);
		#else
		unitarize(&link_aux);
		#endif
		times(&matrix, &link_aux, &staple);
		double const reEUb = retr(&matrix);
		double const imEUb = imtr(&matrix);

		// unitarize link?
		//unitarize_SuN(&link_aux);

		stop_timer(&(timers.meas_timer));

		// write measures
		fprintf(datafilep, "%10ld % 22.16e % 22.16e % 22.16e % 22.16e % 22.16e % 22.16e \n", count, reE, imE, reEU, imEU, reEUb, imEUb);

		stop_timer(&(timers.step_timer));
		if(wall_time_check(&timers) == 1) break;
		}

	// Monte Carlo end
	stop_timer(&(timers.prog_timer));

	// print simulation details
	print_parameters_local_pt_agf(&param, &timers);

	// close data file
	fclose(datafilep);

	// free hierarchical update parameters
	free_hierarc_params(&param);
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
	if(argc != 2)
		{
		int parallel_tempering = 1;
		int twisted_bc = 1;
		print_authors(parallel_tempering, twisted_bc);

		printf("Usage: %s input_file\n\n", argv[0]);

		print_compilation_details();
		print_template_input();

		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");
	REQUIRE(STDIM == 4 && NCOLOR > 1, "PTBC not implemented for STDIM != 4 or NCOLOR < 2");

	real_main(argv[1]);

	return EXIT_SUCCESS;
	}

#endif
