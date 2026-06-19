#ifndef KP_TEST_C
#define KP_TEST_C

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
	GAUGE_GROUP link, staple, matrix, link_aux;
	GParam param;
	Time_Utils timers;

	long count;
	double reE, imE, reEU, imEU, reEUb, imEUb;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// read input file
	readinput(in_file, &param);
	
	// initialize staple
	FILE *filep = fopen(param.d_conf_file, "r");
	read_from_file_SuN(filep, &staple);
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
	for(count=0; count < param.d_sample; count++)
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
		reE = retr(&matrix);
		imE = imtr(&matrix);
		
		// plaquette with unitarized link
		equal(&link_aux, &link);
		unitarize_SuN(&link_aux);
		times(&matrix, &link_aux, &staple);
		reEU = retr(&matrix);
		imEU = imtr(&matrix);
		
		// plaquette with badly-unitarized link
		equal(&link_aux, &link);
		bad_unitarize_SuN(&link_aux, 1e10, stdout, 1);
		times(&matrix, &link_aux, &staple);
		reEUb = retr(&matrix);
		imEUb = imtr(&matrix);
		
		// unitarize link?
		//unitarize_SuN(&link_aux);
		
		stop_timer(&(timers.meas_timer));
		
		// write measures
		fprintf(datafilep, "%10ld % 22.16e % 22.16e % 22.16e % 22.16e % 22.16e % 22.16e \n", count, reE, imE, reEU, imEU, reEUb, imEUb);
		
		stop_timer(&(timers.step_timer));
		if (wall_time_check(&timers) == 1) break;
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
	FILE *fp;

	fp=fopen("template_input.example", "w");

	if(fp==NULL)
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

int main (int argc, char **argv)
	{
	char in_file[STD_STRING_LENGTH];

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
	else
		{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
			{
			fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in /include/macro.h\n");
			return EXIT_SUCCESS;
			}
		else
			{
			#if(STDIM==4 && NCOLOR>1)
				strcpy(in_file, argv[1]);
				real_main(in_file);
				return EXIT_SUCCESS;
			#else
				fprintf(stderr, "Parallel tempering of volume defect not implemented for STDIM =/= 4 and N_color < 2.\n");
				return EXIT_SUCCESS;
			#endif
			}
		}
	}

#endif
