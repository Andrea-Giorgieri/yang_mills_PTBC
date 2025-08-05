#ifndef YM_GF_C
#define YM_GF_C

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

	int count;
	double gftime, chi_prime, tch;

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

	// this code has to start from saved conf.
	param.d_start=2;

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);

	// initialize gauge configurations
	init_gauge_conf(&GC, &geo, &param);

	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);

	stop_timer(&(timers.init_timer));

	gftime=0.0;
	// count starts from 1 to avoid problems with %
	for(count=1; count < (param.d_ngfsteps+1); count++)
		{
		start_timer(&(timers.step_timer));

		gradflow_RKstep(&GC, &geo, &param, param.d_gfstep, &meas_aux);
		gftime+=param.d_gfstep;

		if ( (count % param.d_gf_meas_each) == 0)
			{
			tch=topcharge(&GC, &geo, &param);
			chi_prime=topo_chi_prime(&GC, &geo, &param);
			fprintf(meas_aux.datafilep, "%ld	%.16lf	%.16lf	%16lf\n", GC.update_index, gftime, tch, chi_prime);
			fflush(meas_aux.datafilep);
			}

		stop_timer(&(timers.step_timer));
		if (wall_time_check(&timers) == 1) break;
		}

	stop_timer(&(timers.prog_timer));

	// print simulation details
	print_parameters_gf(&param, &timers);

	// free gauge configurations
	free_gauge_conf(&GC, &param);

	// free geometry
	free_geometry(&geo, &param);

	// free meas utils
	free_meas_utils(meas_aux, &param, 0);
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
		print_template_gradflow_parameters(fp);
		fprintf(fp, "# Output files\n");
		fprintf(fp, "conf_file  conf.dat\n");
		fprintf(fp, "twist_file twist.dat\n");
		fprintf(fp, "data_file  dati.dat\n");
		fprintf(fp, "log_file   log.dat\n");
		fprintf(fp, "\n");
		fprintf(fp, "randseed 0 #(0=time)\n");
		fclose(fp);
		}
	}


int main (int argc, char **argv)
	{
	char in_file[500];

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
	else
		{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
			{
			fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
			}
		else
			{
			strcpy(in_file, argv[1]);
			}
		}

	real_main(in_file);

	return EXIT_SUCCESS;
	}

#endif

