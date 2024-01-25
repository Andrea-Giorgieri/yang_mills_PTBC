#ifndef YM_T0_C
#define YM_T0_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
#include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
	{
	Gauge_Conf GC, help1, help2;
	Geometry geo;
	GParam param;
	
	long count;
	const long max_count=10000;
	double gftime, energy_clover, energy_clover_old, tch, ris;
	
	FILE *datafilep;
	time_t time1, time2;
	
	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif
	
	// read input file
	readinput(in_file, &param);
	
	// this code has to start from saved conf.
	param.d_start=2;
	
	// initialize random generator
	initrand(param.d_randseed);
	
	// open data_file
	datafilep=fopen(param.d_data_file, "a");
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// initialize gauge configurations
	init_gauge_conf(&GC, &param);
	init_gauge_conf_from_gauge_conf(&help1, &GC, &param);
	init_gauge_conf_from_gauge_conf(&help2, &GC, &param);
	
	time(&time1);
	gftime=0.0;
	count=0;
	energy_clover_old=0.0;
	while(count<max_count)
		{
		gradflow_RKstep(&GC, &help1, &help2, &geo, &param, param.d_gfstep);
		gftime+=param.d_gfstep;
		
		clover_disc_energy(&GC, &geo, &param, &energy_clover);
		tch=topcharge(&GC, &geo, &param);
		
	 	fprintf(datafilep, "%.13lf	%.13lf	%.13lf	%.13lf\n", gftime, energy_clover, energy_clover*gftime*gftime, tch);
		if(energy_clover*gftime*gftime>0.3)
			{
			ris = gftime - param.d_gfstep + (0.3-energy_clover_old*gftime*gftime)*param.d_gfstep/(energy_clover*gftime*gftime-energy_clover_old*gftime*gftime);
			fprintf(datafilep, "%.13lf\n\n", ris);
			count=(max_count+10);
			}
		fflush(datafilep);
		
		count++;
		energy_clover_old=energy_clover;
		}
	time(&time2);
	
	if(count==max_count)
		{
		fprintf(stderr, "max_count reached in (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	// close data file
	fclose(datafilep);
	
	// print simulation details
	print_parameters_t0(&param, time1, time2);
	
	// free gauge configurations
	free_gauge_conf(&GC, &param);
	free_gauge_conf(&help1, &param);
	free_gauge_conf(&help2, &param);
	
	// free geometry
	free_geometry(&geo, &param);
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
		fprintf(fp, "size 12 4 4 12\n");
		fprintf(fp,"\n");
		fprintf(fp, "# For gradient flow evolution\n");
		fprintf(fp, "gfstep 0.01  # integration step for gradient flow\n");
		fprintf(fp, "\n");
		fprintf(fp, "# Output files\n");
		fprintf(fp, "conf_file  conf.dat\n");
		fprintf(fp, "data_file  dati.dat\n");
		fprintf(fp, "log_file   log.dat\n");
		fprintf(fp, "\n");
		fprintf(fp, "randseed 0 #(0=time)\n");
		fclose(fp);
		}
	}


int main (int argc, char **argv)
	{
	char in_file[50];

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

