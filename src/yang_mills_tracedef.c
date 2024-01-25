#ifndef YM_TRACEDEF_C
#define YM_TRACEDEF_C

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
	Gauge_Conf GC;
	Geometry geo;
	GParam param;
	
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	int count;
	double acc, acc_local;
	FILE *datafilep;
	time_t time1, time2;
	
	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif
	
	// read input file
	readinput(in_file, &param);
	
	// initialize random generator
	initrand(param.d_randseed);
	
	// open data_file
	init_data_file(&datafilep, &datafilep, &datafilep, &param);
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// initialize gauge configuration
	init_gauge_conf(&GC, &param);
	
	// acceptance of the metropolis update
	acc=0.0;
	
	// montecarlo
	time(&time1);
	// count starts from 1 to avoid problems using %
	for(count=1; count < param.d_sample + 1; count++)
		{
		update_with_trace_def(&GC, &geo, &param, &acc_local);
		acc+=acc_local;
		
		if(count % param.d_measevery ==0 && count >= param.d_thermal)
			{
			perform_measures_localobs_with_tracedef(&GC, &geo, &param, datafilep);
			}
		
		// save configuration for backup
		if(param.d_saveconf_back_every!=0)
			{
			if(count % param.d_saveconf_back_every == 0 )
				{
				// simple
				write_conf_on_file(&GC, &param);
				
				// backup copy
				write_conf_on_file_back(&GC, &param);
				}
			}
		
		// save configuration for offline analysis
		if(param.d_saveconf_analysis_every!=0)
			{
			if(count % param.d_saveconf_analysis_every == 0 )
				{
				strcpy(name, param.d_conf_file);
				sprintf(aux, "%ld", GC.update_index);
				strcat(name, aux);
				write_conf_on_file_with_name(&GC, &param, name);
				}
			}
		}
	time(&time2);
	// montecarlo end
	
	acc/=(double)param.d_sample;
	
	// close data file
	fclose(datafilep);
	
	// save configuration
	if(param.d_saveconf_back_every!=0)
		{
		write_conf_on_file(&GC, &param);
		}
	
	// print simulation details
	print_parameters_tracedef(&param, time1, time2, acc);
	
	// free gauge configuration
	free_gauge_conf(&GC, &param);
	
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
		print_template_volume_parameters(fp);
		fprintf(fp, "htracedef  1.1\n");
		fprintf(fp,"\n");
		print_template_simul_parameters(fp);
		print_template_metro_parameters(fp);
		print_template_cooling_parameters(fp);
		print_template_output_parameters(fp);
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
		else strcpy(in_file, argv[1]);
		}
	
	real_main(in_file);
	
	return EXIT_SUCCESS;
	}

#endif
