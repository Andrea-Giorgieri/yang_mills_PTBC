#ifndef YM_LOCAL_C
#define YM_LOCAL_C

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
	Acc_Utils acc_counters;
	Meas_Utils meas_aux;
	
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	int count;
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
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// if measure-only mode is active conf must be read from file => start=2 ignoring the value found in input file
	if (param.d_sample == 0)
		{
		fprintf(stdout, "MEASURE-ONLY MODE: performing measures on configuration read from file %s, no update will be performed\n", param.d_conf_file);
		param.d_start=2;
		}
	
	// initialize gauge configuration
	init_gauge_conf(&GC, &geo, &param);
	
	// init meas utils
	init_meas_utils(&meas_aux, &param, 0);
	
	// --- TO BE REMOVED ---
/*	double* ratio;
	int errs=posix_memalign((void**)&ratio, (size_t)DOUBLE_ALIGN, (size_t) param.d_coolrepeat * sizeof(double));
 	if(errs!=0)
	{
		fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	check_correlation_decay_cooling(&GC, &geo, &param, ratio);
	FILE* ratio_fp;
	ratio_fp=fopen("CHECK_charge_cool.dat", "w");
	if (ratio_fp == NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", "CHECK_charge_cool.dat", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for (int iter=0; iter<param.d_coolrepeat; iter++) fprintf(ratio_fp, "%d  % 18.12e\n", (iter+1)*param.d_coolsteps, ratio[iter]);
	free(ratio);
	fclose(ratio_fp);
*/
	//---------------------------
	
	// Monte Carlo begins
	time(&time1);
	if (param.d_sample == 0) // no update is done, only measures are performed on read configuration
		{
		perform_measures_localobs(&GC, &geo, &param, &meas_aux);
		}
	else
		{
		for(count=0; count < param.d_sample; count++)
			{
			// update conf
			update(&GC, &geo, &param, &acc_counters);
			
			// measure local observables
			if(GC.update_index % param.d_measevery == 0 && GC.update_index >= param.d_thermal)
				{
				perform_measures_localobs(&GC, &geo, &param, &meas_aux);
				}
			
			// save configuration for backup
			if(param.d_saveconf_back_every!=0)
				{
				if(GC.update_index % param.d_saveconf_back_every == 0 )
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
				if(GC.update_index % param.d_saveconf_analysis_every == 0 )
					{
					strcpy(name, param.d_conf_file);
					strcat(name, "_step_");
					sprintf(aux, "%ld", GC.update_index);
					strcat(name, aux);
					write_conf_on_file_with_name(&GC, &param, name);
					}
				}
			}
		}
	time(&time2);
	// Monte Carlo ends
	
	// free meas utils
	free_meas_utils(meas_aux, &param, 0);
	
	// save configuration
	if(param.d_saveconf_back_every!=0)
		{
		write_conf_on_file(&GC, &param);
		}
	
	// print simulation details
	print_parameters_local(&param, time1, time2);
	
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
		print_template_simul_parameters(fp);
		print_template_cooling_parameters(fp);
		print_template_output_parameters(fp);
		fclose(fp);
		}
	}

int main (int argc, char **argv)
	{
	char in_file[STD_STRING_LENGTH];
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
