#ifndef YM_LOCAL_AGF_C
#define YM_LOCAL_AGF_C

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
	FILE *datafilep, *chiprimefilep, *swaptrackfilep, *topchar_tprof_filep;
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
	init_data_file(&datafilep, &chiprimefilep, &topchar_tprof_filep, &param);
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// initialize gauge configurations replica and volume defects
	init_gauge_conf(&GC, &param);
	
	// Monte Carlo begin
	time(&time1);
	if(param.d_sample == 0) // no update is done, only measures are performed on read configuration
		{
		perform_measures_localobs_with_adaptive_gradflow(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
		}
	for(count=0; count < param.d_sample; count++)
		{
		// update conf
		update(&GC, &geo, &param);

		// perform measures
		if(GC.update_index % param.d_measevery == 0 && GC.update_index >= param.d_thermal)
			{
			perform_measures_localobs_with_adaptive_gradflow(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
			}
		
		// save configuration for backup
		if(param.d_saveconf_back_every != 0)
			{
			if(GC.update_index % param.d_saveconf_back_every == 0)
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
			if(GC.update_index % param.d_saveconf_analysis_every == 0)
				{
				strcpy(name, param.d_conf_file);
				strcat(name, "_step_");
				sprintf(aux, "%ld", GC.update_index);
				strcat(name, aux);
				write_conf_on_file_with_name(&GC, &param, name);
				}
			}
		}
	
	time(&time2);
	// Monte Carlo end
	
	// close data file
	fclose(datafilep);
	if (param.d_chi_prime_meas==1) fclose(chiprimefilep);
	if (param.d_topcharge_tprof_meas==1) fclose(topchar_tprof_filep);
	
	// save configuration
	if (param.d_saveconf_back_every!=0)
		{
		write_conf_on_file(&GC, &param);
		}

	// print simulation details
	print_parameters_local_agf(&param, time1, time2);
	
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
		int parallel_tempering = 0;
		print_authors(parallel_tempering);
		
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
			strcpy(in_file, argv[1]);
			real_main(in_file);
			return EXIT_SUCCESS;
			}
		}
	}

#endif
