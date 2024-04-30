#ifndef YM_LOCAL_PT_GF_C
#define YM_LOCAL_PT_GF_C

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
	Gauge_Conf *GC;
	Geometry geo;
	GParam param;
	Rect_Utils rect_aux;
	Acc_Utils acc_counters;
	Meas_Utils meas_aux_agf, meas_aux_gf;
	
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH], name_gf[STD_STRING_LENGTH], name_agf[STD_STRING_LENGTH];
	int count;
	FILE *swaptrackfilep;
	time_t time1, time2, time3, time4, time5, time6, gf_time, agf_time, dagf_time;
	
	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif
	
	// read input file
	readinput(in_file, &param);
	
	// initialize random generator
	initrand(param.d_randseed);

	// open swap tracking file
	init_swap_track_file(&swaptrackfilep, &param);
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// initialize gauge configurations replica and volume defects
	init_gauge_conf_replica(&GC, &geo, &param);
	
	// initialize rectangles for hierarchical update and swap
	init_rect_utils(&rect_aux, &param);
	
	// init acceptances array
	init_acc_utils(&acc_counters, &param);

	// init meas utils
	strcpy(aux, param.d_data_file);
	strcpy(name_gf, aux);
	strcat(name_gf, "_gf");
	strcpy(name_agf, aux);
	strcat(name_agf, "_agf");
	strcpy(param.d_data_file, name_gf);
	init_meas_utils(&meas_aux_agf, &param, 0);
	strcpy(param.d_data_file, name_agf);
	init_meas_utils(&meas_aux_gf, &param, 0);
	strcpy(param.d_data_file, aux);
	//step_filep = fopen("./data/step_file.dat", "w");
	
	// Monte Carlo begin
	time(&time1);
	gf_time = 0;
	agf_time = 0;
	dagf_time = 0;
	for(count=0; count < param.d_sample; count++)
		{
		// perform a single step of parallel tempering wth hierarchical update and print state of replica swaps
		parallel_tempering_with_hierarchical_update(GC, &geo, &param, &rect_aux, &acc_counters);
		print_conf_labels(swaptrackfilep, GC, &param);

		// perform measures only on homogeneous configuration
		if(GC[0].update_index % param.d_measevery == 0 && GC[0].update_index >= param.d_thermal)
			{
			time(&time3);
			perform_measures_localobs_with_adaptive_gradflow(&(GC[0]), &geo, &param, &meas_aux_agf);
			time(&time4);
			//perform_measures_localobs_with_adaptive_gradflow_debug(&(GC[0]), &geo, &param, datafilep, chiprimefilep, topchar_tcorr_filep, step_filep);
			time(&time5);
			perform_measures_localobs_with_gradflow(&(GC[0]), &geo, &param, &meas_aux_gf);
			time(&time6);
			agf_time += time4 - time3;
			dagf_time += time5 - time4;
			gf_time += time6 - time5;
			}

		// save configurations for backup
		if(param.d_saveconf_back_every!=0)
			{
			if(GC[0].update_index % param.d_saveconf_back_every == 0 )
				{
				// simple
				write_replica_on_file(GC, &param);
				// backup copy
				write_replica_on_file_back(GC, &param);
				}
			}

		// save homogeneous configuration for offline analysis
		if(param.d_saveconf_analysis_every!=0)
			{
			if(GC[0].update_index % param.d_saveconf_analysis_every == 0 )
				{
				strcpy(name, param.d_conf_file);
				strcat(name, "_step_");
				sprintf(aux, "%ld", GC[0].update_index);
				strcat(name, aux);
				write_conf_on_file_with_name(&(GC[0]), &param, name);
				
				strcpy(name, param.d_twist_file);
				strcat(name, "_step_");
				strcat(name, aux);
				write_twist_on_file_with_name(&(GC[0]), &param, name);
				}
			}
		}
	
	time(&time2);
	// Monte Carlo end
	
	// close data file
	free_meas_utils(meas_aux_agf, &param, 0);
	free_meas_utils(meas_aux_gf, &param, 0);
	//fclose(step_filep);
	
	// close swap tracking file
	if (param.d_N_replica_pt > 1) fclose(swaptrackfilep);
	
	// save configurations
	if (param.d_saveconf_back_every!=0)
		{
		write_replica_on_file(GC, &param);
		}
	
	// print simulation details
	print_parameters_debug_agf_vs_gf(&param, time1, time2, agf_time, dagf_time, gf_time);
	
	// print acceptances of parallel tempering
	print_acceptances(&acc_counters, &param);
	
	// free gauge configurations
	free_replica(GC, &param);
	
	// free geometry
	free_geometry(&geo, &param);
	
	// free rectangles for hierarchical update and swap
	free_rect_utils(&rect_aux, &param);
	
	// free acceptances array
	free_acc_utils(&acc_counters, &param);
	
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
		print_template_simul_parameters(fp);
		print_template_adaptive_gradflow_parameters(fp);
		print_template_gradflow_parameters(fp);
		print_template_output_parameters(fp);
		fclose(fp);
		}
	}

int main (int argc, char **argv)
	{
	char in_file[STD_STRING_LENGTH];
	
	if(argc != 2)
		{
		printf("\nSU(N) Hasenbusch Parallel Tempering implemented by Claudio Bonanno (claudiobonanno93@gmail.com) within yang-mills package\n");
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
