#ifndef DEBUG_AGF_VS_DELTA_C
#define DEBUG_AGF_VS_DELTA_C

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
	Rectangle swap_rectangle;
	Rectangle *most_update, *clover_rectangle;
	Acc_Utils acc_counters;
	int L_R_swap=1;
	
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	char name0[STD_STRING_LENGTH], name1[STD_STRING_LENGTH], name2[STD_STRING_LENGTH], name3[STD_STRING_LENGTH];
	int count;
	FILE *datafilep0, *datafilep1, *datafilep2, *datafilep3, *chiprimefilep, *swaptrackfilep, *topchar_tcorr_filep;
	FILE *step_filep0, *step_filep1, *step_filep2, *step_filep3;
	time_t time_mc_start, time_mc_end, time1, time2, time3, time4, time5, time_agf0, time_agf1, time_agf2, time_agf3;
	double delta0, delta1, delta2, delta3;
	
	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif
	
	// read input file
	readinput(in_file, &param);
	delta0 = param.d_agf_delta;
	delta1 = delta0/10.0;
	delta2 = delta1/10.0;
	delta3 = delta2/10.0;
	delta3 = delta3/1.0;
	
	// initialize random generator
	initrand(param.d_randseed);
	
	// open data_file
	strcpy(aux, param.d_data_file);
	strcpy(name0, aux);
	strcpy(name1, aux);
	strcpy(name2, aux);
	strcpy(name3, aux);
	strcat(name0, "_delta0");
	strcat(name1, "_delta1");
	strcat(name2, "_delta2");
	strcat(name3, "_delta3");
	
	strcpy(param.d_data_file, name0);
	init_data_file(&datafilep0, &chiprimefilep, &topchar_tcorr_filep, &param);
	strcpy(param.d_data_file, name1);
	init_data_file(&datafilep1, &chiprimefilep, &topchar_tcorr_filep, &param);
	strcpy(param.d_data_file, name2);
	init_data_file(&datafilep2, &chiprimefilep, &topchar_tcorr_filep, &param);
	strcpy(param.d_data_file, name3);
	init_data_file(&datafilep3, &chiprimefilep, &topchar_tcorr_filep, &param);
	strcpy(param.d_data_file, aux);
	
	step_filep0 = fopen("step_file0.dat", "a");
	if (step_filep0 == NULL) step_filep0 = fopen("step_file0.dat", "w");
	step_filep1 = fopen("step_file1.dat", "a");
	if (step_filep1 == NULL) step_filep0 = fopen("step_file1.dat", "w");
	step_filep2 = fopen("step_file2.dat", "a");
	if (step_filep2 == NULL) step_filep0 = fopen("step_file2.dat", "w");
	step_filep3 = fopen("step_file3.dat", "a");
	if (step_filep3 == NULL) step_filep0 = fopen("step_file3.dat", "w");
	
	// open swap tracking file
	init_swap_track_file(&swaptrackfilep, &param);
	
	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	
	// initialize gauge configurations replica and volume defects
	init_gauge_conf_replica(&GC, &param);
	
	// initialize rectangles for hierarchical update
	init_rect_hierarc(&most_update, &clover_rectangle, &param);
	
	// initialize rectangle for swap probability evaluation (L_R_swap = 1)
	init_rect(&swap_rectangle, L_R_swap, &param);
	
	// init acceptances array
	init_swap_acc_arrays(&acc_counters, &param);
	
	// Monte Carlo begin
	time(&time_mc_start);
	time_agf0 = 0;
	time_agf1 = 0;
	time_agf2 = 0;
	time_agf3 = 0;
	for(count=0; count < param.d_sample; count++)
		{
		// perform a single step of parallel tempering wth hierarchical update and print state of replica swaps
		parallel_tempering_with_hierarchical_update(GC, &geo, &param, most_update, clover_rectangle, &swap_rectangle, &acc_counters);
		print_conf_labels(swaptrackfilep, GC, &param);
		
		// perform measures only on homogeneous configuration
		if(GC[0].update_index % param.d_measevery == 0 && GC[0].update_index >= param.d_thermal)
			{
			param.d_agf_delta = delta0;
			time(&time1);
			perform_measures_localobs_with_adaptive_gradflow_debug(&(GC[0]), &geo, &param, datafilep0, chiprimefilep, topchar_tcorr_filep, step_filep0);
			//param.d_agf_delta = delta1;
			time(&time2);
			perform_measures_localobs_with_adaptive_gradflow_debug2(&(GC[0]), &geo, &param, datafilep1, chiprimefilep, topchar_tcorr_filep, step_filep1);
			//param.d_agf_delta = delta2;
			time(&time3);
			//perform_measures_localobs_with_adaptive_gradflow_debug(&(GC[0]), &geo, &param, datafilep2, chiprimefilep, topchar_tcorr_filep, step_filep2);
			//param.d_agf_delta = delta3;
			time(&time4);
			//perform_measures_localobs_with_adaptive_gradflow_debug(&(GC[0]), &geo, &param, datafilep3, chiprimefilep, topchar_tcorr_filep, step_filep3);
			time(&time5);
			time_agf0 += time2 - time1;
			time_agf1 += time3 - time2;
			time_agf2 += time4 - time3;
			time_agf3 += time5 - time4;
			}
		param.d_agf_delta = delta0;
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
	
	time(&time_mc_end);
	// Monte Carlo end
	
	// close data file
	fclose(datafilep0);
	fclose(datafilep1);
	fclose(datafilep2);
	fclose(datafilep3);
	if(param.d_chi_prime_meas==1) fclose(chiprimefilep);
	if(param.d_topcharge_tcorr_meas==1) fclose(topchar_tcorr_filep);
	fclose(step_filep0);
	fclose(step_filep1);
	fclose(step_filep2);
	fclose(step_filep3);
	
	// close swap tracking file
	if (param.d_N_replica_pt > 1) fclose(swaptrackfilep);
	
	// save configurations
	if (param.d_saveconf_back_every!=0)
		{
		write_replica_on_file(GC, &param);
		}
	
	// print simulation details
	print_parameters_debug_agf_vs_delta(&param, time_mc_end-time_mc_start, time_agf0, time_agf1, time_agf2, time_agf3);
	
	// print acceptances of parallel tempering
	print_acceptances(&acc_counters, &param);
	
	// free gauge configurations
	free_replica(GC, &param);
	
	// free geometry
	free_geometry(&geo, &param);
	
	// free rectangles for hierarchical update
	free_rect_hierarc(most_update, clover_rectangle, &param);
	
	// free rectangle for swap probability evaluation
	free_rect(&swap_rectangle);
	
	// free acceptances array
	end_swap_acc_arrays(&acc_counters, &param);
	
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
