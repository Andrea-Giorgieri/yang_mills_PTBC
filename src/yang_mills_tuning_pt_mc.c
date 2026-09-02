#ifndef YM_TUNING_PT_MC_C
#define YM_TUNING_PT_MC_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/timing.h"

void real_main(char const *in_file)
	{
	Gauge_Conf *GC;
	Geometry geo;
	GParam param;
	Rect_Utils rect_aux;
	Acc_Utils acc_counters;
	Meas_Utils *meas_aux;
	Tune_Utils tune_utils;
	Time_Utils timers;

	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	int count, tune_check;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// check if program was compiled in multicanonical mode
	#ifndef MULTICANONICAL_MODE
	REQUIRE(0, "this program can be used only in MULTICANONICAL_MODE");
	#endif

	// read input file
	readinput(in_file, &param);

	// initialize timers
	init_time_utils(&timers, param.d_walltime);
	start_timer(&(timers.prog_timer));
	start_timer(&(timers.init_timer));

	// initialize random generator
	initrand(param.d_randseed);

	// initialize geometry
	init_geometry(&geo, &param);

	// initialize arrays and flags for tuning of topo_potential
	init_tune_utils(&tune_utils, &param);

	// initialize rectangles for hierarchical update and swap
	init_rect_utils(&rect_aux, &param);

	// initialize swap acceptance and multicanonic Metropolis acceptance arrays, open multicanonic acceptance file
	init_acc_utils(&acc_counters, &param);

	// initialize auxiliary arrays and lattices for measurements, open data files
	init_meas_utils_replica(&meas_aux, &param);

	// initialize gauge configuration of all replica and volume defects
	init_gauge_conf_replica(&GC, &geo, &param);

	stop_timer(&(timers.init_timer));

	// if number of MC updates is set to 0, perform measures on initialized configuration
	if(param.d_sample == 0)
		{
		start_timer(&(timers.step_timer));
		start_timer(&(timers.meas_timer));
		perform_measures_localobs(&(GC[0]), &geo, &param, &(meas_aux[0]));
		#ifdef REPLICA_MEAS_MODE
		for(int i = 1; i < param.d_N_replica_pt; i++)
			perform_measures_localobs(&(GC[i]), &geo, &param, &(meas_aux[i]));
		#endif
		start_timer(&(timers.meas_timer));
		stop_timer(&(timers.step_timer));
		}

	// Monte Carlo begin
	for(count = 0; count < param.d_sample; count++)
		{
		start_timer(&(timers.step_timer));

		// perform a single step of parallel tempering wth hierarchical update and print state of replica swaps
		start_timer(&(timers.update_timer));
		parallel_tempering_with_hierarchical_update(GC, &geo, &param, &rect_aux, &acc_counters);
		stop_timer(&(timers.update_timer));
		print_conf_labels(GC, &param, &acc_counters);

		// perform measures only on physical replica, unless REPLICA_MEAS_MODE is enabled
		if(GC[0].update_index % param.d_measevery == 0 && GC[0].update_index >= param.d_thermal)
			{
			start_timer(&(timers.meas_timer));
			perform_measures_localobs(&(GC[0]), &geo, &param, &(meas_aux[0]));
			#ifdef REPLICA_MEAS_MODE
			for(int i = 1; i < param.d_N_replica_pt; i++)
				perform_measures_localobs(&(GC[i]), &geo, &param, &(meas_aux[i]));
			#endif
			stop_timer(&(timers.meas_timer));
			}

		// save configuration of all replica for backup
		if(param.d_saveconf_back_every != 0)
			{
			if(GC[0].update_index % param.d_saveconf_back_every == 0)
				{
				// simple
				write_replica_on_file(GC, &param);
				// backup copy
				write_replica_on_file_back(GC, &param);
				}
			}

		// save configuration of physical replica for offline analysis
		if(param.d_saveconf_analysis_every != 0)
			{
			if(GC[0].update_index % param.d_saveconf_analysis_every == 0)
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

		// tune topo_potential
		tune_topo_potential(GC, &param, &tune_utils);

		// save current potentials
		if(param.d_topo_tuning_save_every != 0)
			{
			if(GC[0].update_index % param.d_topo_tuning_save_every == 0)
				{
				strcpy(name, param.d_topo_potential_file);
				strcat(name, "_step_");
				sprintf(aux, "%ld", GC[0].update_index);
				strcat(name, aux);
				write_topo_potential(&param, name);
				}
			}

		// update tuning_stp and check if tuning is completed
		tune_check = update_tuning_stp(&tune_utils, &param);

		stop_timer(&(timers.step_timer));
		if(tune_check == 1) print_tuning_stp(GC[0].update_index, &tune_utils, &param);
		if(tune_check == -1) break;
		if(wall_time_check(&timers) == 1) break;
		}

	// Monte Carlo end
	stop_timer(&(timers.prog_timer));

	// save topo potential
	write_topo_potential(&param, param.d_topo_potential_file);

	// save configurations
	if(param.d_saveconf_back_every != 0) write_replica_on_file(GC, &param);

	// print simulation details
	print_parameters_tuning_pt_mc(&param, &timers, count);

	// print acceptances of parallel tempering
	print_acceptances(&acc_counters, &param);

	// free gauge configurations
	free_replica(GC, &param);

	// free geometry
	free_geometry(&geo, &param);

	// free tune_utils
	free_tune_utils(&tune_utils);

	// free rectangles for hierarchical update and swap
	free_rect_utils(&rect_aux, &param);

	// free swap acceptance and multicanonic Metropolis acceptance arrays, close multicanonic file
	free_acc_utils(&acc_counters, &param);

	// free auxiliary arrays and lattices for measurements, close data files
	free_meas_utils_replica(meas_aux, &param);

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
	print_template_multicanonic_parameters(fp);
	print_template_multicanonic_tuning_parameters(fp);
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
