#ifndef YM_AGF_C
#define YM_AGF_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/function_pointers.h"
#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"

void real_main(char *in_file)
	{
	Gauge_Conf GC, GC_old, help1, help2, help3;
	Geometry geo;
	GParam param;

	int meas_count, gradflowrepeat, accepted;
	double gftime, gftime_step;
	double 	*meanplaq, *clover_energy, *chi_prime, *charge, **charge_prime, *sum_q_timeslices;

	FILE *datafilep, *chiprimefilep, *topchar_tcorr_filep;
	time_t time1, time2;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
	// omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// read input file
	readinput(in_file, &param);

	// this code has to start from saved conf.
	param.d_start = 2;

	// initialize random generator
	initrand(param.d_randseed);

	// open data files
	init_data_file(&datafilep, &chiprimefilep, &topchar_tcorr_filep, &param);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);

	// initialize gauge configurations
	init_gauge_conf(&GC, &param);
	init_gauge_conf_from_gauge_conf(&GC_old, &GC, &param);
	init_gauge_conf_from_gauge_conf(&help1, &GC, &param);
	init_gauge_conf_from_gauge_conf(&help2, &GC, &param);
	init_gauge_conf_from_gauge_conf(&help3, &GC, &param);
	
	// allocate meas arrays
	gradflowrepeat = (int)floor((param.d_agf_length+MIN_VALUE)/param.d_agf_meas_each); // number of meas to perform
	allocate_measures_arrays(gradflowrepeat, &param, &meanplaq, &clover_energy, &charge, &sum_q_timeslices, &chi_prime, &charge_prime);
	
	// meas no gradflow
	perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tcorr_filep);

	// gradflow starts
	time(&time1);
	gftime = 0.0;						// gradient flow time
	gftime_step = param.d_agf_step; 	// initial integration time step
	meas_count = 0; 					// meas counter
	
	while(meas_count < gradflowrepeat)
		{
		// integration step
		gradflow_RKstep_adaptive(&GC, &GC_old, &help1, &help2, &help3, &geo, &param, &gftime, &gftime_step, &accepted);
		
		// step accepted, perform measures if it is time to do so
		if (accepted == 1 && fabs(gftime - param.d_agf_meas_each*(meas_count+1)) - param.d_agf_time_bin < MIN_VALUE )
			{
			perform_measures_aux(&GC, &geo, &param, meas_count, meanplaq, clover_energy, charge, sum_q_timeslices, chi_prime, charge_prime, topchar_tcorr_filep);
			meas_count = meas_count + 1;
			}
		
		// adapt step to the time of next measure
		if ((gftime + gftime_step - param.d_agf_meas_each*(meas_count+1)) > param.d_agf_time_bin )
			{
			gftime_step = param.d_agf_meas_each*(meas_count+1) - gftime;
			}
		}
	time(&time2);
	// gradflow ends

	// print meas gradflow, close files
	print_measures_arrays(gradflowrepeat, GC.update_index, &param, meanplaq, clover_energy, charge, chi_prime,charge_prime, datafilep, chiprimefilep);
	fprintf(datafilep, "\n");
	fclose(datafilep);
	if (param.d_chi_prime_meas == 1 ) fclose(chiprimefilep);
	if (param.d_topcharge_tcorr_meas == 1 ) fclose(topchar_tcorr_filep);

	// free memory
	free_measures_arrays(gradflowrepeat, &param, meanplaq, clover_energy, charge, sum_q_timeslices, chi_prime, charge_prime);
	free_gauge_conf(&GC, &param);
	free_gauge_conf(&GC_old, &param);
	free_gauge_conf(&help1, &param);
	free_gauge_conf(&help2, &param);
	free_gauge_conf(&help3, &param);
	
	// free geometry
	free_geometry(&geo, &param);

	// print simulation details
	print_parameters_agf(&param, time1, time2);
	}

void print_template_input(void)
	{
	FILE *fp;

	fp = fopen("template_input.example", "w");

	if(fp == NULL)
		{
		fprintf(stderr, "Error in opening the file template_input.example (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		print_template_volume_parameters(fp);
		print_template_adaptive_gradflow_parameters(fp);
		
		fprintf(fp, "# Observables to measure\n");
		fprintf(fp, "plaquette_meas        0  # 1=YES, 0=NO\n");
		fprintf(fp, "clover_energy_meas    1  # 1=YES, 0=NO\n");
		fprintf(fp, "charge_meas           1  # 1=YES, 0=NO\n");
		fprintf(fp, "polyakov_meas         0  # 1=YES, 0=NO\n");
		fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
		fprintf(fp, "topcharge_tcorr_meas  0  # 1=YES, 0=NO\n");
		fprintf(fp,"\n");
		
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

int main(int argc, char **argv)
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