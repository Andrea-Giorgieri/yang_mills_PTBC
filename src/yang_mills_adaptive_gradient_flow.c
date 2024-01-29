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
	double *meanplaq, *clover_energy, *chi_prime, *charge, *sum_q_timeslices;

	FILE *datafilep, *chiprimefilep, *topchar_tprof_filep;
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
	init_data_file(&datafilep, &chiprimefilep, &topchar_tprof_filep, &param);

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
	gradflowrepeat = (int)floor((param.d_agf_length + MIN_VALUE) / param.d_agf_meas_each); // number of meas to perform
	if(param.d_plaquette_meas == 1)
		if(posix_memalign((void **)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t)gradflowrepeat * sizeof(double)) != 0)
			{
			fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	if(param.d_clover_energy_meas == 1)
		if(posix_memalign((void **)&clover_energy, (size_t)DOUBLE_ALIGN, (size_t)gradflowrepeat * sizeof(double)) != 0)
			{
			fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	if(param.d_charge_meas == 1)
		if(posix_memalign((void **)&charge, (size_t)DOUBLE_ALIGN, (size_t)gradflowrepeat * sizeof(double)) != 0)
			{
			fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	if(param.d_topcharge_tprof_meas == 1) // notice the size, measures are handled differently
		if(posix_memalign((void **)&sum_q_timeslices, (size_t)DOUBLE_ALIGN, (size_t)param.d_size[0] * sizeof(double)) != 0)
			{
			fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	if(param.d_chi_prime_meas == 1)
		if(posix_memalign((void **)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t)gradflowrepeat * sizeof(double)) != 0)
			{
			fprintf(stderr, "Problems allocating an array of doubles! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}

	// meas no gradflow
	perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);

	// gradflow starts
	time(&time1);
	gftime = 0.0;					// gradient flow time
	gftime_step = param.d_agf_step; // initial integration time step
	meas_count = 0;					// meas counter

	while(meas_count < gradflowrepeat)
		{
		// integration step
		gradflow_RKstep_adaptive(&GC, &GC_old, &help1, &help2, &help3, &geo, &param, &gftime, &gftime_step, &accepted);

		// step accepted, perform measures if it is time to do so
		if(accepted == 1 && fabs(gftime - param.d_agf_meas_each * (meas_count + 1)) - param.d_agf_time_bin < MIN_VALUE)
			{
			if(param.d_plaquette_meas == 1)
				{
				double plaqs, plaqt;
				plaquette(&GC, &geo, &param, &plaqs, &plaqt);
				#if(STDIM == 4)
				meanplaq[meas_count] = 0.5 * (plaqs + plaqt);
				#else
				meanplaq[meas_count] = plaqt;
				#endif
				}
			if(param.d_clover_energy_meas == 1) clover_disc_energy(&GC, &geo, &param, &clover_energy[meas_count]);
			if(param.d_charge_meas == 1) charge[meas_count] = topcharge(&GC, &geo, &param);
			if(param.d_topcharge_tprof_meas == 1) topcharge_timeslices(&GC, &geo, &param, sum_q_timeslices, meas_count + 1, topchar_tprof_filep);
			if(param.d_chi_prime_meas == 1) chi_prime[meas_count] = topo_chi_prime(&GC, &geo, &param);
			meas_count = meas_count + 1;
			}

		// adapt step to the time of next measure
		if((gftime + gftime_step - param.d_agf_meas_each * (meas_count + 1)) > param.d_agf_time_bin) { gftime_step = param.d_agf_meas_each * (meas_count + 1) - gftime; }
		}
	time(&time2);
	// gradflow ends

	// print meas gradflow, close files
	for(int i = 0; i < gradflowrepeat; i++)
		{
		if(param.d_plaquette_meas == 1) fprintf(datafilep, "%.12g ", meanplaq[i]);
		if(param.d_clover_energy_meas == 1) fprintf(datafilep, "%.12g ", clover_energy[i]);
		if(param.d_charge_meas == 1) fprintf(datafilep, "%.12g ", charge[i]);
		if(param.d_chi_prime_meas == 1) fprintf(chiprimefilep, "%ld %.12lg %.12lg\n", GC.update_index, (i + 1) * param.d_agf_meas_each, chi_prime[i]);
		}
	fprintf(datafilep, "\n");
	fclose(datafilep);
	if(param.d_chi_prime_meas == 1) fclose(chiprimefilep);
	if(param.d_topcharge_tprof_meas == 1) fclose(topchar_tprof_filep);

	// free memory
	if(param.d_plaquette_meas == 1) free(meanplaq);
	if(param.d_clover_energy_meas == 1) free(clover_energy);
	if(param.d_charge_meas == 1) free(charge);
	if(param.d_topcharge_tprof_meas == 1) free(sum_q_timeslices);
	if(param.d_chi_prime_meas == 1) free(chi_prime);
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

	fp = fopen("template_input.in", "w");

	if(fp == NULL)
		{
		fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		fprintf(fp, "size 12 4 4 12  # Nt Nx Ny Nz\n");
		fprintf(fp, "\n");

		fprintf(fp, "# For adaptive gradient flow evolution\n");
		fprintf(fp, "agf_length       10    # total integration time for adaptive gradient flow\n");
		fprintf(fp, "agf_step       0.01    # initial integration step for adaptive gradient flow\n");
		fprintf(fp, "agf_meas_each     1    # time interval between measures during adaptive gradient flow\n");
		fprintf(fp, "agf_delta     0.001    # error threshold on gauge links for adaptive gradient flow\n");
		fprintf(fp, "agf_time_bin  0.001    # error threshold on time of measures for adaptive gradient flow\n");
		fprintf(fp, "\n");

		fprintf(fp, "# Observables to measure\n");
		fprintf(fp, "plaquette_meas        0  # 1=YES, 0=NO\n");
		fprintf(fp, "clover_energy_meas    1  # 1=YES, 0=NO\n");
		fprintf(fp, "charge_meas           1  # 1=YES, 0=NO\n");
		fprintf(fp, "polyakov_meas         0  # 1=YES, 0=NO\n");
		fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
		fprintf(fp, "topcharge_tprof_meas  0  # 1=YES, 0=NO\n");
		fprintf(fp, "\n");
		
		fprintf(fp, "# Input files\n");
		fprintf(fp, "conf_file             conf.dat\n");
		fprintf(fp, "\n");

		fprintf(fp, "# Output files\n");
		fprintf(fp, "data_file             dati.dat\n");
		fprintf(fp, "chiprime_data_file    chi_prime_cool.dat\n");
		fprintf(fp, "topcharge_tprof_file  topo_tprof_cool.dat\n");
		fprintf(fp, "log_file              log.dat\n");
		fprintf(fp, "\n");

		fprintf(fp, "randseed 0    #(0=time)\n");
		fclose(fp);
		}
	}

int main(int argc, char **argv)
	{
	char in_file[500];

	if(argc != 2)
		{
		printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
		printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
		printf("Usage: %s input_file\n\n", argv[0]);

		printf("Compilation details:\n");
		printf("\tN_c (number of colors): %d\n", NCOLOR);
		printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
		printf("\tNum_levels (number of levels): %d\n", NLEVELS);
		printf("\n");
		printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
		printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

		#ifdef DEBUG
		printf("\n\tDEBUG mode\n");
		#endif

		#ifdef OPENMP_MODE
		printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
		#endif

		printf("\n");

		#ifdef __INTEL_COMPILER
		printf("\tcompiled with icc\n");
		#elif defined(__clang__)
		printf("\tcompiled with clang\n");
		#elif defined(__GNUC__)
		printf("\tcompiled with gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		#endif

		print_template_input();

		return EXIT_SUCCESS;
		}
	else
		{
		if(strlen(argv[1]) >= STD_STRING_LENGTH) { fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n"); }
		else
			{
			strcpy(in_file, argv[1]);
			}
		}

	real_main(in_file);

	return EXIT_SUCCESS;
	}

#endif
