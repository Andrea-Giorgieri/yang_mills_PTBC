#ifndef DEBUG_TWIST_C
#define DEBUG_TWIST_C

#include "../include/macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/gauge_group.h"
#include "../include/gauge_conf.h"
#include "../include/gparam.h"
#include "../include/random.h"
#include "../include/sun.h"

void print_complex(double complex z);

void print_twist_factor(int const i, int const j, long const r, double complex const z);

void print_plaquette(int const i, int const j, long const r, double complex const z);

void conf_translation_dir(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int dir); //dir not random, debug only

double action_with_defect(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int a, int b, long r, int i, int j);


void real_main(char const *in_file)
	{
	long const r_fixed = 0;
	double complex z, Z_a, Z_b;

	Gauge_Conf *GC;
	GAUGE_GROUP M, N;
	Geometry geo;
	GParam param;
	Rectangle swap_rectangle;
	Rectangle *most_update, *clover_rectangle;
	Acc_Utils acc_counters;
	int L_R_swap = 1;

	readinput(in_file, &param);
	if(param.d_N_replica_pt < 2) param.d_N_replica_pt = 2;
	param.d_plaquette_meas = 1;
	param.d_clover_energy_meas = 1;
	param.d_charge_meas = 1;
	param.d_chi_prime_meas = 1;
	param.d_polyakov_meas = 1;
	param.d_energy_slices_meas = 1;
	param.d_charge_slices_meas = 1;

	initrand(param.d_randseed);
	init_geometry(&geo, &param);
	init_gauge_conf_replica(&GC, &geo, &param);
	init_rect_hierarc(&most_update, &clover_rectangle, &param);
	init_rect(&swap_rectangle, L_R_swap, &param);
	init_acc_utils(&acc_counters, &param);

	printf("\n***********************************\n");
	printf("PROGRAM FOR THE DEBUGGING OF TWIST\n");
	printf("***********************************\n\n");

	printf("N = %s", QUOTEME(NCOLOR));

	printf("\n\n");
	printf("PRINT TWIST FACTORS AT ORIGIN ...\n\n");

	for(int i = 0; i < STDIM; i++)
		{
		for(int j = i + 1; j < STDIM; j++)
			{
			z = GC[0].Z[0][dirs_to_si(i, j)];
			print_twist_factor(i, j, r_fixed, z);
			z = GC[0].Z[0][dirs_to_si(j, i)];
			print_twist_factor(j, i, r_fixed, z);
			}
		}

	printf("\n\n");
	printf("VERIFY THAT TWISTED COLD START GIVES ZERO ENERGY ...\n\n");

	param.d_start = 3;
	int err = 0;
	init_gauge_conf_replica(&GC, &geo, &param);
	for(long r = 0; r < param.d_volume; r++)
		{
		double energy = 0.0;
		for(int i = 0; i < STDIM; i++)
			for(int j = i + 1; j < STDIM; j++)
				{
				clover(&(GC[0]), &geo, &param, r, i, j, &M);
				ta(&M);
				equal(&N, &M);
				times_equal(&M, &N);
				energy += -NCOLOR * retr(&M) / 16.0;
				}
		if(fabs(energy) > MIN_VALUE)
			{
			err = 1;
			printf("    ERROR: E_cl[%ld] = %g\n\n", r, energy);
			}
		}
	if(err == 0) printf("Vacuum Clover Energy OK\n");
	err = 0;

	for(long r = 0; r < (param.d_volume); r++)
		{
		for(int i = 0; i < STDIM; i++)
			for(int j = i + 1; j < STDIM; j++)
				{
				double complex trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r, i, j);
				if(cabs(trace_plaquettep - 1.0) > MIN_VALUE)
					{
					err = 1;
					printf("    ERROR: ");
					print_plaquette(i, j, r, z);
					printf("\n");
					}
				}
		}
	if(err == 0) printf("Vacuum Plaquette OK\n");
	free_replica(GC, &param);
	printf("\n\n");
	printf("VERIFY THAT plaquettep GETS CONJUGATED SWAPPING THE DIRS OF A PLAQUETTE AT ORIGIN...\n\n");

	param.d_start = 1;
	init_gauge_conf_replica(&GC, &geo, &param);
	for(int i = 0; i < STDIM; i++)
		{
		for(int j = i + 1; j < STDIM; j++)
			{
			double complex trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r_fixed, i, j);
			double complex trace_plaquettep_swap = plaquettep_complex(&(GC[0]), &geo, &param, r_fixed, j, i);

			if(cabs(trace_plaquettep_swap - conj(trace_plaquettep)) < MIN_VALUE)
				{
				printf("Trace on plane (%d, %d) OK\n", i, j);
				}
			else
				{
				printf("    ERROR:\n    ");
				print_plaquette(i, j, r_fixed, trace_plaquettep);
				printf("    ");
				print_plaquette(j, i, r_fixed, trace_plaquettep_swap);
				printf("\n");
				}
			}
		}

	printf("\n\n");
	printf("VERIFY THAT calcstaples AND plaquettep GIVE THE SAME RESULT FOR THE TRACES OF A LINK AT ORIGIN ...\n\n");

	for(int i = 0; i < STDIM; i++)
		{
		calcstaples_wilson(&(GC[0]), &geo, &param, r_fixed, i, &M);
		times_equal(&M, &(GC[0].lattice[r_fixed][i]));
		double complex trace_calcstaples = retr(&M) + I * imtr(&M);

		double complex trace_plaquettep = 0.0 + I * 0.0;
		for(int j = 0; j < STDIM; j++)
			{
			if(j != i)
				{
				trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, r_fixed, i, j);
				trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r_fixed, j), j, i);
				}
			}

		if(cabs(trace_calcstaples - trace_plaquettep) < MIN_VALUE)
			{
			printf("Trace dir %d OK\n", i);
			}
		else
			{
			printf("    ERROR:\n    ");
			printf("	dir %d calcstaples = ", i);
			print_complex(trace_calcstaples);
			printf(", plaquettep  = ");
			print_complex(trace_plaquettep);
			printf("\n");
			}
		}

	printf("\n\n");
	printf("VERIFY THAT clover AND plaquettep GIVE THE SAME RESULT FOR THE TRACE OF A CLOVER AT ORIGIN ...\n\n");

	for(int i = 0; i < STDIM; i++)
		{
		for(int j = i + 1; j < STDIM; j++)
			{
			clover(&(GC[0]), &geo, &param, r_fixed, i, j, &M);
			double complex trace_clover = retr(&M) + I * imtr(&M);

			double complex trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r_fixed, i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r_fixed, j), i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r_fixed, i), i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, nnm(&geo, r_fixed, i), j), i, j);

			if(cabs(trace_clover - trace_plaquettep) < MIN_VALUE)
				{
				printf("Trace plane (%d, %d) OK\n", i, j);
				}
			else
				{
				printf("    ERROR:\n    ");
				printf("plane (%d, %d) clover = ", i, j);
				print_complex(trace_clover);
				printf(", plaquettep = ");
				print_complex(trace_plaquettep);
				printf("\n\n");
				}
			}
		}

	printf("\n\n");
	printf("VERIFY THAT TRANSLATIONS DON'T CHANGE THE MEAN PLAQUETTE ...\n\n");

	for(int i = 0; i < STDIM; i++)
		{
		double plaqs, plaqt, plaqs_new, plaqt_new;
		plaquette(&(GC[0]), &geo, &param, &plaqs, &plaqt);         //mean plaquette
		conf_translation_dir(&(GC[0]), &geo, &param, i);           // translation
		plaquette(&(GC[0]), &geo, &param, &plaqs_new, &plaqt_new); //mean plaquette new

		if(fabs(plaqs_new - plaqs) < MIN_VALUE)
			{
			printf("Space plaquette dir %d OK\n", i);
			}
		else
			{
			printf("    ERROR:\n    ");
			printf("dir %d DeltaPlaqs = %g, Plaqs = %g\n\n", i, plaqs_new - plaqs, plaqs);
			}

		if(fabs(plaqt_new - plaqt) < MIN_VALUE)
			{
			printf("Time plaquette dir %d OK\n", i);
			}
		else
			{
			printf("    ERROR:\n    ");
			printf("dir %d DeltaPlaqt = %g, Plaqt = %g\n\n", i, plaqt_new - plaqt, plaqt);
			}
		}

	printf("\n\n");
	printf("VERIFY THAT metropolis_single_swap SWAPS THE TWIST FACTORS ...\n\n");

	err = 0;
	for(long r = 0; r < param.d_volume; r++)
		{
		for(int j = 0; j < param.d_n_planes; j++)
			{
			Z_a = GC[0].Z[r][j];
			Z_b = GC[1].Z[r][j];
			metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);

			if(cabs(GC[1].Z[r][j] - Z_a) > MIN_VALUE || cabs(GC[0].Z[r][j] - Z_b) > MIN_VALUE)
				{
				err = 1;
				printf("    ERROR:\n    ");
				printf("Error swapping Z[%ld][%d]\n\n", r, j);
				}
			}
		}
	if(err == 0) printf("Z[r][j] swaps OK\n");
	metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);

	printf("\n\n");
	printf("VERIFY THAT delta_action_swap IS CORRECT ...\n\n");

	err = 0;
	param.d_start = 1;
	init_gauge_conf_replica(&GC, &geo, &param);
	for(int i = 0; i < STDIM; i++)
		for(int j = 0; j < STDIM; j++)
			if (i != j)
				for(long r = 0; r < param.d_volume; r++)
					{
					double delta_action = delta_action_swap(GC, &geo, &param, r, i, j, 0, 1);
					double action = action_with_defect(GC, &geo, &param, 0, 1, r, i, j);
					metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);
					double swapped_action = action_with_defect(GC, &geo, &param, 0, 1, r, i, j);

					if(fabs(delta_action - (swapped_action - action)) > MIN_VALUE)
						{
						err = 1;
						printf("    ERROR:\n    ");
						printf("site %ld, plane (%d,%d):\n    ", r, i, j);
						printf("Delta S manual: %lf, function: %lf\n", swapped_action - action, delta_action);
						}
					}
	if(err == 0) printf("delta_action_swap OK\n");

	free_replica(GC, &param);
	free_geometry(&geo, &param);
	free_rect_hierarc(most_update, clover_rectangle, &param);
	free_rect(&swap_rectangle);
	free_acc_utils(&acc_counters, &param);
	free_hierarc_params(&param);

	printf("\nALL TESTS ENDED\n\n");
	}

void conf_translation_dir(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int const dir)
	{
	// copy the conf in an auxiliary one (should be defined outside and passed to the function?), including the twist factors
	Gauge_Conf aux_conf;
	init_gauge_conf_from_gauge_conf(&aux_conf, GC, param); // now aux_conf=GC

	// translation in direction +dir, including the twist factors
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < (param->d_n_planes + 1) * (param->d_volume); s++)
		{
		// s = j * volume + r
		long const r = s % (param->d_volume);
		int const j = (int) ((s - r) / (param->d_volume));
		if(j < STDIM) equal(&(GC->lattice[r][j]), &(aux_conf.lattice[nnm(geo, r, dir)][j]));
		GC->Z[r][j] = aux_conf.Z[nnm(geo, r, dir)][j];
		}

	// free auxiliary conf, including the twist factors
	free_gauge_conf(&aux_conf, param);
	}

double action_with_defect(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param, int a, int b, long r, int i, int j)
	{
	double K_a = (GC[a].C[r][i]) * (GC[a].C[nnp(geo, r, i)][j]) * (GC[a].C[nnp(geo, r, j)][i]) * (GC[a].C[r][j]);
	double K_b = (GC[b].C[r][i]) * (GC[b].C[nnp(geo, r, i)][j]) * (GC[b].C[nnp(geo, r, j)][i]) * (GC[b].C[r][j]);
	double re_tr_plaq_a = plaquettep(&(GC[a]), geo, param, r, i, j);
	double re_tr_plaq_b = plaquettep(&(GC[b]), geo, param, r, i, j);
	return -param->d_beta * (K_a * re_tr_plaq_a + K_b * re_tr_plaq_b);
	}

void print_complex(double complex const z)
	{
	printf("% .12f%+.12fi", creal(z), cimag(z));
	}

void print_twist_factor(int const i, int const j, long const r, double complex const z)
	{
	printf("Z_%d%d[%ld] = % .12e %+.12ei\n", i, j, r, creal(z), cimag(z));
	}

void print_plaquette(int const i, int const j, long const r, double complex const z)
	{
	printf("P_%d%d[%ld] = % .12e %+.12ei\n", i, j, r, creal(z), cimag(z));
	}

int main(int argc, char **argv)
	{
	if(argc != 2)
		{
		printf("Usage: %s input_file\n\n", argv[0]);

		print_compilation_details();

		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");

	real_main(argv[1]);

	return EXIT_SUCCESS;
	}

#endif
