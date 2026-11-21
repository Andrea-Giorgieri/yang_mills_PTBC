#ifndef DEBUG_ACC_ROUTINES_C
#define DEBUG_ACC_ROUTINES_C

// Check of the OpenACC device routines.
// nvc generates the device version of a 'static inline' function marked with ACC_ROUTINE_SEQ only if it is
// called inside a compute region: this program calls all of them inside one compute region, so that a marked
// function calling a non-marked one (e.g. I/O or a host-only generator) is a compilation error.
// When run on a GPU it also executes the kernel and checks that the results are finite.

#include "../include/macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"

int main(void)
	{
	#ifndef OPENACC_MODE
	printf("This program checks the OpenACC device routines: configure with --enable-use-openacc\n");
	return EXIT_SUCCESS;
	#else
	GParam param = {0};
	param.d_beta = 6.0;
	param.d_volume = 1;
	for(int i = 0; i < STDIM; i++)
		{
		param.d_size[i] = 4;
		param.d_volume *= param.d_size[i];
		}

	long const n = 256;
	double *res = (double *) malloc((size_t) n * sizeof(double));
	REQUIRE(res != NULL, "failed to allocate memory");

	#pragma acc parallel loop copyin(param) copyout(res[0:n])
	for(long s = 0; s < n; s++)
		{
		rng_t rng;
		GAUGE_GROUP A, B, C;
		double x = 0.0;

		// random numbers
		#ifdef RNG_PHILOX_MODE
		rng_init_stream(&rng, 1u, 0u, (uint32_t) s);
		#endif
		x += casuale(&rng) + gauss1(&rng);
		double g1, g2;
		gauss2(&g1, &g2, &rng);
		x += g1 + g2;

		// algebra of the gauge group
		rand_matrix(&A, &rng);
		rand_matrix(&B, &rng);
		one(&C);
		zero(&C);
		equal(&C, &A);
		equal_dag(&C, &A);
		plus_equal(&C, &B);
		plus_equal_dag(&C, &B);
		minus_equal(&C, &B);
		minus_equal_times_real(&C, &B, 0.5);
		minus_equal_dag(&C, &B);
		times_equal_real(&C, 0.5);
		times_equal_complex(&C, 1.0 + 0.0 * I);
		times_equal(&C, &A);
		times_equal_dag(&C, &A);
		times(&C, &A, &B);
		times_dag1(&C, &A, &B);
		times_dag2(&C, &A, &B);
		times_dag12(&C, &A, &B);
		lin_comb(&C, 0.5, &A, 0.5, &B);
		lin_comb_dag1(&C, 0.5, &A, 0.5, &B);
		lin_comb_dag2(&C, 0.5, &A, 0.5, &B);
		lin_comb_dag12(&C, 0.5, &A, 0.5, &B);
		x += norm(&C) + relative_dist(&A, &B) + retr(&C) + imtr(&C);
		unitarize(&C);
		ta(&C);
		taexp(&C);
		x += retr(&C);

		// updates
		single_heatbath(&A, &B, &param, &rng);
		single_overrelaxation(&A, &B, &rng);
		cool(&C, &B);
		x += retr(&A) + retr(&C);

		// functions of the specific group not covered by the wrappers above
		#if NCOLOR > 2
		x += argtr_SuN(&A) + creal(tr_times_dag_SuN(&A, &B)) + creal(det_SuN(&A));
		x += (double) scheck_SuN(&A);
		bad_cool_SuN(&C, &B);
		rand_algebra_gauss_matrix_SuN(&C, &rng);
		x += (double) ta_check_SuN(&C);
		exp_of_ta_SuN(&C);
		single_heatbath_aux_SuN(&A, &B, param.d_beta, &rng);
		x += retr(&A) + retr(&C);
		#endif

		// indexing functions that need only the lattice sizes
		int cartcoord[STDIM];
		long const r = s % param.d_volume;
		si_to_cart(cartcoord, r, &param);
		x += (double) cart_to_si(cartcoord, &param) + (double) lex_to_si(r, &param) + (double) si_to_lex(r, &param);

		res[s] = x;
		}

	int finite = 1;
	double sum = 0.0;
	for(long s = 0; s < n; s++)
		{
		if(!isfinite(res[s])) finite = 0;
		sum += res[s];
		}
	printf("device routines executed: sum of results = %.12g, all finite: %s\n", sum, finite ? "yes" : "NO");

	free(res);
	return finite ? EXIT_SUCCESS : EXIT_FAILURE;
	#endif
	}

#endif
