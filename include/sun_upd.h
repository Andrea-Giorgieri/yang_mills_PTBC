#ifndef SUN_UPD_H
#define SUN_UPD_H

#include "macro.h"

#include "gparam.h"
#include "su2.h"
#include "su2_upd.h"
#include "sun.h"

#include <complex.h>
#include <math.h>


// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
static inline void single_heatbath_SuN(SuN *link, SuN const *const staple, GParam const *const param)
	{
	SuN aux;
	Su2 u, v, w;
	double xi, p0;
	double complex temp0, temp1;
	double complex fii, fij, fji, fjj;

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			equal_SuN(&aux, staple);     // aux=staple
			times_equal_SuN(&aux, link); // aux=staple*link
			ennetodue(&aux, i, j, &xi, &u);

			xi *= param->d_beta * 2.0 / (double) NCOLOR;

			if(xi > MIN_VALUE)
				{
				randheat_Su2(xi, &p0);

				equal_dag_Su2(&w, &u); // w=u^{dag}
				rand_matrix_p0_Su2(p0, &v);
				times_equal_Su2(&w, &v); // w*=v

				fii = w.comp[0] + w.comp[3] * I;
				fij = w.comp[2] + w.comp[1] * I;
				fji = -w.comp[2] + w.comp[1] * I;
				fjj = w.comp[0] - w.comp[3] * I;

				// link*=final
				for(int k = 0; k < NCOLOR; k++)
					{
					temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
					temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
					link->comp[m(k, i)] = temp0;
					link->comp[m(k, j)] = temp1;
					}
				}
			else
				{
				rand_matrix_Su2(&w);

				fii = w.comp[0] + w.comp[3] * I;
				fij = w.comp[2] + w.comp[1] * I;
				fji = -w.comp[2] + w.comp[1] * I;
				fjj = w.comp[0] - w.comp[3] * I;

				// link*=final
				for(int k = 0; k < NCOLOR; k++)
					{
					temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
					temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
					link->comp[m(k, i)] = temp0;
					link->comp[m(k, j)] = temp1;
					}
				}
			}
		}
	}


// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
static inline void single_heatbath_aux_SuN(SuN *link, SuN const *const staple, double const beta)
	{
	SuN aux;
	Su2 u, v, w;
	double xi, p0;

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			equal_SuN(&aux, staple);     // aux=staple
			times_equal_SuN(&aux, link); // aux=staple*link
			ennetodue(&aux, i, j, &xi, &u);

			xi *= beta * 2.0 / (double) NCOLOR;

			if(xi > MIN_VALUE)
				{
				randheat_Su2(xi, &p0);

				equal_dag_Su2(&w, &u); // w=u^{dag}
				rand_matrix_p0_Su2(p0, &v);
				times_equal_Su2(&w, &v); // w*=v

				double complex const fii = w.comp[0] + w.comp[3] * I;
				double complex const fij = w.comp[2] + w.comp[1] * I;
				double complex const fji = -w.comp[2] + w.comp[1] * I;
				double complex const fjj = w.comp[0] - w.comp[3] * I;

				// link*=final
				for(int k = 0; k < NCOLOR; k++)
					{
					double complex const temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
					double complex const temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
					link->comp[m(k, i)] = temp0;
					link->comp[m(k, j)] = temp1;
					}
				}
			}
		}
	}


// Pseudo-overrelaxation by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
static inline void single_overrelaxation_SuN(SuN *link, SuN const *const staple)
	{
	SuN prod;
	Su2 u, v;
	double xi;

	for(int i = 0; i < NCOLOR - 1; i++)
		{
		for(int j = i + 1; j < NCOLOR; j++)
			{
			equal_SuN(&prod, staple);     // prod=staple
			times_equal_SuN(&prod, link); // prod=staple*link
			ennetodue(&prod, i, j, &xi, &u);

			if(xi > MIN_VALUE)
				{
				equal_dag_Su2(&v, &u);       // v=u^{dag}
				times_equal_dag_Su2(&v, &u); // v=(u^{dag})^2

				double complex const fii = v.comp[0] + v.comp[3] * I;
				double complex const fij = v.comp[2] + v.comp[1] * I;
				double complex const fji = -v.comp[2] + v.comp[1] * I;
				double complex const fjj = v.comp[0] - v.comp[3] * I;

				//link*=final
				for(int k = 0; k < NCOLOR; k++)
					{
					double complex const temp0 = link->comp[m(k, i)] * fii + link->comp[m(k, j)] * fji;
					double complex const temp1 = link->comp[m(k, i)] * fij + link->comp[m(k, j)] * fjj;
					link->comp[m(k, i)] = temp0;
					link->comp[m(k, j)] = temp1;
					}
				}
			}
		}
	}


// TODO: bad unitarize for testing, remove
// unitarize with heatbath at beta O(1e20) saving bad links to file fp if print_flag != 0,
// returns phase gained during unitarization (bad links if != 0)
static inline double bad_unitarize_SuN(SuN *restrict A, double const beta, FILE *fp, int const print_flag)
	{
	double center_element = 0.0; // k such that gained phase is C = exp(-i 2pi/N k)
	SuN F;                       // F = A^{dag}, force to unitarize A by cooling
	SuN G, G_old;                // current and previous guess for unitarized A
	SuN H, H_copy, H_square;     // helpers to check convergence of unitarization
	SuN prod;

	// check if A needs re-unitarization: check_SuN(A) passes (=0) if
	// |A * A^{dag} - 1| < MIN_VALUE and |det(A) - 1| < MIN_VALUE
	if(scheck_SuN(A) == 1)
		{
		// use A^{dag} as force
		equal_dag_SuN(&F, A);

		// guess initialized to identity
		one_SuN(&G);
		double check = 1.0;
		while(check > MIN_VALUE)
			{
			// store old guess
			equal_SuN(&G_old, &G);

			// get new guess by heatbath
			single_heatbath_aux_SuN(&G, &F, beta); // maximize Tr(G*F) in large-beta limit

			// calculate the distance between old guess G_old and new guess G:
			// check = sqrt(|ReTr[(G-G_old)^2]|/N^2)
			equal_SuN(&H, &G);
			minus_equal_SuN(&H, &G_old);
			equal_SuN(&H_copy, &H);
			times_SuN(&H_square, &H, &H_copy);
			check = sqrt(fabs(retr_SuN(&H_square)) / (double) NCOLOR);
			}

		// Maximize ReTr[staple * C *link] for C \in Z(SU(N)) and update link *= C
		// \phi \equiv carg(Tr[staple * link]) => C = \exp{-i * 2\pi/N * round(\phi / (2*\pi/N))}
		equal_SuN(&prod, &F);                                   // prod=staple
		times_equal_SuN(&prod, &G);                             // prod=staple*link
		center_element = argtr_SuN(&prod);                      // center_element = phi
		center_element = round(center_element / PI2_N) * PI2_N; // round center_element to nearest center phase (PI2_N = 2*pi/N in marco.h)
		if(print_flag != 0 && fabs(center_element) > MIN_VALUE) // bad link: center_element != 0
			print_on_file_SuN(fp, A);

		// replace A with G
		equal_SuN(A, &G);
		}
	return center_element;
	}



#endif
