#ifndef U1_UPD
#define U1_UPD

#include "macro.h"

#include <complex.h>
#include <math.h>

#include "gparam.h"
#include "random.h"
#include "u1.h"


// heatbath (see Moriarty Phys. Rev. D 25, p2185 (1982) )
static inline void single_heatbath_U1(U1 *link, U1 const *const staple, GParam const *const param)
	{
	double const alpha = param->d_beta * norm_U1(staple);

	if(alpha > MIN_VALUE)
		{
		double const exp2am1 = expm1(2.0 * alpha);
		double const q_max = exp(0.210513662353018684327769 * alpha); // see later for this number
		double const theta_staple = atan2(-cimag(staple->comp), creal(staple->comp));
		double x, q;
		do
			{
			x = -1.0 + log1p(exp2am1 * casuale()) / alpha;
			q = exp(alpha * (cos(HALF_PI * (1.0 - x)) - x));
			} while(q <= q_max * casuale());

		double theta = (1.0 - x) * HALF_PI;

		if(casuale() > 0.5)
			{
			theta = -theta;
			}
		theta += theta_staple;

		link->comp = cos(theta) + I * sin(theta);
		}
	else
		{
		rand_matrix_U1(link);
		}
	}


// overrelaxation
static inline void single_overrelaxation_U1(U1 *link, U1 const *const staple)
	{
	U1 newlink, helper;

	double const p = norm_U1(staple);

	if(p > MIN_VALUE)
		{
		equal_U1(&helper, staple);
		times_equal_real_U1(&helper, 1.0 / p);

		times_dag12_U1(&newlink, &helper, link);
		times_equal_dag_U1(&newlink, &helper);

		equal_U1(link, &newlink);
		}
	else
		{
		rand_matrix_U1(link);
		}
	}


// cooling
static inline void cool_U1(U1 *link, U1 const *const staple)
	{
	U1 aux;

	equal_U1(&aux, staple);
	unitarize_U1(&aux);
	equal_dag_U1(link, &aux);
	}


/*
WorkingPrecision->1000;

Q[x_]:=Exp[Cos[Pi/2*(1-x)]-x]

FindRoot[D[Q[x],x]==0, {x,1}, WorkingPrecision->100]

Out[52]= {x -> 0.5606641805798867176366776048997096707812104519411362714885751166519976969907076829844764496162691569}


N[Log[Q[x]/.%41], 100]

Out[53]= 0.210513662353018684327769435155832317434879346989632455087165428289411464536813734686515674410131331
*/


#endif
