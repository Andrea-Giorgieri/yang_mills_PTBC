#ifndef GAUGE_CONF_UPDATE_MULTICANONIC_C
#define GAUGE_CONF_UPDATE_MULTICANONIC_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
#include<omp.h>
#endif
#include<stdlib.h>

#include"../include/memalign.h"
#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// initialize stored topcharge and lattice_copy_cold
void init_multicanonic_gauge_conf(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param)
	{
	long s;
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice_cold), param->d_volume, __FILE__, __LINE__);
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice_copy_cold), param->d_volume, __FILE__, __LINE__);
	for(s=0; s<(param->d_volume); s++)
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice_cold[s]), STDIM, __FILE__, __LINE__);
		allocate_array_GAUGE_GROUP(&(GC->lattice_copy_cold[s]), STDIM, __FILE__, __LINE__);
		}
	
	// TODO: debug cold topcharge multicanonic
	// initialize lattice_cold = lattice_copy_cold = lattice
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s)
	#endif
	for(s=0; s<STDIM*(param->d_volume); s++)
		{
		long r = s % (param->d_volume);
		int i = (int) ( (s - r) / (param->d_volume) );
		equal(&(GC->lattice_cold[r][i]), &(GC->lattice[r][i]));
		equal(&(GC->lattice_copy_cold[r][i]), &(GC->lattice[r][i]));
		}
	
	// cool lattice_cold
	switch (param->d_topo_cooling)
		{
		case 0:	// topcharge not cooled
			GC->stored_topcharge = topcharge(GC, geo, param);
			break;
		case 1:	// topcharge cooled
			GC->stored_topcharge = multicanonic_topcharge_cooling(GC, geo, param);
			break;
		default:
			fprintf(stderr, "Undefined cooling method %d! (%s, %d)\n", param->d_topo_cooling, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

	// move lattice_cold to lattice_copy_cold, lattice_cold initialized before topcharge calculation
	GAUGE_GROUP **aux;
	aux = GC->lattice_cold;
	GC->lattice_cold = GC->lattice_copy_cold;
	GC->lattice_copy_cold = aux;
	}

// initialize multicanonic acceptances arrays and file
void init_multicanonic_acc_utils(Acc_Utils *acc_counters, GParam const * const param)
	{
	int i;

	// init acc arrays
	allocate_array_long(&(acc_counters->num_accepted_metro_multicanonic), param->d_N_replica_pt, __FILE__, __LINE__);
	allocate_array_long(&(acc_counters->num_metro_multicanonic), param->d_N_replica_pt, __FILE__, __LINE__);
	for(i=0;i<(param->d_N_replica_pt);i++) 
		{
		acc_counters->num_accepted_metro_multicanonic[i]=0;
		acc_counters->num_metro_multicanonic[i]=0;
		}
		
	// init acc file
	acc_counters->multicanonic_acc_filep=fopen(param->d_multicanonic_acc_file, "r");
	if((acc_counters->multicanonic_acc_filep)!=NULL) // file exists
		{
		fclose(acc_counters->multicanonic_acc_filep);
		acc_counters->multicanonic_acc_filep=fopen(param->d_multicanonic_acc_file, "a");
		}
	else // file doesn't exist, write first line
		{
		acc_counters->multicanonic_acc_filep=fopen(param->d_multicanonic_acc_file, "w");
		fprintf(acc_counters->multicanonic_acc_filep, "# %f ", param->d_beta);
		for(i=0; i<STDIM; i++) fprintf(acc_counters->multicanonic_acc_filep, "%d ", param->d_size[i]);
		fprintf(acc_counters->multicanonic_acc_filep, "\n");
		}
	fflush(acc_counters->multicanonic_acc_filep);
	}
		
void free_multicanonic_acc_utils(Acc_Utils *acc_counters)
	{
	free(acc_counters->num_accepted_metro_multicanonic);
	free(acc_counters->num_metro_multicanonic);
	fclose(acc_counters->multicanonic_acc_filep);
	}

// print metropolis acceptance of multicanonic algorithm on file
void print_multicanonic_acceptance(Gauge_Conf const * const GC, GParam const * const param, Acc_Utils const * const acc_counters)
	{
	double mean_acc;
	if(GC->update_index % param->d_measevery == 0)
		{
		fprintf(acc_counters->multicanonic_acc_filep, "%9ld ", GC->update_index);
		for(int i=0; i<param->d_N_replica_pt; i++)
			{
			if (acc_counters->num_metro_multicanonic[i] != 0) mean_acc = ( (double) acc_counters->num_accepted_metro_multicanonic[i] ) / ( (double) acc_counters->num_metro_multicanonic[i] );
			else mean_acc = 0.0;
			fprintf(acc_counters->multicanonic_acc_filep, "%7.3f ", 100.0*mean_acc);
			}
		fprintf(acc_counters->multicanonic_acc_filep, "\n");
		fflush(acc_counters->multicanonic_acc_filep);
		}
	}

// compute topo potential V_a(x)
double compute_topo_potential(int const a, double x, GParam const * const param)
	{		
	int i_grid;
	double x0, m, q;

	//alpha-rounding of topcharge
	if(param->d_topo_alpha > MIN_VALUE)  x = round(param->d_topo_alpha*x);
	
	// find index of nearest grid point to x
	i_grid=(int)(floor((x+param->d_grid_max)/param->d_grid_step));

	if(i_grid>=0 && i_grid<param->d_n_grid) // if x inside the barriers compute V_a(x) with a linear interpolation
		{
		// perform linear interpolation: V_a(x) = V_a(x0) + [dV_a/dx|(x0)] (x-x0)
		x0 = i_grid * param->d_grid_step - param->d_grid_max;
		m = (param->d_grid[a][i_grid+1] - param->d_grid[a][i_grid]) / param->d_grid_step; // dV/dx|(x0) = [ V(x0+step) - V(x0) ] / step
		q = param->d_grid[a][i_grid]-m*x0; // V(x0) - [dV/dx|(x0)] x0
		return q + m * x; // V(x0) + [dV/dx|(x0)] (x-x0) 
		}
	else // if x outside the barriers just saturate to extreme values
		{
		if(i_grid<0) return param->d_grid[a][0];
		else return param->d_grid[a][param->d_n_grid-1];
		}
	}

// compute the difference of the local topological charge for two gauge confs
double multicanonic_delta_loc_topcharge(Gauge_Conf const * const GC1,
							Gauge_Conf const * const GC2,
							Geometry const * const geo,
							GParam const * const param,
							long r)
	{
	double ris;
	
	#if (STDIM==4 && NCOLOR>1)
	GAUGE_GROUP aux1, aux2, aux3;
	double real1, real2, loc_charge;
	const double chnorm=1.0/(128.0*PI*PI);
	int i, dir[4][3], sign;
	
	dir[0][0] = 0;
	dir[0][1] = 0;
	dir[0][2] = 0;
	
	dir[1][0] = 1;
	dir[1][1] = 2;
	dir[1][2] = 3;
	
	dir[2][0] = 2;
	dir[2][1] = 1;
	dir[2][2] = 1;
	
	dir[3][0] = 3;
	dir[3][1] = 3;
	dir[3][2] = 2;
	
	sign=-1;
	loc_charge=0.0;
	for(i=0; i<3; i++)
		{
		// GC1 contribution (with +)
		clover(GC1, geo, param, r, dir[0][i], dir[1][i], &aux1);
		clover(GC1, geo, param, r, dir[2][i], dir[3][i], &aux2);
		
		times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
		real1=retr(&aux3)*NCOLOR;
		
		times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
		real2=retr(&aux3)*NCOLOR;
		
		loc_charge+=((double) sign*(real1-real2));
	
		// GC2 contribution (with -)
		clover(GC2, geo, param, r, dir[0][i], dir[1][i], &aux1);
		clover(GC2, geo, param, r, dir[2][i], dir[3][i], &aux2);
		
		times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
		real1=retr(&aux3)*NCOLOR;
		
		times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
		real2=retr(&aux3)*NCOLOR;
		
		loc_charge-=((double) sign*(real1-real2));
		
		sign=-sign;
		}
	ris = (loc_charge*chnorm);
	
	// TODO: debug, remove
	//fprintf(stderr, "%ld % f ", r, ris);
	#endif
	
	#if (STDIM==2 && NCOLOR==1)
	GAUGE_GROUP u1matrix;
	double angle;
	
	// GC1 contribution (with +)
	plaquettep_matrix(GC1, geo, param, r, 0, 1, &u1matrix);
	angle = atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;
	
	// GC2 contribution (with -)
	plaquettep_matrix(GC2, geo, param, r, 0, 1, &u1matrix);
	angle -= atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;
	
	ris = angle;
	#endif
	
	return ris;
	}

// compute the difference of the topological charge of a rectangle
double multicanonic_delta_topcharge_rectangle(Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							Rectangle const * const topcharge_rect)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	Gauge_Conf helper;
	double ris;
	long n;
	
	helper.lattice = GC->lattice_copy;
	helper.Z = GC->Z_copy;
	ris=0.0;

	// TODO: benchmark number of threads, often rectangles are small
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(n) reduction(+ : ris)
	#endif
	for(n=0; n<(topcharge_rect->d_vol_rect); n++)
		{
		long r = topcharge_rect->rect_sites[n];
		ris += multicanonic_delta_loc_topcharge(GC, &helper, geo, param, r);
		}

	// TODO: remove, debug only
	/*
	for(n=0; n<(param->d_volume); n++)
		{
		ris += multicanonic_delta_loc_topcharge(GC, &helper, geo, param, n);
		int is_on_rect=0;
		long r;
		for(r=0; r<(topcharge_rect->d_vol_rect); r++)
			{
			if (topcharge_rect->rect_sites[r] == n) is_on_rect = 1;
			}
		fprintf(stderr, "%d\n", is_on_rect);
		}
	*/

	return ris;
	}

// compute the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
double multicanonic_topcharge_cooling(Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	if(param->d_topo_coolsteps>0)
		{
		Gauge_Conf helper;
		
		// TODO: gcc gives warning
		//equal_lattice(GC->lattice_cold, (GAUGE_GROUP const * const * const)GC->lattice, param);
		equal_lattice(GC->lattice_cold, GC->lattice, param);
		// TODO: unitarize before cooling?
		helper.lattice = GC->lattice_cold;
		helper.Z = GC->Z;
		cooling(&helper, geo, param, param->d_topo_coolsteps);
		return topcharge(&helper, geo, param);
		}
	else
		{
		return topcharge(GC, geo, param);
		}
	}

// compute the difference of the topological charge after some cooling only on a rectangle
// in the cooling procedure the action at theta=0 is minimized
double multicanonic_delta_topcharge_cooling_rectangle(Gauge_Conf * const GC, 
										Geometry const * const geo,
										GParam const * const param,
										int const hierarc_level,
										Rect_Utils const * const rect_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	// TODO: debug
	if(param->d_topo_coolsteps>0)
		{
		Gauge_Conf helper;
		
		// TODO: gcc gives warning
		//equal_lattice(GC->lattice_cold, (GAUGE_GROUP const * const * const)GC->lattice, param);
		equal_lattice(GC->lattice_cold, GC->lattice, param);
		
		helper.lattice = GC->lattice_cold;
		helper.Z = GC->Z;
		
		helper.lattice_copy = GC->lattice_copy_cold;
		helper.Z_copy = GC->Z_copy;
		
		//hierarchical_cooling(&helper, geo, param, rect_aux->cooling_rect[hierarc_level]);
		cooling(&helper, geo, param, param->d_topo_coolsteps);
		return multicanonic_delta_topcharge_rectangle(&helper, geo, param, &(rect_aux->topcharge_rect[hierarc_level]));
		}
	else
		{
		return multicanonic_delta_topcharge_rectangle(GC, geo, param, &(rect_aux->topcharge_rect[hierarc_level]));
		}
	}

// compute the topological charge after some gradflow time
// TODO: refactor or remove, very inefficient, currently unused
double topcharge_agf_multicanonic(Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
		{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	if (param->d_topo_agf_time > 0.0)	//if using gradient flow
		{
		int accepted;
		double gftime, gftime_step, charge;
		
		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		while(gftime < param->d_topo_agf_time)
			{
			gradflow_RKstep_adaptive(GC, geo, param, &gftime, &gftime_step, &accepted, meas_aux);
			
			// adapt step to the time of the measure
			if ((gftime + gftime_step - param->d_topo_agf_time) > param->d_agf_time_bin + MIN_VALUE)
				{
				gftime_step = param->d_topo_agf_time - gftime;
				}
			}
		charge = topcharge(GC, geo, param);
		restore_gauge_conf(GC, param);
		
		return charge; 
		}
	else	// no gradient flow
		{
		return topcharge(GC, geo, param);
		}
	}

// compute the multicanonic Metropolis probability p=exp(delta V_a) where V_a is the topo potential for replica a
double metropolis_prob_multicanonic(int const a, double const Q_new, double const Q_old, GParam const * const param)
	{
	double V_old, V_new;
	
	V_old = compute_topo_potential(a, Q_old, param);
	V_new = compute_topo_potential(a, Q_new, param);
	//return 1.1;
	return exp(V_old-V_new);
	}

// compute variation of topo potential swapping two replicas
double delta_topo_potential_swap(Gauge_Conf const * const GC, int const a, int const b, GParam const * const param)
	{
	double V_noswap, V_swap, Q_a, Q_b;
	
	Q_a = GC[a].stored_topcharge;
	Q_b = GC[b].stored_topcharge;
	
	// before swap V_a(Q_a) + V_b(Q_b) , after swap V_a(Q_b) + V_b(Q_a)
	V_noswap = compute_topo_potential(a, Q_a, param) + compute_topo_potential(b, Q_b, param);
	V_swap   = compute_topo_potential(a, Q_b, param) + compute_topo_potential(b, Q_a, param);

	// swapped potential - unswapped potential
	return V_swap - V_noswap;
	}


//	Metropolis test with p_metro=exp(- delta topo_potential)
int multicanonic_metropolis_step_all_links(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param, int const idx)
	{
	// perform multicanonic Metropolis test
	double Q_old, Q_new, p;
	
	Q_old = GC->stored_topcharge;
	switch (param->d_topo_cooling)
		{
		case 0:	// topcharge not cooled
			Q_new = topcharge(GC, geo, param);
			break;
		case 1:	// topcharge cooled
			Q_new = multicanonic_topcharge_cooling(GC, geo, param);
			break;
		default:
			fprintf(stderr, "Undefined cooling method %d! (%s, %d)\n", param->d_topo_cooling, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	
	p = metropolis_prob_multicanonic(idx, Q_new, Q_old, param);
	
	// Metropolis test: p < 1 --> acc=1 with probability p, else --> acc=1
	int acc = 1;
	if(p < 1)
		{
		double random_number=casuale();
		if(random_number > p) acc = 0;
		}
	
	// if Metropolis is accepted, store the new topological charge
	if (acc == 1) GC->stored_topcharge = Q_new;
	return acc;
	}

//	Metropolis test with p_metro=exp(- delta topo_potential) only on a given rectangle when using cold charge
int multicanonic_metropolis_step_rectangle(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
											int const hierarc_level, Rect_Utils const * const rect_aux, int const idx)
	{
	// perform multicanonic Metropolis test
	double Q_old, Q_new, p;
	//double Q_new_debug;
	
	Q_old = GC->stored_topcharge;
	switch (param->d_topo_cooling)
		{
		case 0:	// topcharge not cooled
			Q_new = Q_old + multicanonic_delta_topcharge_rectangle(GC, geo, param, &(rect_aux->topcharge_rect[hierarc_level]));
			//Q_new_debug = topcharge(GC, geo, param);
			break;
		case 1:	// topcharge cooled
			Q_new = Q_old + multicanonic_delta_topcharge_cooling_rectangle(GC, geo, param, hierarc_level, rect_aux);
			//Q_new_debug = multicanonic_topcharge_cooling(GC, geo, param);
			break;
		default:
			fprintf(stderr, "Undefined cooling method %d! (%s, %d)\n", param->d_topo_cooling, __FILE__, __LINE__);
			//exit(EXIT_FAILURE);
		}
	// TODO: debug, remove
	//fprintf(stdout, "%ld %d %f %f\n", GC->update_index, hierarc_level, Q_new, Q_new_debug);
	
	p = metropolis_prob_multicanonic(idx, Q_new, Q_old, param);
	
	// Metropolis test: p < 1 --> acc=1 with probability p, else --> acc=1
	int acc = 1;
	if(p < 1)
		{
		double random_number=casuale();
		if(random_number > p) acc = 0;
		}
	
	// if Metropolis is accepted, store the new topological charge
	if (acc == 1) GC->stored_topcharge = Q_new;
	return acc;
	}


/*********************************************************************************************************************/
/* OLD IMPLEMENTATION, WAITING TO BE REMOVED */
/*********************************************************************************************************************/




// compute single clover insertion and store it in <clover_insertion>
void compute_single_clover_insertion(GAUGE_GROUP * clover_insertion, Gauge_Conf const * const GC,
							Geometry const * const geo, GParam const * const param, long r, int i, int j)
	{
	GAUGE_GROUP aux;
	clover(GC, geo, param, r, i, j, &aux); // aux = clover[r][i][j]
	equal(clover_insertion, &aux); // clover_insertion = aux
	minus_equal_dag(clover_insertion, &aux);	// clover_insertion -= aux^{dag}
	// ==> clover_insertion = clover[r][i][j] - clover[r][i][j]^dag = 2i Im{ clover[r][i][j] }
	}

// compute staple of topological charge relative to link (r,i)
void compute_topostaple_alone(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, long r, int i, GAUGE_GROUP * topo_stap)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#endif
	
	if(STDIM!=4)
		{
		fprintf(stderr, "Error: topological charge can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	GAUGE_GROUP link1, link2, link3, link12, stap, aux, clover_insertion;
	const double coeff = 1.0/(128.0*PI*PI);
	long k;
	int j, l;
	int i0, j0;
	int sood1[4][4], sood2[4][4]; // Signed Ordered Orthogonal Directions (SOOD)
	
	zero(topo_stap); // topo_stap=0
	
	// the clover topological charge is written as
	// -1/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
	// where Q_{\mu\nu} here stands for the clover on plane (\mu \nu)
	// the independent permutations are 3, here we use: 0123 0231 0312
	
	sood1[0][1] = 2;
	sood2[0][1] = 3;
	sood1[1][0] = 3;
	sood2[1][0] = 2;
	
	sood1[0][2] = 3;
	sood2[0][2] = 1;
	sood1[2][0] = 1;
	sood2[2][0] = 3;
	
	sood1[0][3] = 1;
	sood2[0][3] = 2;
	sood1[3][0] = 2;
	sood2[3][0] = 1;
	
	sood1[1][2] = 0;
	sood2[1][2] = 3;
	sood1[2][1] = 3;
	sood2[2][1] = 0;
	
	sood1[1][3] = 2;
	sood2[1][3] = 0;
	sood1[3][1] = 0;
	sood2[3][1] = 2;
	
	sood1[2][3] = 0;
	sood2[2][3] = 1;
	sood1[3][2] = 1;
	sood2[3][2] = 0;
	
	for(l=i+1; l< i + STDIM; l++)
		{
		j = (l % STDIM);
		
		i0=sood1[i][j];
		j0=sood2[i][j];
		
//
//		 i  ^
//			|	(1)
//		(b) +----->-----+ (c)
//			|			|
//			|			|
//			|			V (2)
//			|			|
//			|			|
//		(a) +-----<-----+--> j
//          r    (3)    (d)
//
			
		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));	// link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));				// link3 = (3)
		
		times_dag2(&link12, &link1, &link2);	// link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);		// stap=link12*link3^{dag}
		
		//twist (clockwise plaquette) modification
		times_equal_complex(&stap, GC->Z[r][dirs_to_si(i,j)]);	// Z_\mu\nu(x) * staple
		
		// clover insertion in (a)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, r, i0, j0);
		times(&aux, &stap, &clover_insertion); // stap*clover
		plus_equal(topo_stap, &aux);
		
		// clover insertion in (b)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, i), i0, j0);
		times(&aux, &clover_insertion, &stap);	// clover*stap
		plus_equal(topo_stap, &aux);
		
		// clover insertion in (c)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, nnp(geo, r, i), j), i0, j0);
		times(&aux, &link1, &clover_insertion);	// link1*clover
		times_equal_dag(&aux, &link2);		 // *=link2^{dag}
		times_equal_dag(&aux, &link3);		 // *=link3^{dag}
		plus_equal(topo_stap, &aux);
		
		// clover insertion in (d)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, j), i0, j0);
		times(&aux, &link12, &clover_insertion);	// link1*link2*quadri
		times_equal_dag(&aux, &link3);			 // *=link3^{dag}
		plus_equal(topo_stap, &aux);

//
//		 i  ^
//			|	(1)
//		(d) +----<------+ (a)
//			|			|
//			|			|
//		(2) V			|
//			|			|
//			|			| (b)
//		(c) +------>----+--->j
//			k     (3)   r
//
		
		k=nnm(geo, r, j); // k = r - j
		
		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));	// link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));				// link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));				// link3 = (3)
		
		times_dag12(&link12, &link1, &link2); 	// link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3); 			// stap=link12*link3
		
		//twist (anticlockwise plaquette) modification
		times_equal_complex(&stap, GC->Z[k][dirs_to_si(j,i)]); // Z_\mu\nu(x-\nu) * staple
		
		// clover insertion in (a)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, i), i0, j0);
		times(&aux, &clover_insertion, &stap); 	// clover*stap
		minus_equal(topo_stap, &aux);
		
		// clover insertion in (b)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, r, i0, j0);
		times(&aux, &stap, &clover_insertion); 	// stap*clover
		minus_equal(topo_stap, &aux);
		
		// clover insertion in (c)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, k, i0, j0);
		times(&aux, &link12, &clover_insertion); 	// link1^{dag}*link2^{dag}*clover
		times_equal(&aux, &link3); 					// *=link3
		minus_equal(topo_stap, &aux);
		
		// clover insertion in (d)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, k, i), i0, j0);
		times_dag1(&aux, &link1, &clover_insertion); 	// link1^{dag}*clover
		times_equal_dag(&aux, &link2); 					// *=link2^{dag}
		times_equal(&aux, &link3); 						// *=link3
		minus_equal(topo_stap, &aux);
		}
	times_equal_real(topo_stap, coeff); 	// topo_stap *= coeff
	}

// compute the variation of the clover topological charge when the link (r,i) is updated starting from <old_link>
double delta_Q_upd(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, long r, int i, GAUGE_GROUP old_link)
	{
	// compute delta_Q
	GAUGE_GROUP q, topo_stap;
	double delta_Q = 0.0;
	
	// compute topo staple and store it in topo_stap
	compute_topostaple_alone(GC, geo, param, r, i, &topo_stap);
	
	// compute contribution of new link
	times(&q, &topo_stap, &(GC->lattice[r][i])); 	// q = (topo_stap * new_link) 
	delta_Q += retr(&q)*((double) NCOLOR); 			// delta_Q += [ retr(topo_stap * new_link) / NCOLOR ] * NCOLOR (retr automatically adds a 1/NCOLOR factor)
	
	// compute contribution of old_link
	times(&q, &topo_stap, &old_link); // q = (topo_stap * old_link)
	delta_Q -= retr(&q)*((double) NCOLOR); // delta_Q -= [ retr(topo_stap * old_link) / NCOLOR ] * NCOLOR
	
	return delta_Q;
	}


// perform a single step of parallel tempering with hierarchic update with a multicanonical approach
void multicanonic_parallel_tempering_with_hierarchical_update(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
												 Rectangle const * const most_update, Rectangle const * const clover_rectangle,
												 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)												 
	{
	// first hierarc level is 0
	int start_hierarc=0;
	
	// set multicanonic Metropolis acceptance counters to zero to compute mean acc over single updating step
	for(int i=0;i<param->d_N_replica_pt; i++)
		{
		acc_counters->num_accepted_metro_multicanonic[i] = 0;
		acc_counters->num_metro_multicanonic[i] = 0;
		}
	
	// Parallel tempering updating step: full update + hierarchical update + swaps and translations after every sweep for every replica
	multicanonic_update_with_defect(GC, geo, param, acc_counters); // full update of all replicas
	if(param->d_N_replica_pt>1)
		{
		swap(GC, geo, param, swap_rectangle, acc_counters); // swap all replicas
		conf_translation(&(GC[0]), geo, param); // translation of periodic replica (GC[0])
		if(param->d_N_hierarc_levels>0)
			multicanonic_hierarchical_update_rectangle_with_defect(GC, geo, param, start_hierarc, most_update, clover_rectangle,
					swap_rectangle, acc_counters); // hierarchic update
		}
	
	// increase update index of all replicas
	for(int i=0;i<param->d_N_replica_pt; i++) GC[i].update_index++;
	
	// print mean multicanonic acceptance over a single updating step
	print_multicanonic_acceptance(GC, param, acc_counters);
	}

// perform a single step of parallel tempering with hierarchic update with a multicanonical approach using gradflowed charge
void multicanonic_agf_parallel_tempering_with_hierarchical_update(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
												 Rectangle const * const most_update, Rectangle const * const clover_rectangle,
												 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)												 
	{
	// first hierarc level is 0
	int start_hierarc=0;
	
	// set multicanonic Metropolis acceptance counters to zero to compute mean acc over single updating step
	for(int i=0;i<param->d_N_replica_pt; i++)
		{
		acc_counters->num_accepted_metro_multicanonic[i] = 0;
		acc_counters->num_metro_multicanonic[i] = 0;
		}
	
	// Parallel tempering updating step: full update + hierarchical update + swaps and translations after every sweep for every replica
	multicanonic_agf_update_with_defect(GC, geo, param, acc_counters); // full update of all replicas
	if(param->d_N_replica_pt>1)
		{
		swap(GC, geo, param, swap_rectangle, acc_counters); // swap all replicas
		conf_translation(&(GC[0]), geo, param); // translation of periodic replica (GC[0])
		if(param->d_N_hierarc_levels>0)
			multicanonic_agf_hierarchical_update_rectangle_with_defect(GC, geo, param, start_hierarc, most_update, clover_rectangle,
					swap_rectangle, acc_counters); // hierarchic update
		}
	
	// increase update index of all replicas
	for(int i=0;i<param->d_N_replica_pt; i++) GC[i].update_index++;
	
	// print mean multicanonic acceptance over a single updating step
	print_multicanonic_acceptance(GC, param, acc_counters);
	}

// update all replica only on a given rectangle in the presence of a defect (MODIFICARE)
void multicanonic_update_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param, Acc_Utils *acc_counters)
	{
	long s, num_even, num_odd;
	int j, dir;
	int num_replica = param->d_N_replica_pt; // just an auxiliary variable
	long *sum_acc, *count_metro;
	
	for(j=0; j<STDIM; j++)
		{
		if(param->d_size[j]==1)
			{
			fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		}
	
	// init aux variables to compute mean multicanonic acc
	allocate_array_long(&sum_acc, num_replica, __FILE__, __LINE__);
	allocate_array_long(&count_metro, num_replica, __FILE__, __LINE__);
	for(j=0; j<param->d_N_replica_pt; j++)
		{
		sum_acc[j]=0;
		count_metro[j]=0;
		}
	
	num_even = (param->d_volume + (param->d_volume % 2)) / 2;
	num_odd  = (param->d_volume - (param->d_volume % 2)) / 2;
	
	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif
		for(s=0; s<((param->d_N_replica_pt)*num_even); s++)
			{
			// s = i * num_even + r
			long r = s % num_even; 				// site index
			int i = (int) ( (s-r) / num_even );	// replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			}
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif 
		for(s=0; s<((param->d_N_replica_pt)*num_odd); s++)
			{
			// s = i * num_odd + aux ; aux = r - num_even
			long aux = s % num_odd;
			long r = num_even + aux;				// site index
			int i = (int) ( (s-aux) / num_odd );	// replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			}
		}
	
	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif
		
		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<((param->d_N_replica_pt)*num_even); s++)
				{
				// s = i * num_even + r
				long r = s % num_even;				// site index
				int i = (int) ( (s-r) / num_even );	// replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<((param->d_N_replica_pt)*num_odd); s++)
				{
				// s = i * num_odd + aux ; aux = r - num_even
				long aux = s % num_odd;
				long r = num_even + aux;				// site index
				int i = (int) ( (s-aux) / num_odd );	// replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}
			}
		}
			
	// add number of accepted and number of proposed Metropolis multicanonic tests
	for(int i=0; i<param->d_N_replica_pt; i++)
		{
		acc_counters->num_accepted_metro_multicanonic[i] += sum_acc[i];
		acc_counters->num_metro_multicanonic[i] += count_metro[i];
		}
		
	// free aux arrays
	free(sum_acc);
	free(count_metro);
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)); s++)
		{
		// s = i * volume + r
		long r = s % param->d_volume;
		int i = (int) ( (s-r) / (param->d_volume) );
		for(dir=0; dir<STDIM; dir++)
			{
			unitarize(&(GC[i].lattice[r][dir]));
			}
		}
	}

// update all replica in the presence of a defect using gradflowed charge
void multicanonic_agf_update_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param, Acc_Utils *acc_counters)
	{
	long s, acc;
	int j, dir;
	
	for(j=0; j<STDIM; j++)
		{
		if(param->d_size[j]==1)
			{
			fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		}
	
	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif
		for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
			{
			// s = i * volume/2 + r
			long r = s % ( (param->d_volume)/2 ); // site index
			int i = (int) ( (s-r) / ( (param->d_volume)/2 ) ); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
			{
			// s = i * volume/2 + aux ; aux = r - volume/2
			long aux = s % ( (param->d_volume)/2 );
			long r = (param->d_volume/2) + aux; // site index
			int i = (int) ( (s-aux) / ( (param->d_volume)/2 ) ); // replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		}
	
	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica(GC, geo, param, dir);
		#endif
		
		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
				{
				// s = i * volume/2 + r
				long r = s % ( (param->d_volume)/2 ); // site index
				int i = (int) ( (s-r) / ( (param->d_volume)/2 ) ); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
				{
				// s = i * volume/2 + aux ; aux = r - volume/2
				long aux = s % ( (param->d_volume)/2 );
				long r = (param->d_volume/2) + aux; // site index
				int i = (int) ( (s-aux) / ( (param->d_volume)/2 ) ); // replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)); s++)
		{
		// s = i * volume + r
		long r = s % param->d_volume;
		int i = (int) ( (s-r) / (param->d_volume) );
		for(dir=0; dir<STDIM; dir++)
			{
			unitarize(&(GC[i].lattice[r][dir]));
			}
		}
	
	// multicanonic Metropolis tests and add number of accepted and number of proposed ones
	for(j=0; j<param->d_N_replica_pt; j++)
		{
		acc = multicanonic_metropolis_step_all_links(&(GC[j]), geo, param, j);
		acc_counters->num_accepted_metro_multicanonic[j] += acc;
		acc_counters->num_metro_multicanonic[j] += 1;
		}
	}

// TODO: BUGGED, conf_label -> conf index
//	Metropolis test with p_metro=exp(- delta topo_potential)
int multicanonic_metropolis_step_single_link(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
												long r, int i, GAUGE_GROUP old_link)
	{
	// perform multicanonic Metropolis test
	double Q_old = GC->stored_topcharge;
	double Q_new = Q_old + delta_Q_upd(GC, geo, param, r, i, old_link);
	double p = metropolis_prob_multicanonic(GC->conf_label, Q_new, Q_old, param);
	
	// Metropolis test: p < 1 --> acc=1 with probability p, else --> acc=1
	int acc = 1;
	if(p < 1)
		{
		double random_number=casuale();
		if(random_number > p) acc = 0;
		}
	
	if (acc == 0) GC->lattice[r][i] = old_link; // if Metropolis is refused go back to original link
	else GC->stored_topcharge = Q_new; 			// if Metropolis is accepted store the new topological charge of the conf
	
	return acc;
	}

// hierarchical update functions

// update all replica only on a given rectangle in the presence of a defect (MODIFICARE)
void multicanonic_update_rectangle_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
												Rectangle const * const most_update, Rectangle const * const clover_rectangle, Acc_Utils *acc_counters)
	{
	
	long s, num_even, num_odd;
	int j, dir;
	int num_replica = param->d_N_replica_pt; // just an auxiliary variable
	long *sum_acc, *count_metro;
	
	// init aux variables to compute mean multicanonic acc
	allocate_array_long(&sum_acc, num_replica, __FILE__, __LINE__);
	allocate_array_long(&count_metro, num_replica, __FILE__, __LINE__);
	
	for(int i=0; i<param->d_N_replica_pt; i++)
		{
		sum_acc[i]=0;
		count_metro[i]=0;
		}
	
	#ifndef THETA_MODE
	(void) clover_rectangle; // to avoid compiler warning of unused variable
	#endif	
	
	/* Check if there's at least one even dimension of the rectangle, i.e. check if d_vol_rect is even.
		If there's at least one even dimension: d_vol_rect/2 even sites and d_vol_rect/2 odd sites.
		Otherwise: (d_vol_rect+1)/2 even sites and (d_vol_rect-1)/2 odd sites. */
	
	long is_even = ( most_update->d_vol_rect ) % 2;
	
	num_even = ( most_update->d_vol_rect + is_even ) / 2;	// number of even sites
	num_odd	= ( most_update->d_vol_rect - is_even ) / 2;	// number of odd sites
	
	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif 
		for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_even + n
			long n = s % num_even;					// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-n) / num_even );		// replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			}
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif 
		for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_odd + aux; aux = n - num_even
			long aux = s % num_odd;
			long n = aux + num_even;				// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-aux) / num_odd );	// replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			} 
		}
	
	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif
		
		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_even + n
				long n = s % num_even;					// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-n) / num_even );		// replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}
			
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_odd + aux; aux = n - num_even
				long aux = s % num_odd; 
				long n = aux + num_even;				// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-aux) / num_odd );	// replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}
			}
		}
	
	// add number of accepted and number of proposed Metropolis multicanonic tests
	for(int i=0; i<param->d_N_replica_pt; i++)
		{
		acc_counters->num_accepted_metro_multicanonic[i] += sum_acc[i];
		acc_counters->num_metro_multicanonic[i] += count_metro[i];
		}
	
	// free aux arrays
	free(sum_acc);
	free(count_metro);
	
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((most_update->d_vol_rect)*(param->d_N_replica_pt)); s++)
		{
		for(dir=0; dir<STDIM; dir++)
			{
			// s = i * volume_rect + n
			long n = s % (most_update->d_vol_rect);					// site index on rectangle
			long r = most_update->rect_sites[n];					// site index on lattice
			int i = (int) ( (s-n) / (most_update->d_vol_rect) );	// replica index
			unitarize(&(GC[i].lattice[r][dir]));
			} 
		}
	}

// perform a hierarchical update on all rectangles
void multicanonic_hierarchical_update_rectangle_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
														 int const hierarc_level, 
														 Rectangle const * const most_update,
														 Rectangle const * const clover_rectangle,
														 Rectangle const * const swap_rectangle,
														 Acc_Utils *acc_counters)
	{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
			{
			multicanonic_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]), acc_counters);
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			}
		}
	else
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
			{
			multicanonic_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]), acc_counters);	
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			multicanonic_hierarchical_update_rectangle_with_defect(GC,geo,param,(hierarc_level+1),most_update,clover_rectangle, swap_rectangle,acc_counters);
			}
		}
	}

// update all replica only on a given rectangle in the presence of a defect using gradflowed charge
void multicanonic_agf_update_rectangle_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
													Rectangle const * const most_update,
													Rectangle const * const clover_rectangle,
													Acc_Utils *acc_counters)
	{
	long s, acc, num_even, num_odd;
	int j, dir;
	
	#ifndef THETA_MODE
	(void) clover_rectangle; // to avoid compiler warning of unused variable
	#endif	
	
	/* Check if there's at least one even dimension of the rectangle, i.e. check if d_vol_rect is even.
		If there's at least one even dimension: d_vol_rect/2 even sites and d_vol_rect/2 odd sites.
		Otherwise: (d_vol_rect+1)/2 even sites and (d_vol_rect-1)/2 odd sites. */
	
	long is_even = ( most_update->d_vol_rect ) % 2;
	
	num_even = ( most_update->d_vol_rect + is_even ) / 2;	// number of even sites
	num_odd	= ( most_update->d_vol_rect - is_even ) / 2;	// number of odd sites
	
	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_even + n
			long n = s % num_even;					// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-n) / num_even );		// replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			}
		
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s)
		#endif 
		for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_odd + aux; aux = n - num_even
			long aux = s % num_odd;
			long n = aux + num_even;				// site index on rectangle
			long r = most_update->rect_sites[n];	// site index on lattice
			int i = (int) ( (s-aux) / num_odd );	// replica index
			heatbath_with_defect(&(GC[i]), geo, param, r, dir);
			} 
		}
	
	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
		compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif
		
		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_even + n
				long n = s % num_even;					// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-n) / num_even );		// replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s)
			#endif 
			for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_odd + aux; aux = n - num_even
				long aux = s % num_odd; 
				long n = aux + num_even;				// site index on rectangle
				long r = most_update->rect_sites[n];	// site index on lattice
				int i = (int) ( (s-aux) / num_odd );	// replica index
				overrelaxation_with_defect(&(GC[i]), geo, param, r, dir);
				}
			}
		}

	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((most_update->d_vol_rect)*(param->d_N_replica_pt)); s++)
		{
		for(dir=0; dir<STDIM; dir++)
			{
			// s = i * volume_rect + n
			long n = s % (most_update->d_vol_rect);					// site index on rectangle
			long r = most_update->rect_sites[n];					// site index on lattice
			int i = (int) ( (s-n) / (most_update->d_vol_rect) );	// replica index
			unitarize(&(GC[i].lattice[r][dir]));
			} 
		}

	// multicanonic Metropolis tests and add number of accepted and number of proposed ones
	for(j=0; j<param->d_N_replica_pt; j++)
		{
		// TODO: It should be possible to calculate Q_new as Q_old + delta_Q with delta_Q calculated
		//	  around the rectangle. At each gradflow step the rectangle in which links change needs to
		//	  extend to the first neighbors. delta_Q is non-zero only for those clovers intersecting
		//	  the final rectangle
		acc = multicanonic_metropolis_step_all_links(&(GC[j]), geo, param, j);
		acc_counters->num_accepted_metro_multicanonic[j] += acc;
		acc_counters->num_metro_multicanonic[j] += 1;
		}
	}

// perform a hierarchical update on all rectangles using gradflowed charge
void multicanonic_agf_hierarchical_update_rectangle_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
														 int const hierarc_level, 
														 Rectangle const * const most_update,
														 Rectangle const * const clover_rectangle,
														 Rectangle const * const swap_rectangle,
														 Acc_Utils *acc_counters)
	{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
			{
			multicanonic_agf_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]), acc_counters);
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			}
		} // end if
	else
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
			{
			multicanonic_agf_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]), acc_counters);	
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			multicanonic_agf_hierarchical_update_rectangle_with_defect(GC,geo,param,(hierarc_level+1),most_update,clover_rectangle, swap_rectangle,acc_counters);
			}
		} // end else
	}

// perform an update with overrelaxation and multicanonic Metropolis step
int multicanonic_overrelaxation_with_defect(Gauge_Conf * const GC,
							Geometry const * const geo,
							GParam const * const param,
							long r,
							int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#endif
	(void) param; // just to avoid warnings
	
	GAUGE_GROUP stap;
	
	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
	#endif
	
	// store old link before update
	GAUGE_GROUP old_link = GC->lattice[r][i];
	
	// perform usual single-link over-relaxation update
	single_overrelaxation(&(GC->lattice[r][i]), &stap);
	
	// accept/reject new link with multicanonic Metropolis step
	int acc = multicanonic_metropolis_step_single_link(GC, geo, param, r, i, old_link);
	return acc;
	}

// perform an update with heatbath in the presence of a defect with multicanonic Metropolis step
int multicanonic_heatbath_with_defect(Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
										long r, int i)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
		{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(i >= STDIM)
		{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	#endif
	
	GAUGE_GROUP stap;
	
	#ifndef THETA_MODE
	calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
	#else
	calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
	#endif
	
	// store old link before update
	GAUGE_GROUP old_link = GC->lattice[r][i];
	
	// perform usual single-link over-heat-bath update
	single_heatbath(&(GC->lattice[r][i]), &stap, param);
	
	// accept/reject new link with multicanonic Metropolis step
	int acc = multicanonic_metropolis_step_single_link(GC, geo, param, r, i, old_link);
	return acc;
	}

#endif
