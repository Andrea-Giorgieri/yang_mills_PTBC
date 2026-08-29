#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include "../include/macro.h"
#include "../include/gauge_conf.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OPENMP_MODE
#include <omp.h>
#endif
#ifdef HASH_MODE
#include <openssl/md5.h>
#endif

#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/memalign.h"
#include "../include/tens_prod.h"


void allocate_lattice_with_copy(Gauge_Conf *GC, GParam const *const param)
	{
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice), param->d_volume, __FILE__, __LINE__);
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice_copy), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice[r]), STDIM, __FILE__, __LINE__);
		allocate_array_GAUGE_GROUP(&(GC->lattice_copy[r]), STDIM, __FILE__, __LINE__);
		}
	}


void allocate_lattice_cold_with_copy(Gauge_Conf *GC, GParam const *const param)
	{
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice_cold), param->d_volume, __FILE__, __LINE__);
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice_copy_cold), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice_cold[r]), STDIM, __FILE__, __LINE__);
		allocate_array_GAUGE_GROUP(&(GC->lattice_copy_cold[r]), STDIM, __FILE__, __LINE__);
		}
	}


void allocate_C(Gauge_Conf *GC, GParam const *const param)
	{
	allocate_array_double_pointer(&(GC->C), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_double(&(GC->C[r]), STDIM, __FILE__, __LINE__);
		}
	}


void allocate_Z_with_copy(Gauge_Conf *GC, GParam const *const param)
	{
	allocate_array_double_complex_pointer(&(GC->Z), param->d_volume, __FILE__, __LINE__);
	allocate_array_double_complex_pointer(&(GC->Z_copy), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < (param->d_volume); r++)
		{
		allocate_array_double_complex(&(GC->Z[r]), param->d_n_planes + 1, __FILE__, __LINE__);
		allocate_array_double_complex(&(GC->Z_copy[r]), param->d_n_planes + 1, __FILE__, __LINE__);
		}
	}


void initialize_Z_with_copy(Gauge_Conf *GC, GParam const *const param, int x_mu, int x_nu, int x_obc)
	{
	int const si_bulk = param->d_n_planes;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		int is_on_open_boundary = 0;
		int cut_obc_dir_link = 0;
		int cartcoord[STDIM];
		si_to_cart(cartcoord, r, param);

		//for open boundary conditions: the cut is between x_obc and x_obc + 1
		GC->Z[r][si_bulk] = 1.0 + I * 0.0;
		if(param->d_obc_dir != -1)
			{
			int const dir_obc = param->d_obc_dir;
			int const r_obc = cartcoord[dir_obc];
			int const L_obc = param->d_size[dir_obc];
			int const obc_distance = link_ring_distance(r_obc, x_obc, L_obc);
			double const bulk_border_distance = (L_obc - param->d_obc_bulk) / 2.0 - obc_distance;
			if(r_obc == x_obc)
				{
				is_on_open_boundary = 1;
				cut_obc_dir_link = 1;
				}
			if(r_obc == periodic_condition(x_obc + 1, L_obc))
				{
				is_on_open_boundary = 1;
				cut_obc_dir_link = 0;
				}
			if(bulk_border_distance > 0.01)
				{
				GC->Z[r][si_bulk] = 0.5 + I * 0.0;
				}
			if(bulk_border_distance > 1.01)
				{
				GC->Z[r][si_bulk] = 0.0 + I * 0.0;
				}
			}
		GC->Z_copy[r][si_bulk] = GC->Z[r][si_bulk];

		for(int i = 0; i < STDIM; i++)
			{
			for(int j = i + 1; j < STDIM; j++)
				{
				//for anti-clockwise and clockwise plaquette
				int const si_ij = dirs_to_si(i, j);
				int const si_ji = dirs_to_si(j, i);

				// initialize to 1, 0, 1/2 if pbc, obc(temporal), obc(spatial) respectively,
				// then multiply by twist phase
				GC->Z[r][si_ij] = 1.0 + I * 0.0;
				if(i == param->d_obc_dir || j == param->d_obc_dir)
					{
					if(cut_obc_dir_link == 1)
						{
						GC->Z[r][si_ij] = 0.0 + I * 0.0;
						}
					}
				else
					{
					if(is_on_open_boundary == 1)
						{
						GC->Z[r][si_ij] = 0.5 + I * 0.0;
						}
					}
				if(cartcoord[i] == x_mu && cartcoord[j] == x_nu)
					{
					GC->Z[r][si_ij] *= cexp(I * PI2_N * (param->d_k_twist[si_ij]));
					}
				GC->Z[r][si_ji] = conj(GC->Z[r][si_ij]);
				GC->Z_copy[r][si_ij] = GC->Z[r][si_ij];
				GC->Z_copy[r][si_ji] = GC->Z[r][si_ji];
				}
			}
		}
	}


void equal_lattice(GAUGE_GROUP *const *const lattice1,
                   GAUGE_GROUP const *const *const lattice2,
                   GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		equal(&(lattice1[r][i]), &(lattice2[r][i]));
		}
	}


void equal_equal_lattice(GAUGE_GROUP *const *const lattice1,
                         GAUGE_GROUP *const *const lattice2,
                         GAUGE_GROUP const *const *const lattice3,
                         GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		equal(&(lattice1[r][i]), &(lattice3[r][i]));
		equal(&(lattice2[r][i]), &(lattice3[r][i]));
		}
	}


double lattice_total_dist(GAUGE_GROUP const *const *const lattice1,
                          GAUGE_GROUP const *const *const lattice2,
                          GParam const *const param)
	{
	double res = 0.0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+ : res)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		GAUGE_GROUP aux;

		equal(&aux, &(lattice1[r][i]));
		minus_equal(&aux, &(lattice2[r][i]));
		res += norm(&aux);
		}
	return res / ((double) NCOLOR * (double) NCOLOR);
	}


double lattice_max_dist(GAUGE_GROUP const *const *const lattice1,
                        GAUGE_GROUP const *const *const lattice2,
                        GParam const *const param)
	{
	double local_res[NTHREADS];
	for(int thread_num = 0; thread_num < NTHREADS; thread_num++)
		{
		local_res[thread_num] = 0.0;
		}

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));

		GAUGE_GROUP aux;
		equal(&aux, &(lattice1[r][i]));
		minus_equal(&aux, &(lattice2[r][i]));
		double const res_aux = norm(&aux);

		#ifdef OPENMP_MODE
		int const thread_num = omp_get_thread_num();
		#else
		int const thread_num = 0;
		#endif
		if(res_aux > local_res[thread_num])
			{
			local_res[thread_num] = res_aux;
			}
		}

	double res = 0.0;
	for(int thread_num = 0; thread_num < NTHREADS; thread_num++)
		{
		if(local_res[thread_num] > res)
			{
			res = local_res[thread_num];
			}
		}
	return res / ((double) NCOLOR * (double) NCOLOR);
	}


void equal_gauge_conf(Gauge_Conf *GC1, Gauge_Conf *GC2, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		equal(&(GC1->lattice[r][i]), &(GC2->lattice[r][i]));
		equal(&(GC1->lattice_copy[r][i]), &(GC2->lattice_copy[r][i]));
		#ifdef MULTICANONICAL_MODE
		equal(&(GC1->lattice_cold[r][i]), &(GC2->lattice_cold[r][i]));
		equal(&(GC1->lattice_copy_cold[r][i]), &(GC2->lattice_copy_cold[r][i]));
		#endif
		}
	#ifdef MULTICANONICAL_MODE
	GC1->stored_topcharge = GC2->stored_topcharge;
	#endif
	}


void accept_gauge_conf(Gauge_Conf *const GC, GParam const *const param)
	{
	#ifdef MULTICANONICAL_MODE
	GAUGE_GROUP **aux;
	aux = GC->lattice_cold;
	GC->lattice_cold = GC->lattice_copy_cold;
	GC->lattice_copy_cold = aux;
	#endif

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		// s = i * volume + r
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		unitarize(&(GC->lattice[r][i]));
		equal(&(GC->lattice_copy[r][i]), &(GC->lattice[r][i]));
		}
	}


void restore_gauge_conf(Gauge_Conf *const GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (param->d_volume); s++)
		{
		long const r = s % (param->d_volume);
		int const i = (int) ((s - r) / (param->d_volume));
		equal(&(GC->lattice[r][i]), &(GC->lattice_copy[r][i]));
		}
	}


void accept_gauge_conf_rectangle(Gauge_Conf *const GC, int const hierarc_level, Rect_Utils const *const rect_aux)
	{
	Rectangle const *rect = &(rect_aux->update_rect[hierarc_level]);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (rect->d_vol_rect); s++)
		{
		long const n = s % (rect->d_vol_rect);
		long const r = rect->rect_sites[n];
		int const i = (int) ((s - n) / (rect->d_vol_rect));
		unitarize(&(GC->lattice[r][i]));
		equal(&(GC->lattice_copy[r][i]), &(GC->lattice[r][i]));
		}

	#ifdef MULTICANONICAL_MODE
	GAUGE_GROUP **aux;
	aux = GC->lattice_cold;
	GC->lattice_cold = GC->lattice_copy_cold;
	GC->lattice_copy_cold = aux;
	/*
	rect = &(rect_aux->topcharge_rect[hierarc_level]);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s=0; s<STDIM*(rect->d_vol_rect); s++)
		{
		long const n = s % (rect->d_vol_rect);
		long const r = rect->rect_sites[n];
		int const i = (int) ( (s - n) / (rect->d_vol_rect) );
		equal(&(GC->lattice_cold[r][i]), &(GC->lattice[r][i]));
		}
	*/
	#endif
	}


void restore_gauge_conf_rectangle(Gauge_Conf *const GC, int const hierarc_level, Rect_Utils const *const rect_aux)
	{
	Rectangle const *rect;
	/*
	#ifdef MULTICANONICAL_MODE
	rect = &(rect_aux->topcharge_rect[hierarc_level]);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s=0; s<STDIM*(rect->d_vol_rect); s++)
		{
		long const n = s % (rect->d_vol_rect);
		long const r = rect->rect_sites[n];
		int const i = (int) ( (s - n) / (rect->d_vol_rect) );
		equal(&(GC->lattice_cold[r][i]), &(GC->lattice_copy[r][i]));
		}
	#endif
	*/
	rect = &(rect_aux->update_rect[hierarc_level]);
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long s = 0; s < STDIM * (rect->d_vol_rect); s++)
		{
		long const n = s % (rect->d_vol_rect);
		long const r = rect->rect_sites[n];
		int const i = (int) ((s - n) / (rect->d_vol_rect));
		equal(&(GC->lattice[r][i]), &(GC->lattice_copy[r][i]));
		}
	}


void init_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const *const param, char const *const filename)
	{
	GC->update_index = 0;

	for(int i = 0; i < STDIM; i++)
		{
		GC->translation[i] = 0;
		GC->stdim_shuffle[i] = i;
		}
	GC->parity_shuffle[0][0] = 0;
	GC->parity_shuffle[0][1] = param->d_n_even;
	GC->parity_shuffle[1][0] = param->d_n_even;
	GC->parity_shuffle[1][1] = param->d_even_volume;

	allocate_lattice_with_copy(GC, param);

	#ifdef THETA_MODE
	allocate_clover_array(GC, param);
	#endif

	// initialize lattice:
	// 0 = PBC ordered start
	// 1 = random start
	// 2 = from stored conf
	// 3 = TBC ordered start (only for one twisted plane)
	if(param->d_start == 0)
		{
		GAUGE_GROUP aux1, aux2;
		one(&aux1);
		for(long r = 0; r < (param->d_volume); r++)
			{
			for(int j = 0; j < STDIM; j++)
				{
				rand_matrix(&aux2);
				times_equal_real(&aux2, 0.001);
				plus_equal(&aux2, &aux1);
				unitarize(&aux2);
				equal_dag(&(GC->lattice[r][j]), &aux2);
				}
			}
		}
	else if(param->d_start == 1)
		{
		GAUGE_GROUP aux1;
		for(long r = 0; r < (param->d_volume); r++)
			{
			for(int j = 0; j < STDIM; j++)
				{
				rand_matrix(&aux1);
				equal(&(GC->lattice[r][j]), &aux1);
				}
			}
		}
	else if(param->d_start == 2)
		{
		read_gauge_conf_from_file_with_name(GC, param, filename);
		}
	else if(param->d_start == 3)
		{
		GAUGE_GROUP aux, Pmatrix, Qmatrix;

		int mu = 0;
		int nu = 0;
		int k_mu_nu = 0;

		#if NCOLOR > 2
		#if NCOLOR % 2 == 0
		double complex const zf= cexp(I * PI / (double) NCOLOR);
		#else
		double complex const zf = 1.0 + I * 0.0;
		#endif
		#endif

		// check if twist is non-trivial and save parameters (only last occurrence)
		int twisted_bc = 0;
		for(int i = 0; i < STDIM; i++)
			for(int j = i + 1; j < STDIM; j++)
				if(param->d_k_twist[dirs_to_si(i, j)] != 0)
					{
					twisted_bc = 1;
					mu = i;
					nu = j;
					k_mu_nu = param->d_k_twist[dirs_to_si(i, j)];
					}

		if(twisted_bc == 1)
			{
			// P matrix
			zero(&Pmatrix);
			#if NCOLOR == 1
			Pmatrix.comp = 1.0 + I * 0.0;
			#elif NCOLOR == 2
			Pmatrix.comp[1] = -1.0;
			#else
			for(int i = 0; i < (NCOLOR - 1); i++)
				Pmatrix.comp[m(i+1, i)] = conj(zf);
			Pmatrix.comp[m(0, NCOLOR - 1)] = conj(zf);
			#endif
			unitarize(&Pmatrix);

			// Q matrix
			one(&aux);
			zero(&Qmatrix);
			#if NCOLOR == 1
			Qmatrix.comp = 1.0 + I * 0.0;
			#elif NCOLOR == 2
			Qmatrix.comp[3] = -1.0;
			#else
			for(int i = 0; i < NCOLOR; i++)
				Qmatrix.comp[m(i, i)] = zf * cexp(I * PI2_N * (i + 1));
			#endif
			for(int i = 0; i < k_mu_nu; i++)
				times_equal(&aux, &Qmatrix);
			equal(&Qmatrix, &aux);
			unitarize(&Qmatrix);
			}

		// set all links to one, then overwrite links on the twisted plane
		for(long r = 0; r < (param->d_volume); r++)
			{
			int cartcoord[STDIM];
			si_to_cart(cartcoord, r, param);
			for(int i = 0; i < STDIM; i++)
				{
				one(&(GC->lattice[r][i]));
				}
			if(twisted_bc == 1)
				{
				if(cartcoord[mu] == 0)
					{
					equal(&(GC->lattice[r][mu]), &Pmatrix);
					}
				if(cartcoord[nu] == 0)
					{
					equal(&(GC->lattice[r][nu]), &Qmatrix);
					}
				}
			}
		}
	else
		{
		REQUIRE(0, "invalid value of d_start parameter: %d", param->d_start);
		}
	equal_lattice(GC->lattice_copy, (GAUGE_GROUP const * const *) GC->lattice, param);
	}


void init_gauge_conf(Gauge_Conf *GC, Geometry const *const geo, GParam const *const param)
	{
	GC->conf_label = 0;
	GC->replica_index = 0;

	init_gauge_conf_from_file_with_name(GC, param, param->d_conf_file);
	init_twist_cond_from_file_with_name(GC, param, param->d_twist_file);
	init_bound_cond(GC, param);
	#ifdef MULTICANONICAL_MODE
	init_multicanonic_gauge_conf(GC, geo, param);
	#else
	(void) geo;
	#endif
	}


int init_gauge_conf_step(Gauge_Conf *GC, GParam const *const param, long step)
	{
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	FILE *fp;

	GC->replica_index = 0;

	// gauge conf filename at step
	strcpy(name, param->d_conf_file);
	strcat(name, "_step_");
	sprintf(aux, "%ld", step);
	strcat(name, aux);

	// if conf file exists, initialize conf; otherwise return 0 (interpreted as no conf found)
	fp = fopen(name, "r");
	if(fp == NULL)
		return 0;
	fclose(fp);
	init_gauge_conf_from_file_with_name(GC, param, name);

	// twist filename at step
	strcpy(name, param->d_twist_file);
	strcat(name, "_step_");
	strcat(name, aux);

	// if twist file exists, initialize twist; otherwise return 0 (interpreted as no conf found)
	fp = fopen(name, "r");
	if(fp == NULL)
		return 0;
	fclose(fp);
	init_twist_cond_from_file_with_name(GC, param, name);

	return 1;
	}


int read_gauge_conf_step(Gauge_Conf *GC, GParam const *const param, long step)
	{
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	FILE *fp;

	// gauge conf filename at step
	strcpy(name, param->d_conf_file);
	strcat(name, "_step_");
	sprintf(aux, "%ld", step);
	strcat(name, aux);

	// if conf file exists, read conf; otherwise return 0 (interpreted as no conf found)
	fp = fopen(name, "r");
	if(fp == NULL)
		return 0;
	fclose(fp);
	read_gauge_conf_from_file_with_name(GC, param, name);

	// twist filename at step
	strcpy(name, param->d_twist_file);
	strcat(name, "_step_");
	strcat(name, aux);

	// if twist file exists, read twist; otherwise return 0 (interpreted as no conf found)
	fp = fopen(name, "r");
	if(fp == NULL)
		return 0;
	fclose(fp);
	int x_mu = 0;
	int x_nu = 0;
	int x_obc = param->d_obc_default_pos;
	read_twist_cond_from_file_with_name(&x_mu, &x_nu, &x_obc, param, name);
	initialize_Z_with_copy(GC, param, x_mu, x_nu, x_obc);

	return 1;
	}


// used to allocate all replicas in the parallel tempering
void init_gauge_conf_replica(Gauge_Conf **GC, Geometry const *const geo, GParam const *const param)
	{
	// allocate the vector to store replicas
	allocate_array_Gauge_Conf(GC, param->d_N_replica_pt, __FILE__, __LINE__);

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int i = 0; i < param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index_str[STD_STRING_LENGTH];
		strcpy(filename, param->d_conf_file); // filename = param->d_conf_file
		strcat(filename, "_replica_");
		sprintf(replica_index_str, "%d", i);
		strcat(filename, replica_index_str); // filename = param->d_conf_file + "_replica_${i}"
		init_gauge_conf_from_file_with_name(&((*GC)[i]), param, filename);

		strcpy(filename, param->d_twist_file); // filename = param->d_twist_file
		strcat(filename, "_replica_");
		strcat(filename, replica_index_str); // filename = param->d_twist_file + "_replica_${i}"
		init_twist_cond_from_file_with_name(&((*GC)[i]), param, filename);

		((*GC)[i]).conf_label = i;
		((*GC)[i]).replica_index = i;

		init_bound_cond(&((*GC)[i]), param);

		#ifdef MULTICANONICAL_MODE
		init_multicanonic_gauge_conf(&((*GC)[i]), geo, param);
		#else
		(void) geo;
		#endif
		}
	}


// initialization of the defect for a single replica
void init_bound_cond(Gauge_Conf *GC, GParam const *const param)
	{
	// allocation of C[r][j]
	allocate_C(GC, param);

	// initialization of C[r][j]
	for(long r = 0; r < param->d_volume; r++)
		for(int j = 0; j < STDIM; j++)
			{
			if(j == param->d_defect_dir && is_on_defect(r, param) == 1)
				GC->C[r][j] = param->d_pt_bound_cond_coeff[GC->replica_index];
			else
				GC->C[r][j] = 1.0;
			}
	}


// initialization of the twist factors
void init_twist_cond_from_file_with_name(Gauge_Conf *GC, GParam const *const param, char const *const filename)
	{
	//allocation of Z[r][j]
	allocate_Z_with_copy(GC, param);

	// default twist position
	int x_mu = 0;
	int x_nu = 0;

	// default open boundary position
	int x_obc = param->d_obc_default_pos;

	// update twist and open boundary position if starting from stored conf
	if(param->d_start == 2)
		{
		read_twist_cond_from_file_with_name(&x_mu, &x_nu, &x_obc, param, filename);
		}

	// assign Z on positions x_mu, x_nu,
	initialize_Z_with_copy(GC, param, x_mu, x_nu, x_obc);
	}


void read_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const *const param, char const *const filename)
	{
	FILE *fp;
	int err;
	GAUGE_GROUP matrix;
	#ifdef HASH_MODE
	char md5sum_new[2 * MD5_DIGEST_LENGTH + 1];
	char md5sum_old[2 * MD5_DIGEST_LENGTH + 1];
	#endif

	// open the configuration file in txt to read the header
	fp = fopen(filename, "r");
	REQUIRE(fp != NULL, "failed to open the configuration file %s in text mode", filename);

	int dimension;
	err = fscanf(fp, "%d", &dimension);
	REQUIRE(err == 1, "failed to read the dimension from the file %s", filename);
	REQUIRE(dimension == STDIM, "the configuration space-time dimension (%d) does not coincide with the macro STDIM (%d)", dimension, STDIM);

	for(int i = 0; i < STDIM; i++)
		{
		int tmp_i;
		err = fscanf(fp, "%d", &tmp_i);
		REQUIRE(err == 1, "failed to read the %d-th size of the configuration from the file %s", i, filename);
		REQUIRE(tmp_i == param->d_size[i], "the %d-th size of the configuration (%d) does not coincide with the size parameter (%d)", i, tmp_i, param->d_size[i]);
		}

	#ifdef HASH_MODE
	err = fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
	REQUIRE(err == 2, "failed to read the update index and md5sum from the file %s", filename);
	#else
	err = fscanf(fp, "%ld \n", &(GC->update_index));
	REQUIRE(err == 1, "failed to read the update index from the file %s", filename);
	#endif
	fclose(fp);

	// open the configuration file in binary to read the links
	fp = fopen(filename, "rb");
	REQUIRE(fp != NULL, "failed to open the configuration file %s in binary mode", filename);

	// read again the header
	err = 0;
	while(err != '\n')
		{
		err = fgetc(fp);
		}

	for(long lex = 0; lex < param->d_volume; lex++)
		{
		long const si = lex_to_si(lex, param);
		for(int mu = 0; mu < STDIM; mu++)
			{
			read_from_binary_file_bigen(fp, &matrix);
			equal(&(GC->lattice[si][mu]), &matrix);
			equal(&(GC->lattice_copy[si][mu]), &matrix);
			}
		}
	fclose(fp);

	#ifdef HASH_MODE
	// compute the new md5sum and check for consistency
	compute_md5sum_conf(md5sum_new, GC, param);
	int aux = strncmp(md5sum_old, md5sum_new, 2 * MD5_DIGEST_LENGTH + 1);
	REQUIRE(aux == 0, "the computed md5sum %s of the configuration file does not match the stored %s", md5sum_new, md5sum_old);
	#endif
	}


void read_twist_cond_from_file_with_name(int *x_mu, int *x_nu, int *x_obc, GParam const *const param, char const *const filename)
	{
	int err;
	FILE *fp = fopen(filename, "r");
	REQUIRE(fp != NULL, "failed to open the configuration file %s in text mode", filename);

	int dimension;
	err = fscanf(fp, "%d", &dimension);
	REQUIRE(err == 1, "failed to read the dimension from the twist file %s", filename);
	REQUIRE(dimension == STDIM, "the configuration space-time dimension (%d) does not coincide with the macro STDIM (%d)", dimension, STDIM);

	for(int i = 0; i < STDIM; i++)
		{
		int tmp_i;
		err = fscanf(fp, "%d", &tmp_i);
		REQUIRE(err == 1, "failed to read the %d-th size of the configuration from the twist file %s", i, filename);
		REQUIRE(tmp_i == param->d_size[i], "the %d-th size of the configuration (%d) does not coincide with the size parameter (%d)", i, tmp_i, param->d_size[i]);
		}

	err = fscanf(fp, "%*d %*d %d %d ", x_mu, x_nu);
	REQUIRE(err == 2, "failed to read the twist positions from the twist file %s", filename);

	err = fscanf(fp, "%*d %d", x_obc);
	if(err != 1)
		{
		REQUIRE(param->d_obc_dir == -1, "failed to read the OBC direction from the twist file %s", filename);
		*x_obc = param->d_obc_default_pos;
		}
	fclose(fp);
	}


void free_gauge_conf(Gauge_Conf *GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		free(GC->lattice[r]);
		free(GC->lattice_copy[r]);
		#ifdef MULTICANONICAL_MODE
		free(GC->lattice_cold[r]);
		free(GC->lattice_copy_cold[r]);
		#endif
		free(GC->Z[r]);
		free(GC->Z_copy[r]);
		}
	free(GC->lattice);
	free(GC->lattice_copy);
	#ifdef MULTICANONICAL_MODE
	free(GC->lattice_cold);
	free(GC->lattice_copy_cold);
	#endif
	free(GC->Z);
	free(GC->Z_copy);
	#ifdef THETA_MODE
	free_clover_array(GC, param);
	#endif
	}


void free_replica(Gauge_Conf *GC, GParam const *const param)
	{
	for(int i = 0; i < param->d_N_replica_pt; i++)
		{
		free_gauge_conf(&(GC[i]), param);
		free_bound_cond(&(GC[i]), param);
		}
	free(GC);
	}


void free_bound_cond(Gauge_Conf *GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		free(GC->C[r]);
		}
	free(GC->C);
	}


void free_twist_cond(Gauge_Conf *GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < param->d_volume; r++)
		{
		free(GC->Z[r]);
		free(GC->Z_copy[r]);
		}
	free(GC->Z);
	free(GC->Z_copy);
	}


// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Gauge_Conf const *const GC,
                                  GParam const *const param,
                                  char const *const filename)
	{
	#ifdef HASH_MODE
	char md5sum[2 * MD5_DIGEST_LENGTH + 1];
	#else
	char md5sum[2 * STD_STRING_LENGTH + 1] = {0};
	#endif

	#ifdef HASH_MODE
	compute_md5sum_conf(md5sum, GC, param);
	#endif

	FILE *fp;

	// open the configuration file in text mode to write the header
	fp = fopen(filename, "w");
	REQUIRE(fp != NULL, "failed to open the configuration file %s in text mode", filename);
	fprintf(fp, "%d ", STDIM);
	for(int i = 0; i < STDIM; i++)
		{
		fprintf(fp, "%d ", param->d_size[i]);
		}
	fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
	fclose(fp);

	// open the configuration file in binary mode to write the links
	fp = fopen(filename, "ab");
	REQUIRE(fp != NULL, "failed to open the configuration file %s in binary mode", filename);
	for(long lex = 0; lex < param->d_volume; lex++)
		{
		long const si = lex_to_si(lex, param);
		for(int mu = 0; mu < STDIM; mu++)
			{
			print_on_binary_file_bigen(fp, &(GC->lattice[si][mu]));
			}
		}
	fclose(fp);
	}


void write_twist_on_file_with_name(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   char const *const filename)
	{
	// open the twist configuration file
	FILE *fp = fopen(filename, "w");
	REQUIRE(fp != NULL, "failed to open the twist file %s in text mode", filename);
	fprintf(fp, "%d ", STDIM);
	for(int i = 0; i < STDIM; i++)
		{
		fprintf(fp, "%d ", param->d_size[i]);
		}
	fprintf(fp, "\n");

	// check if twist non-trivial, save its plane (mu,nu)
	int twisted_bc = 0, mu = 0, nu = 1, si_munu = 0;
	int cartcoord[STDIM];
	for(int i = 0; i < STDIM; i++)
		{
		cartcoord[i] = 0; // initialize to zero for later loop
		for(int j = i + 1; j < STDIM; j++)
			if(param->d_k_twist[dirs_to_si(i, j)] != 0)
				{
				twisted_bc = 1;
				mu = i;
				nu = j;
				si_munu = dirs_to_si(mu, nu);
				}
		}

	if(twisted_bc == 1) // find twist cartcoord on plane cartcoord[i] = 0 for i != mu, nu (no need to read all volume)
		{
		for(int i = 0; i < param->d_size[mu]; i++)
			{
			cartcoord[mu] = i;
			for(int j = 0; j < param->d_size[nu]; j++)
				{
				cartcoord[nu] = j;
				double complex const z = GC->Z[cart_to_si(cartcoord, param)][si_munu];
				if(fabs(carg(z)) > MIN_VALUE)
					{
					fprintf(fp, "%d %d %d %d \n", mu, nu, i, j);
					}
				}
			}
		}
	else // arbitrary twist position to avoid errors
		{
		fprintf(fp, "%d %d %d %d \n", 0, 1, 0, 0);
		}

	int j = param->d_obc_default_pos; // arbitrary boundary position to avoid errors
	if(param->d_obc_dir != -1)        // find position of open boundary
		{
		for(int i = 0; i < STDIM; i++)
			{
			cartcoord[i] = 0;
			if(i != param->d_obc_dir) mu = i; // any direction different from obc_dir
			}
		int const si_obcmu = dirs_to_si(param->d_obc_dir, mu);
		for(int i = 0; i < param->d_size[param->d_obc_dir]; i++)
			{
			cartcoord[param->d_obc_dir] = i;
			double complex const z = GC->Z[cart_to_si(cartcoord, param)][si_obcmu];
			if(cabs(z) < MIN_VALUE)
				{
				j = i;
				break;
				}
			}
		}
	fprintf(fp, "%d %d \n", param->d_obc_dir, j);
	fclose(fp);
	}


void write_conf_on_file(Gauge_Conf const *const GC, GParam const *const param)
	{
	write_conf_on_file_with_name(GC, param, param->d_conf_file);
	write_twist_on_file_with_name(GC, param, param->d_twist_file);
	}


void write_replica_on_file(Gauge_Conf const *const GC, GParam const *const param)
	{
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int i = 0; i < param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index_str[STD_STRING_LENGTH];
		strcpy(filename, param->d_conf_file);
		strcat(filename, "_replica_");
		sprintf(replica_index_str, "%d", i);
		strcat(filename, replica_index_str); // filename = d_conf_file + "_replica_${i}"
		write_conf_on_file_with_name(&(GC[i]), param, filename);

		strcpy(filename, param->d_twist_file);
		strcat(filename, "_replica_");
		strcat(filename, replica_index_str); // filename = d_twist_file + "_replica_${i}"
		write_twist_on_file_with_name(&(GC[i]), param, filename);
		}
	}


void write_replica_on_file_back(Gauge_Conf const *const GC, GParam const *const param)
	{
	static int counter = 0;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(int i = 0; i < param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index_str[STD_STRING_LENGTH], aux_back[STD_STRING_LENGTH];
		if(counter == 0)
			sprintf(aux_back, "_back0");
		else
			sprintf(aux_back, "_back1");
		strcpy(filename, param->d_conf_file);
		strcat(filename, "_replica_");
		sprintf(replica_index_str, "%d", i);
		strcat(filename, replica_index_str);
		strcat(filename, aux_back); // filename = d_conf_file + "_replica_${i}_back${counter}"
		write_conf_on_file_with_name(&(GC[i]), param, filename);

		strcpy(filename, param->d_twist_file);
		strcat(filename, "_replica_");
		strcat(filename, replica_index_str);
		strcat(filename, aux_back); // filename = d_conf_file + "_replica_${i}_back${counter}"
		write_twist_on_file_with_name(&(GC[i]), param, filename);
		}
	counter = 1 - counter;
	}


void write_conf_on_file_back(Gauge_Conf const *const GC, GParam const *const param)
	{
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	static int counter = 0;

	if(counter == 0)
		{
		sprintf(aux, "_back0");
		}
	else
		{
		sprintf(aux, "_back1");
		}

	strcpy(name, param->d_conf_file);
	strcat(name, aux);
	write_conf_on_file_with_name(GC, param, name);

	strcpy(name, param->d_twist_file);
	strcat(name, aux);
	write_twist_on_file_with_name(GC, param, name);

	counter = 1 - counter;
	}


// allocate GC and initialize with GC2, including the twist factors
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const *const GC2, GParam const *const param)
	{
	GC->update_index = GC2->update_index;
	GC->replica_index = GC2->replica_index;

	for(int i = 0; i < STDIM; i++)
		{
		GC->translation[i] = GC2->translation[i];
		GC->stdim_shuffle[i] = i;
		}
	GC->parity_shuffle[0][0] = 0;
	GC->parity_shuffle[0][1] = param->d_n_even;
	GC->parity_shuffle[1][0] = param->d_n_even;
	GC->parity_shuffle[1][1] = param->d_even_volume;

	allocate_lattice_with_copy(GC, param);

	#ifdef MULTICANONICAL_MODE
	allocate_lattice_cold_with_copy(GC, param);
	#endif

	allocate_Z_with_copy(GC, param);

	#ifdef THETA_MODE
	allocate_clover_array(GC, param);
	#endif

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for(long r = 0; r < (param->d_volume); r++)
		{
		for(int j = 0; j < STDIM; j++)
			{
			equal(&(GC->lattice[r][j]), &(GC2->lattice[r][j]));
			equal(&(GC->lattice_copy[r][j]), &(GC2->lattice_copy[r][j]));
			}
		for(int j = 0; j < param->d_n_planes + 1; j++)
			{
			GC->Z[r][j] = GC2->Z[r][j];
			GC->Z_copy[r][j] = GC2->Z_copy[r][j];
			}
		}
	}


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Gauge_Conf const *const GC, GParam const *const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];
	GAUGE_GROUP matrix;

	MD5_Init(&mdContext);
	for(long lex = 0; lex < param->d_volume; lex++)
		{
		long const si = lex_to_si(lex, param);
		for(int mu = 0; mu < STDIM; mu++)
			{
			equal(&matrix, &(GC->lattice[si][mu]));

	#if NCOLOR == 1
	MD5_Update(&mdContext, &(matrix.comp), sizeof(double complex));
	#elif NCOLOR == 2
	for(int k = 0; k < 4; k++)
		{
		MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
		}
	#else
	for(int k = 0; k < NCOLOR * NCOLOR; k++)
		{
		MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double complex)
		)
		;
		}
	#endif
			}
		}
	MD5_Final(c, &mdContext);

	for(int k = 0; k < MD5_DIGEST_LENGTH; k++)
		{
		sprintf(&(res[2 * k]), "%02x", c[k]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


// allocate the ml_polycorr arrays and related stuff
void allocate_polycorr_stuff(Gauge_Conf *GC, GParam const *const param)
	{
	allocate_array_TensProd_pointer_pointer(&(GC->ml_polycorr), NLEVELS, __FILE__, __LINE__);
	for(int i = 0; i < NLEVELS; i++)
		{
		allocate_array_TensProd_pointer(&(GC->ml_polycorr[i]), param->d_size[0] / param->d_ml_step[i], __FILE__, __LINE__);
		for(int j = 0; j < (param->d_size[0] / param->d_ml_step[i]); j++)
			{
			allocate_array_TensProd(&(GC->ml_polycorr[i][j]), param->d_space_vol[0], __FILE__, __LINE__);
			}
		}
	allocate_array_GAUGE_GROUP_pointer(&(GC->loc_poly), param->d_size[0] / param->d_ml_step[NLEVELS - 1], __FILE__, __LINE__);
	for(int i = 0; i < param->d_size[0] / param->d_ml_step[NLEVELS - 1]; i++)
		allocate_array_GAUGE_GROUP(&(GC->loc_poly[i]), param->d_space_vol[0], __FILE__, __LINE__);
	}


// free the ml_polycorr arrays and related stuff
void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const *const param)
	{
	for(int i = 0; i < NLEVELS; i++)
		{
		for(int j = 0; j < (param->d_size[0] / param->d_ml_step[i]); j++)
			{
			free(GC->ml_polycorr[i][j]);
			}
		free(GC->ml_polycorr[i]);
		}
	free(GC->ml_polycorr);
	for(int i = 0; i < param->d_size[0] / param->d_ml_step[NLEVELS - 1]; i++)
		{
		free(GC->loc_poly[i]);
		}
	free(GC->loc_poly);
	}


// allocate the ml_polycorr, polyplaq arrays and related stuff
void allocate_tube_disc_stuff(Gauge_Conf *GC,
                              GParam const *const param)
	{
	allocate_polycorr_stuff(GC, param);

	allocate_array_TensProd_pointer(&(GC->ml_polyplaq), NLEVELS, __FILE__, __LINE__);
	for(int i = 0; i < NLEVELS; i++) allocate_array_TensProd(&(GC->ml_polyplaq[i]), param->d_space_vol[0], __FILE__, __LINE__);

	allocate_array_double_complex(&(GC->loc_plaq), param->d_space_vol[0], __FILE__, __LINE__);
	}


// free the ml_polycorr, ml_polyplaq arrays and relates stuff
void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const *const param)
	{
	free_polycorr_stuff(GC, param);
	for(int i = 0; i < NLEVELS; i++) free(GC->ml_polyplaq[i]);
	free(GC->ml_polyplaq);
	free(GC->loc_plaq);
	}


// allocate the polycorr, polyplaq, polyplaqconn arrays and related stuff
void allocate_tube_conn_stuff(Gauge_Conf *GC,
                              GParam const *const param)
	{
	allocate_tube_disc_stuff(GC, param);
	allocate_array_TensProd_pointer(&(GC->ml_polyplaqconn), NLEVELS, __FILE__, __LINE__);
	for(int i = 0; i < NLEVELS; i++) allocate_array_TensProd(&(GC->ml_polyplaqconn[i]), param->d_space_vol[0], __FILE__, __LINE__);
	allocate_array_GAUGE_GROUP(&(GC->loc_plaqconn), param->d_space_vol[0], __FILE__, __LINE__);
	}


// free the polycorr, polyplaq, polyplaqconn arrays and related stuff
void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const *const param)
	{
	free_tube_disc_stuff(GC, param);
	for(int i = 0; i < NLEVELS; i++) free(GC->ml_polyplaqconn[i]);
	free(GC->ml_polyplaqconn);
	free(GC->loc_plaqconn);
	}


void allocate_multilevel_arrays(Gauge_Conf *GC, GParam const *const param, Multilevel_Obs const ml_obs)
	{
	switch(ml_obs)
		{
		case NONE:
			break;
		case POLYCORR:
		case POLYCORR_LONG:
			allocate_polycorr_stuff(GC, param);
			break;
		case TUBE_DISC:
			allocate_tube_disc_stuff(GC, param);
			break;
		case TUBE_CONN:
		case TUBE_CONN_LONG:
			allocate_tube_conn_stuff(GC, param);
			break;
		default:
			REQUIRE(0, "unknown multilevel observable (%d)\n", (int)ml_obs);
		}
	}


void free_multilevel_arrays(Gauge_Conf *GC, GParam const *const param, Multilevel_Obs const ml_obs)
	{
	switch(ml_obs)
		{
		case NONE:
			break;
		case POLYCORR:
		case POLYCORR_LONG:
			free_polycorr_stuff(GC, param);
			break;
		case TUBE_DISC:
			free_tube_disc_stuff(GC, param);
			break;
		case TUBE_CONN:
		case TUBE_CONN_LONG:
			free_tube_conn_stuff(GC, param);
			break;
		default:
			REQUIRE(0, "unknown multilevel observable (%d)\n", (int)ml_obs);
		}
	}


// save ml_polycorr[0] arrays on file
void write_polycorr_on_file(Gauge_Conf const *const GC,
                            GParam const *const param,
                            int const iteration)
	{
	FILE *fp;
	#ifdef HASH_MODE
	char md5sum[2 * MD5_DIGEST_LENGTH + 1];
	compute_md5sum_polycorr(md5sum, GC, param);
	#else
	char md5sum[2 * STD_STRING_LENGTH + 1] = {0};
	#endif

	// open the configuration file in text mode
	fp = fopen(param->d_ml_file, "w");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);
	fprintf(fp, "%ld %d %s\n", param->d_space_vol[0], iteration, md5sum);
	fclose(fp);

	// open the configuration file in binary mode
	fp = fopen(param->d_ml_file, "ab");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);
	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long i = 0; i < (param->d_space_vol[0]); i++)
			{
			print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}
	fclose(fp);
	}


// read ml_polycorr[0] arrays from file
void read_polycorr_stuff_from_file(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   int *iteration)
	{
	long loc_space_vol;
	FILE *fp;
	#ifdef HASH_MODE
	char md5sum_new[2 * MD5_DIGEST_LENGTH + 1];
	char md5sum_old[2 * MD5_DIGEST_LENGTH + 1];
	#endif

	// open the multilevel file in text mode to read the header
	fp = fopen(param->d_ml_file, "r");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);

	#ifdef HASH_MODE
	int err = fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
	REQUIRE(err == 3, "failed to read the header of the multilevel file %s", param->d_ml_file);
	#else
	int err = fscanf(fp, "%ld %d \n", &loc_space_vol, iteration);
	REQUIRE(err == 2, "failed to read the header of the multilevel file %s", param->d_ml_file);
	#endif
	REQUIRE(loc_space_vol == param->d_space_vol[0], "space volume in multilevel file %s does not coincide with the input parameters", param->d_ml_file);
	fclose(fp);

	// open the multilevel file in binary mode to read the configuration
	fp = fopen(param->d_ml_file, "rb");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);

	// read again the header: loc_space_vol, iteration, hash
	int i = 0;
	while(i != '\n')
		{
		i = fgetc(fp);
		}

	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long r = 0; r < (param->d_space_vol[0]); r++)
			{
			read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][r]));
			}
		}
	fclose(fp);

	#ifdef HASH_MODE
	// compute the new md5sum and check for consistency
	compute_md5sum_polycorr(md5sum_new, GC, param);
	int aux = strncmp(md5sum_old, md5sum_new, 2 * MD5_DIGEST_LENGTH + 1);
	REQUIRE(aux == 0, "the computed md5sum %s of the multilevel file does not match the stored %s", md5sum_new, md5sum_old);
	#endif
	}


// compute the md5sum of the ml_polycorr[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorr(char *res, Gauge_Conf const *const GC, GParam const *const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];

	MD5_Init(&mdContext);
	size_t const size = sizeof(double complex);
	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long i = 0; i < (param->d_space_vol[0]); i++)
			{
			for(int n1 = 0; n1 < NCOLOR; n1++)
				{
				for(int n2 = 0; n2 < NCOLOR; n2++)
					{
					for(int n3 = 0; n3 < NCOLOR; n3++)
						{
						for(int n4 = 0; n4 < NCOLOR; n4++)
							{
							MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), size);
							}
						}
					}
				}
			}
		}

	MD5_Final(c, &mdContext);

	for(long i = 0; i < MD5_DIGEST_LENGTH; i++)
		{
		sprintf(&(res[2 * i]), "%02x", c[i]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


// save ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays on file
void write_tube_conn_stuff_on_file(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   int const iteration)
	{
	#ifdef HASH_MODE
	char md5sum[2 * MD5_DIGEST_LENGTH + 1];
	compute_md5sum_tube_conn_stuff(md5sum, GC, param);
	#else
	char md5sum[2 * STD_STRING_LENGTH + 1] = {0};
	#endif
	FILE *fp;

	// open the configuration file in text mode
	fp = fopen(param->d_ml_file, "w");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);
	fprintf(fp, "%ld %d %s\n", param->d_space_vol[0], iteration, md5sum);
	fclose(fp);

	// open the configuration file in binary mode
	fp = fopen(param->d_ml_file, "ab");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);
	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long i = 0; i < (param->d_space_vol[0]); i++)
			{
			print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}
	for(long i = 0; i < (param->d_space_vol[0]); i++)
		{
		print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
		}
	for(long i = 0; i < (param->d_space_vol[0]); i++)
		{
		print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
		}
	fclose(fp);
	}


// read ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays from file
void read_tube_conn_stuff_from_file(Gauge_Conf const *const GC,
                                    GParam const *const param,
                                    int *iteration)
	{
	FILE *fp;
	#ifdef HASH_MODE
	char md5sum_new[2 * MD5_DIGEST_LENGTH + 1];
	char md5sum_old[2 * MD5_DIGEST_LENGTH + 1];
	#endif

	// open the multilevel file in text mode
	fp = fopen(param->d_ml_file, "r");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);

	long loc_space_vol;
	#ifdef HASH_MODE
	int err = fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
	REQUIRE(err == 3, "failed to read the header of the multilevel file %s", param->d_ml_file);
	#else
	int err = fscanf(fp, "%ld %d \n", &loc_space_vol, iteration);
	REQUIRE(err == 2, "failed to read the header of the multilevel file %s", param->d_ml_file);
	#endif
	REQUIRE(loc_space_vol == param->d_space_vol[0], "space volume in the multilevel file %s does not coincide with input parameters", param->d_ml_file);
	fclose(fp);

	// open the multilevel file in binary mode
	fp = fopen(param->d_ml_file, "rb");
	REQUIRE(fp != NULL, "failed to open the file %s", param->d_ml_file);

	// read again the header: loc_space_vol, iteration, hash
	int i = 0;
	while(i != '\n')
		{
		i = fgetc(fp);
		}

	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long r = 0; r < (param->d_space_vol[0]); r++)
			{
			read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][r]));
			}
		}
	for(long r = 0; r < (param->d_space_vol[0]); r++)
		{
		read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][r]));
		}
	for(long r = 0; r < (param->d_space_vol[0]); r++)
		{
		read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][r]));
		}

	fclose(fp);

	#ifdef HASH_MODE
	// compute the new md5sum and check for consistency
	compute_md5sum_tube_conn_stuff(md5sum_new, GC, param);
	int aux = strncmp(md5sum_old, md5sum_new, 2 * MD5_DIGEST_LENGTH + 1);
	REQUIRE(aux == 0, "the computed md5sum %s of the multilevel file does not match the stored %s", md5sum_new, md5sum_old);
	#endif
	}


// compute the md5sum of the ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const *const GC, GParam const *const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];

	MD5_Init(&mdContext);

	size_t const size = sizeof(double complex);

	for(int j = 0; j < param->d_size[0] / param->d_ml_step[0]; j++)
		{
		for(long i = 0; i < (param->d_space_vol[0]); i++)
			{
			for(int n1 = 0; n1 < NCOLOR; n1++)
				{
				for(int n2 = 0; n2 < NCOLOR; n2++)
					{
					for(int n3 = 0; n3 < NCOLOR; n3++)
						{
						for(int n4 = 0; n4 < NCOLOR; n4++)
							{
							MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), size);
							}
						}
					}
				}
			}
		}

	for(long i = 0; i < (param->d_space_vol[0]); i++)
		{
		for(int n1 = 0; n1 < NCOLOR; n1++)
			{
			for(int n2 = 0; n2 < NCOLOR; n2++)
				{
				for(int n3 = 0; n3 < NCOLOR; n3++)
					{
					for(int n4 = 0; n4 < NCOLOR; n4++)
						{
						MD5_Update(&mdContext, &((GC->ml_polyplaq[0][i]).comp[n1][n2][n3][n4]), size);
						}
					}
				}
			}
		}

	for(long i = 0; i < (param->d_space_vol[0]); i++)
		{
		for(int n1 = 0; n1 < NCOLOR; n1++)
			{
			for(int n2 = 0; n2 < NCOLOR; n2++)
				{
				for(int n3 = 0; n3 < NCOLOR; n3++)
					{
					for(int n4 = 0; n4 < NCOLOR; n4++)
						{
						MD5_Update(&mdContext, &((GC->ml_polyplaqconn[0][i]).comp[n1][n2][n3][n4]), size);
						}
					}
				}
			}
		}

	MD5_Final(c, &mdContext);

	for(long i = 0; i < MD5_DIGEST_LENGTH; i++)
		{
		sprintf(&(res[2 * i]), "%02x", c[i]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


void write_multilevel_status_on_file(Gauge_Conf const *const GC,
                                     GParam const *const param,
                                     int const iteration,
                                     Multilevel_Obs const ml_obs)
	{
	switch(ml_obs)
		{
		case NONE:
			break;
		case POLYCORR_LONG:
			write_polycorr_on_file(GC, param, iteration);
			break;
		case TUBE_CONN_LONG:
			write_tube_conn_stuff_on_file(GC, param, iteration);
			break;
		default:
			REQUIRE(0, "unknown multilevel observable (%d)\n", (int)ml_obs);
		}
	}


void read_multilevel_status_from_file(Gauge_Conf const *const GC,
                                      GParam const *const param,
                                      int *iteration,
                                      Multilevel_Obs const ml_obs)
	{
	switch(ml_obs)
		{
		case NONE:
			break;
		case POLYCORR_LONG:
			read_polycorr_stuff_from_file(GC, param, iteration);
			break;
		case TUBE_CONN_LONG:
			read_tube_conn_stuff_from_file(GC, param, iteration);
			break;
		default:
			REQUIRE(0, "unknown multilevel observable (%d)\n", (int)ml_obs);
		}
	}


// allocate the clovers arrays
void allocate_clover_array(Gauge_Conf *GC,
                           GParam const *const param)
	{
	allocate_array_GAUGE_GROUP_pointer_pointer(&(GC->clover_array), param->d_volume, __FILE__, __LINE__);
	for(long r = 0; r < param->d_volume; r++)
		{
		allocate_array_GAUGE_GROUP_pointer(&(GC->clover_array[r]), STDIM, __FILE__, __LINE__);
		for(int i = 0; i < STDIM; i++) allocate_array_GAUGE_GROUP(&(GC->clover_array[r][i]), STDIM, __FILE__, __LINE__);
		}
	}


// free the clovers arrays
void free_clover_array(Gauge_Conf *GC,
                       GParam const *const param)
	{
	for(long r = 0; r < param->d_volume; r++)
		{
		for(int i = 0; i < STDIM; i++)
			{
			free(GC->clover_array[r][i]);
			}
		free(GC->clover_array[r]);
		}
	free(GC->clover_array);
	}

#endif
