#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include"../include/macro.h"

#ifdef HASH_MODE
#include<openssl/md5.h>
#endif
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"../include/function_pointers.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/memalign.h"

void equal_gauge_conf(Gauge_Conf *GC1, Gauge_Conf *GC2, GParam const * const param)
	{
	long r;
	for(int i=0; i<STDIM; i++)
		{
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(r)
		#endif
		for(r=0; r<param->d_volume; r++) equal(&(GC1->lattice[r][i]), &(GC2->lattice[r][i]));
		}
	}

void init_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const * const param, char const * const filename)
	{
	long r, j;
	
	// allocate the local lattice
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice), param->d_volume, __FILE__, __LINE__);
	for(r=0; r<(param->d_volume); r++) allocate_array_GAUGE_GROUP(&(GC->lattice[r]), STDIM, __FILE__, __LINE__);
	
	#ifdef THETA_MODE
	alloc_clover_array(GC, param);
	#endif
	
	// initialize lattice
	if(param->d_start==0) // ordered start
		{
		GAUGE_GROUP aux1, aux2;
		one(&aux1);
		
		GC->update_index=0;
		
		for(r=0; r<(param->d_volume); r++)
			{
			for(j=0; j<STDIM; j++)
				{
				rand_matrix(&aux2);
				times_equal_real(&aux2, 0.001);
				plus_equal(&aux2, &aux1);
				unitarize(&aux2);
				equal_dag(&(GC->lattice[r][j]), &aux2);
				}
			}
		}
	if(param->d_start==1)	// random start
		{
		GAUGE_GROUP aux1;
		
		GC->update_index=0;
		
		for(r=0; r<(param->d_volume); r++)
			{
			for(j=0; j<STDIM; j++)
				{
				rand_matrix(&aux1);
				equal(&(GC->lattice[r][j]), &aux1);
				}
			}
		}
	
	if(param->d_start==2) // initialize from stored conf
		{
		read_gauge_conf_from_file_with_name(GC, param, filename);
		}
	
	if(param->d_start==3) // cold twisted start (only for one twist parameter)
		{
		int i, j, k_mu_nu, mu, nu, twisted_bc, cartcoord[STDIM];
		double complex zf;
		GAUGE_GROUP aux, Pmatrix, Qmatrix;
		
		// to suppress gcc warning "maybe-uninitialized"
		mu = 0;
		nu = 0;
		
		if ((NCOLOR%2)==0) zf = cexp(I*PI/(double)NCOLOR);
		else zf = 1.0 + 0.0*I;
		
		twisted_bc = 0;	// check if twist is non trivial and save parameters (only last occurence)
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
				if (param->d_k_twist[dirs_to_si(i,j)] != 0)
					{
					twisted_bc = 1;
					mu = i;
					nu = j;
					k_mu_nu = param->d_k_twist[dirs_to_si(i,j)];
					}
		
		if (twisted_bc == 1)
			{
			// P matrix
			zero(&Pmatrix);
			for(i=0; i<(NCOLOR-1); i++) Pmatrix.comp[m(i+1,i)] = conj(zf);
			Pmatrix.comp[m(0,2)] = conj(zf);
			unitarize(&Pmatrix);
		
			// Q matrix
			one(&aux);
			zero(&Qmatrix);
			for(i=0; i<NCOLOR; i++) Qmatrix.comp[m(i,i)] = zf*cexp(PI2*I*(i+1)/(double)NCOLOR);
			for(i=0; i<k_mu_nu; i++) times_equal(&aux,&Qmatrix);
			equal(&Qmatrix,&aux);
			unitarize(&Qmatrix);
			}
		
		GC->update_index=0;
		for(r=0; r<(param->d_volume); r++)
			{
			si_to_cart(cartcoord, r, param);
			for(i=0; i<STDIM; i++) one(&(GC->lattice[r][i])); // all links set to one
			if (twisted_bc == 1) //overwrite links on the twisted plane
				{
				if (cartcoord[mu] == 0) equal(&(GC->lattice[r][mu]), &Pmatrix);
				if (cartcoord[nu] == 0) equal(&(GC->lattice[r][nu]), &Qmatrix);
				}
			}
		}
	}

void init_gauge_conf(Gauge_Conf *GC, GParam const * const param)
	{
	init_gauge_conf_from_file_with_name(GC, param, param->d_conf_file);
	init_twist_cond_from_file_with_name(GC, param, param->d_twist_file);
	}

void init_gauge_conf_step(Gauge_Conf *GC, GParam const * const param, long step, int *stop)
	{
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	FILE *file;
	
	//gauge conf filename at step
	strcpy(name, param->d_conf_file);
	strcat(name, "_step_");
	sprintf(aux, "%ld", step);
	strcat(name, aux);
	
	//check if conf file exists and initialize conf if it does
	if ((file = fopen(name, "r")))
		{
		fclose(file);
		init_gauge_conf_from_file_with_name(GC, param, name);
		}
	else *stop = 1;
	
	//twist filename at step
	strcpy(name, param->d_twist_file);
	strcat(name, "_step_");
	strcat(name, aux);
	
	//check if twist file exists and initialize cond if it does
	if ((file = fopen(name, "r")))
		{
		fclose(file);
		init_twist_cond_from_file_with_name(GC, param, name);
		}
	else *stop = 1;
	}

// used to allocate all replicas in the parallel tempering
void init_gauge_conf_replica(Gauge_Conf **GC, GParam const * const param)
	{
	int i;
	
	// allocate the vector to store replicas
	allocate_array_Gauge_Conf(GC, param->d_N_replica_pt, __FILE__, __LINE__);
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif
	for(i=0; i<param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH];
		strcpy(filename,param->d_conf_file); // filename = param->d_conf_file
		strcat(filename,"_replica_");			
		sprintf(replica_index, "%d", i);
		strcat(filename, replica_index); // filename = param->d_conf_file + "_replica_${i}"
		init_gauge_conf_from_file_with_name(&((*GC)[i]), param, filename);
		
		strcpy(filename,param->d_twist_file); // filename = param->d_twist_file
		strcat(filename,"_replica_");			
		strcat(filename, replica_index); // filename = param->d_twist_file + "_replica_${i}"
		init_twist_cond_from_file_with_name(&((*GC)[i]), param, filename);
		
		init_bound_cond(&((*GC)[i]), param, i);
		
		((*GC)[i]).conf_label=i;		
		}
	}

// initialization of the defect for the replica with label i
void init_bound_cond(Gauge_Conf *GC, GParam const * const param, int const i)
	{
	const int TRUE	= 1;
	const int FALSE = 0;
	long r;
	int j, is_on_defect;
	int cartcoord[STDIM];
	
	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };
	
	//allocation of C[r][j]
	allocate_array_double_pointer(&(GC->C), param->d_volume, __FILE__, __LINE__);
	for(r=0; r<(param->d_volume); r++) allocate_array_double(&(GC->C[r]), STDIM, __FILE__, __LINE__);
	
	// initialization of C[r][j]
	
	// start initializing them to 1
	for(r=0; r<param->d_volume; r++)
		for(j=0; j<STDIM; j++)
			{
			GC->C[r][j]=1.0;
			}
	
	// if there is more than 1 replica, initialize c(r) != 1 on the defect and along direction <defect_dir>
	if (param->d_N_replica_pt > 1)
		{
		j=param->d_defect_dir;
		for(r=0; r<param->d_volume; r++)
			{
			// check if r is on the defect or not
			si_to_cart(cartcoord, r, param);
			if( (cartcoord[param->d_defect_dir]==((param->d_size[param->d_defect_dir])-1)) &&
				(cartcoord[perp_dir[param->d_defect_dir][0]]<(param->d_L_defect[0])) &&
				(cartcoord[perp_dir[param->d_defect_dir][1]]<(param->d_L_defect[1])) &&
				(cartcoord[perp_dir[param->d_defect_dir][2]]<(param->d_L_defect[2])) ) is_on_defect=TRUE;
			else is_on_defect=FALSE;
			
			// if r is on defect assign bound cond
			if (is_on_defect==TRUE) {GC->C[r][j]=param->d_pt_bound_cond_coeff[i];}
			}
		}
	}

// initialization of the twist factors
void init_twist_cond_from_file_with_name(Gauge_Conf *GC, GParam const * const param, char const * const filename)
	{
	long r;
	int i, j, x_mu, x_nu;
	int cartcoord[STDIM];
	
	//allocation of Z[r][j]
	allocate_array_double_complex_pointer(&(GC->Z), param->d_volume, __FILE__, __LINE__);
	for(r=0; r<(param->d_volume); r++) allocate_array_double_complex(&(GC->Z[r]), param->d_n_planes, __FILE__, __LINE__);
	
	// initialization of Z[r][j]
	
	// start initializing them to 1
	for(r=0; r<param->d_volume; r++)
		for(j=0; j<param->d_n_planes; j++)
				GC->Z[r][j]=1.0+0.0*I;
	
	x_mu = 0;
	x_nu = 0;
	
	if(param->d_start==2) // initialize from stored conf
		{
		read_twist_cond_from_file_with_name(&x_mu, &x_nu, param, filename);
		}
	
	// assign Z on positions x_mu, x_nu, 
	for(r=0; r<param->d_volume; r++)
		{
		si_to_cart(cartcoord, r, param);
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
				if (cartcoord[i] == x_mu && cartcoord[j] == x_nu)
					{
					GC->Z[r][dirs_to_si(i,j)] = cexp(I*PI2*(param->d_k_twist[dirs_to_si(i,j)])/(double)NCOLOR);	//for clockwise plaquette
					GC->Z[r][dirs_to_si(j,i)] = conj(GC->Z[r][dirs_to_si(i,j)]);								//for anticlockwise plaquette
					}
		}
}

void read_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const * const param, char const * const filename)
	{
	FILE *fp;
	int i, dimension, tmp_i;
	int err, mu;
	long lex, si;
	GAUGE_GROUP matrix;
	#ifdef HASH_MODE
	char md5sum_new[2*MD5_DIGEST_LENGTH+1];
	char md5sum_old[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum_old[2*STD_STRING_LENGTH+1]={0};
	#endif
	
	fp=fopen(filename, "r"); // open the configuration file
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else // read the txt header of the configuration
		{
		err=fscanf(fp, "%d", &dimension);
		if(err!=1)
			{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(dimension != STDIM)
			{
			fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n", dimension, STDIM);
			exit(EXIT_FAILURE);
			}
		
		for(i=0; i<STDIM; i++)
			{
			err=fscanf(fp, "%d", &tmp_i);
			if(err!=1)
				{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			if(tmp_i != param->d_size[i])
				{
				fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
				exit(EXIT_FAILURE);
				}
			}
		
		err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
		if(err!=2)
			{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		
		fclose(fp);
		}
	
	fp=fopen(filename, "rb"); // open the configuration file in binary
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		// read again the header
		err=0;
		while(err!='\n')
				{
				err=fgetc(fp);
				}
		
		for(lex=0; lex<param->d_volume; lex++)
			{
			si=lex_to_si(lex, param);
			for(mu=0; mu<STDIM; mu++)
				{
				read_from_binary_file_bigen(fp, &matrix);
				equal(&(GC->lattice[si][mu]), &matrix);
				}
			}
		fclose(fp);
		
		#ifdef HASH_MODE
		// compute the new md5sum and check for consistency
		compute_md5sum_conf(md5sum_new, GC, param);
		if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
			{
			fprintf(stderr, "The computed md5sum %s does not match the stored %s (%s, %d)\n", md5sum_new, md5sum_old, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		#endif
		}
	}

void read_twist_cond_from_file_with_name(int *x_mu, int *x_nu, GParam const * const param, char const * const filename)
	{
	FILE *fp;
	int err, i, tmp_i, dimension;

	fp=fopen(filename, "r"); // open the configuration file
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else // read the configuration
		{
		err=fscanf(fp, "%d", &dimension);
		if(err!=1)
			{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(dimension != STDIM)
			{
			fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
				dimension, STDIM);
			exit(EXIT_FAILURE);
			}
	
		for(i=0; i<STDIM; i++)
			{
			err=fscanf(fp, "%d", &tmp_i);
			if(err!=1)
				{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			if(tmp_i != param->d_size[i])
				{
				fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
				exit(EXIT_FAILURE);
				}
			}
		
		err=fscanf(fp, "%*d %*d %d %d ", x_mu, x_nu);
		if(err!=2)
			{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		fclose(fp);
		}
	}

void free_gauge_conf(Gauge_Conf *GC, GParam const * const param)
	{
	long i;
	
	for(i=0; i<(param->d_volume); i++)
		{
		free(GC->lattice[i]);
		free(GC->Z[i]);
		}
	free(GC->lattice);
	free(GC->Z);

	#ifdef THETA_MODE
	end_clover_array(GC, param);
	#endif
	}
	
void free_replica(Gauge_Conf *GC, GParam const * const param)
	{
	int i;
	for(i=0; i<param->d_N_replica_pt; i++)
		{
		free_gauge_conf(&(GC[i]), param);
		free_bound_cond(&(GC[i]), param);
		}
	free(GC);
	}
	
void free_bound_cond(Gauge_Conf *GC, GParam const * const param)
	{
	long r;
	for(r=0; r<(param->d_volume); r++)
		{
		free(GC->C[r]);
		}
	free(GC->C);
	}

void free_twist_cond(Gauge_Conf *GC, GParam const * const param)
	{
	long r;
	for(r=0; r<(param->d_volume); r++)
		{
		free(GC->Z[r]);
		}
	free(GC->Z);
	}
	
// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
									GParam const * const param,
									char const * const namefile)
	{
	long si, lex;
	int i, mu;
	#ifdef HASH_MODE
	char md5sum[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum[2*STD_STRING_LENGTH+1]={0};
	#endif
	FILE *fp;
	
	#ifdef HASH_MODE
	compute_md5sum_conf(md5sum, GC, param);
	#endif
	
	fp=fopen(namefile, "w"); // open the configuration file
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		fprintf(fp, "%d ", STDIM);
		for(i=0; i<STDIM; i++)
			{
			fprintf(fp, "%d ", param->d_size[i]);
			}
		fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
		}
	fclose(fp);
	
	fp=fopen(namefile, "ab"); // open the configuration file in binary mode
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		for(lex=0; lex<param->d_volume; lex++)
			{
			si=lex_to_si(lex, param);
			for(mu=0; mu<STDIM; mu++)
				{
				print_on_binary_file_bigen(fp, &(GC->lattice[si][mu]) );
				}
			}
		fclose(fp);
		}
	}

void write_twist_on_file_with_name(Gauge_Conf const * const GC,
									GParam const * const param,
									char const * const namefile)
	{
	int i, j, mu, nu, cartcoord[STDIM], twisted_bc;
	FILE *fp;

	fp=fopen(namefile, "w"); // open the twist configuration file
	if(fp==NULL)
		{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	else
		{
		fprintf(fp, "%d ", STDIM);
		for(i=0; i<STDIM; i++)
			{
			fprintf(fp, "%d ", param->d_size[i]);
			}
		fprintf(fp, "\n");
		
		twisted_bc = 0;	// check if twist non trivial, save its plane (mu,nu)
		for(i=0; i<STDIM; i++)
			{
			cartcoord[i] = 0; // initialize to zero for later loop
			for(j=i+1; j<STDIM; j++)
				if (param->d_k_twist[dirs_to_si(i,j)] != 0)
						{
						twisted_bc = 1;
						mu = i;
						nu = j;
						}
			}
		
		if (twisted_bc == 1) // find twist cartcoord on plane cartcoord[i] = 0 for i != mu, nu (no need to read all volume)
			{
			for(i=0; i<param->d_size[mu]; i++)
				{
				cartcoord[mu] = i;
				for(j=0; j<param->d_size[nu]; j++)
					{
					cartcoord[nu] = j;
					if (cabs(GC->Z[cart_to_si(cartcoord, param)][dirs_to_si(mu,nu)] - (1.0+0.0*I))> MIN_VALUE)
						{
						fprintf(fp, "%d %d %d %d \n", mu, nu, i, j);
						}
					}
				}
			}
		fclose(fp);
		}
	}

void write_conf_on_file(Gauge_Conf const * const GC, GParam const * const param)
	{
	write_conf_on_file_with_name(GC, param, param->d_conf_file);
	write_twist_on_file_with_name(GC, param, param->d_twist_file);
	}

void write_replica_on_file(Gauge_Conf const * const GC, GParam const * const param)
	{
	int i;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif	
	for(i=0; i<param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH];
		strcpy(filename,param->d_conf_file);
		strcat(filename,"_replica_");
		sprintf(replica_index, "%d", i);
		strcat(filename,replica_index); // filename = d_conf_file + "_replica_${i}"
		write_conf_on_file_with_name(&(GC[i]),param,filename);
		
		strcpy(filename,param->d_twist_file);
		strcat(filename,"_replica_");
		strcat(filename,replica_index); // filename = d_twist_file + "_replica_${i}"
		write_twist_on_file_with_name(&(GC[i]),param,filename);
		}
	}
	
void write_replica_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
	{
	int i;
	static int counter=0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif
	for(i=0; i<param->d_N_replica_pt; i++)
		{
		char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH], aux_back[STD_STRING_LENGTH];
		if(counter==0)
			sprintf(aux_back, "_back0");
		else
			sprintf(aux_back, "_back1");
		strcpy(filename,param->d_conf_file);
		strcat(filename,"_replica_");
		sprintf(replica_index, "%d", i);
		strcat(filename, replica_index);
		strcat(filename, aux_back); // filename = d_conf_file + "_replica_${i}_back${counter}"
		write_conf_on_file_with_name(&(GC[i]),param,filename);
		
		strcpy(filename,param->d_twist_file);
		strcat(filename,"_replica_");
		strcat(filename, replica_index);
		strcat(filename, aux_back); // filename = d_conf_file + "_replica_${i}_back${counter}"
		write_twist_on_file_with_name(&(GC[i]),param,filename);
		}
	counter=1-counter;
	}

void write_conf_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
	{
	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	static int counter=0;

	if(counter==0)
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
	
	counter=1-counter;
	}


// allocate GC and initialize with GC2, including the twist factors
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param) 
	{
	long r;
	int j;

	// allocate the lattice and Z[r][j]
	allocate_array_GAUGE_GROUP_pointer(&(GC->lattice), param->d_volume, __FILE__, __LINE__);
	allocate_array_double_complex_pointer(&(GC->Z), param->d_volume, __FILE__, __LINE__);
	for(r=0; r<(param->d_volume); r++)
		{
		allocate_array_GAUGE_GROUP(&(GC->lattice[r]), STDIM, __FILE__, __LINE__);
		allocate_array_double_complex(&(GC->Z[r]), param->d_n_planes, __FILE__, __LINE__);
		}
	
	#ifdef THETA_MODE
	alloc_clover_array(GC, param);
	#endif
	
	// initialize GC and Z[r][j]
	for(r=0; r<(param->d_volume); r++)
		{
		for(j=0; j<STDIM; j++)
			{
			equal(&(GC->lattice[r][j]), &(GC2->lattice[r][j]) );
			}
		
		for(j=0; j<param->d_n_planes; j++)
			{
			GC->Z[r][j]=GC2->Z[r][j];
			}
		}

	GC->update_index=GC2->update_index;
	}


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Gauge_Conf const * const GC, GParam const * const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];
	long si, lex;
	GAUGE_GROUP matrix;
	int mu, k;

	MD5_Init(&mdContext);
	for(lex=0; lex<param->d_volume; lex++)
		{
		si=lex_to_si(lex, param);
		for(mu=0; mu<STDIM; mu++)
			{
			equal(&matrix, &(GC->lattice[si][mu]));

			#if NCOLOR==1
				MD5_Update(&mdContext, &(matrix.comp), sizeof(double complex));
			#elif NCOLOR==2
				for(k=0; k<4; k++)
					{
					MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
					}
			#else
				for(k=0; k<NCOLOR*NCOLOR; k++)
					{
					MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double complex));
					}
			#endif
			}
		}
	MD5_Final(c, &mdContext);

	for(k = 0; k < MD5_DIGEST_LENGTH; k++)
		{
		sprintf(&(res[2*k]), "%02x", c[k]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


// allocate the ml_polycorr arrays and related stuff
void alloc_polycorr_stuff(Gauge_Conf *GC,
								GParam const * const param)
	{
	int i, j;
	
	allocate_array_TensProd_pointer_pointer(&(GC->ml_polycorr), NLEVELS, __FILE__, __LINE__);
	for(i=0; i<NLEVELS; i++) 
		{
		allocate_array_TensProd_pointer(&(GC->ml_polycorr[i]), param->d_size[0]/param->d_ml_step[i], __FILE__, __LINE__);
		for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
			{
			allocate_array_TensProd(&(GC->ml_polycorr[i][j]), param->d_space_vol, __FILE__, __LINE__);
			}
		}
	
	allocate_array_GAUGE_GROUP_pointer(&(GC->loc_poly), param->d_size[0]/param->d_ml_step[NLEVELS-1], __FILE__, __LINE__);
	for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
		allocate_array_GAUGE_GROUP(&(GC->loc_poly[i]), param->d_space_vol, __FILE__, __LINE__);
	}


// free the ml_polycorr arrays and related stuff
void free_polycorr_stuff(Gauge_Conf *GC,
								GParam const * const param)
	{
	int i, j;

	for(i=0; i<NLEVELS; i++)
	{
	for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
		{
		free(GC->ml_polycorr[i][j]);
		}
	free(GC->ml_polycorr[i]);
	}
	free(GC->ml_polycorr);

	for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
	{
	free(GC->loc_poly[i]);
	}
	free(GC->loc_poly);
	}


// save ml_polycorr[0] arrays on file
void write_polycorr_on_file(Gauge_Conf const * const GC,
									GParam const * const param,
									int iteration)
	{
	long i;
	int j;
	#ifdef HASH_MODE
	char md5sum[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum[2*STD_STRING_LENGTH+1]={0};
	#endif
	FILE *fp;

	#ifdef HASH_MODE
	compute_md5sum_polycorr(md5sum, GC, param);
	#endif

	fp=fopen(param->d_ml_file, "w"); // open the configuration file
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
	}
	fclose(fp);

	fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}

	fclose(fp);
	}
	}


// read ml_polycorr[0] arrays from file
void read_polycorr_from_file(Gauge_Conf const * const GC,
									GParam const * const param,
									int *iteration)
	{
	long i, loc_space_vol;
	int j;
	FILE *fp;
	#ifdef HASH_MODE
	char md5sum_new[2*MD5_DIGEST_LENGTH+1];
	char md5sum_old[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum_old[2*STD_STRING_LENGTH+1]={0};
	#endif

	fp=fopen(param->d_ml_file, "r"); // open the multilevel file
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
	if(i!=3)
		{
		fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(loc_space_vol != param->d_space_vol)
		{
		fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
				param->d_ml_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}
	fclose(fp);

	fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	// read again the header: loc_space_vol, iteration, hash
	i=0;
	while(i!='\n')
			{
			i=fgetc(fp);
			}

	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}

	fclose(fp);
	}

	#ifdef HASH_MODE
	// compute the new md5sum and check for consistency
	compute_md5sum_polycorr(md5sum_new, GC, param);
	if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
		{
		fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
		exit(EXIT_FAILURE);
		}
	#endif
	}


// compute the md5sum of the ml_polycorr[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorr(char *res, Gauge_Conf const * const GC, GParam const * const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];
	long i;
	int j;
	int n1, n2, n3, n4;

	MD5_Init(&mdContext);

	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			for(n1=0; n1<NCOLOR; n1++)
				{
				for(n2=0; n2<NCOLOR; n2++)
					{
					for(n3=0; n3<NCOLOR; n3++)
						{
						for(n4=0; n4<NCOLOR; n4++)
							{
							MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
							}
						}
					}
				}
			}
		}

	MD5_Final(c, &mdContext);

	for(i = 0; i < MD5_DIGEST_LENGTH; i++)
		{
		sprintf(&(res[2*i]), "%02x", c[i]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


// allocate the ml_polycorradj arrays
void alloc_polycorradj(Gauge_Conf *GC,
							GParam const * const param)
	{
	int i, j;
	
	allocate_array_TensProdAdj_pointer_pointer(&(GC->ml_polycorradj), NLEVELS, __FILE__, __LINE__);
	for(i=0; i<NLEVELS; i++)
		{
		allocate_array_TensProdAdj_pointer(&(GC->ml_polycorradj[i]), param->d_size[0]/param->d_ml_step[i], __FILE__, __LINE__);
		for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
			{
			allocate_array_TensProdAdj(&(GC->ml_polycorradj[i][j]), param->d_space_vol, __FILE__, __LINE__);
			}
		}
	}


// free the ml_polycorradj arrays
void free_polycorradj(Gauge_Conf *GC,
							GParam const * const param)
	{
	int i, j;

	for(i=0; i<NLEVELS; i++)
	{
	for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
		{
		free(GC->ml_polycorradj[i][j]);
		}
	free(GC->ml_polycorradj[i]);
	}
	free(GC->ml_polycorradj);
	}


// allocate the ml_polycorr, polyplaq arrays and related stuff
void alloc_tube_disc_stuff(Gauge_Conf *GC,
									GParam const * const param)
	{
	int i;
	
	alloc_polycorr_stuff(GC, param);
	
	allocate_array_TensProd_pointer(&(GC->ml_polyplaq), NLEVELS, __FILE__, __LINE__);
	for(i=0; i<NLEVELS; i++) allocate_array_TensProd(&(GC->ml_polyplaq[i]), param->d_space_vol, __FILE__, __LINE__);
	
	allocate_array_double_complex(&(GC->loc_plaq), param->d_space_vol, __FILE__, __LINE__);
	}


// free the ml_polycorr, ml_polyplaq arrays and relates stuff
void free_tube_disc_stuff(Gauge_Conf *GC,
								GParam const * const param)
	{
	int i;
	
	free_polycorr_stuff(GC, param);
	for(i=0; i<NLEVELS; i++) free(GC->ml_polyplaq[i]);
	free(GC->ml_polyplaq);
	free(GC->loc_plaq);
	}


// allocate the polycorr, polyplaq, polyplaqconn arrays and related stuff
void alloc_tube_conn_stuff(Gauge_Conf *GC,
									GParam const * const param)
	{
	int i;
	
	alloc_tube_disc_stuff(GC, param);
	allocate_array_TensProd_pointer(&(GC->ml_polyplaqconn), NLEVELS, __FILE__, __LINE__);
	for(i=0; i<NLEVELS; i++) allocate_array_TensProd(&(GC->ml_polyplaqconn[i]), param->d_space_vol, __FILE__, __LINE__);
	allocate_array_GAUGE_GROUP(&(GC->loc_plaqconn), param->d_space_vol, __FILE__, __LINE__);
	}


// free the polycorr, polyplaq, polyplaqconn arrays and related stuff
void free_tube_conn_stuff(Gauge_Conf *GC,
								GParam const * const param)
	{
	int i, j;

	for(i=0; i<NLEVELS; i++)
	{
	for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
		{
		free(GC->ml_polycorr[i][j]);
		}
	free(GC->ml_polycorr[i]);
	}
	free(GC->ml_polycorr);

	free(GC->loc_plaqconn);
	}


// save ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays on file
void write_tube_conn_stuff_on_file(Gauge_Conf const * const GC,
											GParam const * const param,
											int iteration)
	{
	long i;
	int j;
	#ifdef HASH_MODE
	char md5sum[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum[2*STD_STRING_LENGTH+1]={0};
	#endif
	FILE *fp;

	#ifdef HASH_MODE
	compute_md5sum_tube_conn_stuff(md5sum, GC, param);
	#endif

	fp=fopen(param->d_ml_file, "w"); // open the configuration file
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
	}
	fclose(fp);

	fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}
	for(i=0; i<(param->d_space_vol); i++)
		{
		print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
		}
	for(i=0; i<(param->d_space_vol); i++)
		{
		print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
		}

	fclose(fp);
	}
	}


// read ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays from file
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
												GParam const * const param,
												int *iteration)
	{
	long i, loc_space_vol;
	int j;
	FILE *fp;
	#ifdef HASH_MODE
	char md5sum_new[2*MD5_DIGEST_LENGTH+1];
	char md5sum_old[2*MD5_DIGEST_LENGTH+1];
	#else
	char md5sum_old[2*STD_STRING_LENGTH+1]={0};
	#endif

	fp=fopen(param->d_ml_file, "r"); // open the multilevel file
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
	if(i!=3)
		{
		fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(loc_space_vol != param->d_space_vol)
		{
		fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
				param->d_ml_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}
	fclose(fp);

	fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
	if(fp==NULL)
	{
	fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	else
	{
	// read again the header: loc_space_vol, iteration, hash
	i=0;
	while(i!='\n')
			{
			i=fgetc(fp);
			}

	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
			}
		}
	for(i=0; i<(param->d_space_vol); i++)
		{
		read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
		}
	for(i=0; i<(param->d_space_vol); i++)
		{
		read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
		}

	fclose(fp);
	}

	#ifdef HASH_MODE
	// compute the new md5sum and check for consistency
	compute_md5sum_tube_conn_stuff(md5sum_new, GC, param);
	if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
		{
		fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
		exit(EXIT_FAILURE);
		}
	#endif
	}


// compute the md5sum of the ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param)
	{
	#ifdef HASH_MODE
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];
	long i;
	int j;
	int n1, n2, n3, n4;

	MD5_Init(&mdContext);

	for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
		{
		for(i=0; i<(param->d_space_vol); i++)
			{
			for(n1=0; n1<NCOLOR; n1++)
				{
				for(n2=0; n2<NCOLOR; n2++)
					{
					for(n3=0; n3<NCOLOR; n3++)
						{
						for(n4=0; n4<NCOLOR; n4++)
							{
							MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
							}
						}
					}
				}
			}
		}

	for(i=0; i<(param->d_space_vol); i++)
		{
		for(n1=0; n1<NCOLOR; n1++)
			{
			for(n2=0; n2<NCOLOR; n2++)
				{
				for(n3=0; n3<NCOLOR; n3++)
					{
					for(n4=0; n4<NCOLOR; n4++)
						{
						MD5_Update(&mdContext, &((GC->ml_polyplaq[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
						}
					}
				}
			}
		}

	for(i=0; i<(param->d_space_vol); i++)
		{
		for(n1=0; n1<NCOLOR; n1++)
			{
			for(n2=0; n2<NCOLOR; n2++)
				{
				for(n3=0; n3<NCOLOR; n3++)
					{
					for(n4=0; n4<NCOLOR; n4++)
						{
						MD5_Update(&mdContext, &((GC->ml_polyplaqconn[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
						}
					}
				}
			}
		}

	MD5_Final(c, &mdContext);

	for(i = 0; i < MD5_DIGEST_LENGTH; i++)
		{
		sprintf(&(res[2*i]), "%02x", c[i]);
		}
	#else
	// just to avoid warning at compile time
	(void) res;
	(void) GC;
	(void) param;
	#endif
	}


// allocate the clovers arrays
void alloc_clover_array(Gauge_Conf *GC,
							GParam const * const param)
	{
	int i;
	long r;
	
	allocate_array_GAUGE_GROUP_pointer_pointer(&(GC->clover_array), param->d_volume, __FILE__, __LINE__);
	for(r=0; r<param->d_volume; r++) 
		{
		allocate_array_GAUGE_GROUP_pointer(&(GC->clover_array[r]), STDIM, __FILE__, __LINE__);
		for(i=0; i<STDIM; i++) allocate_array_GAUGE_GROUP(&(GC->clover_array[r][i]), STDIM, __FILE__, __LINE__);
		}
	}


// free the clovers arrays
void end_clover_array(Gauge_Conf *GC,
							GParam const * const param)
	{
	int i;
	long r;

	for(r=0; r<param->d_volume; r++)
	{
	for(i=0; i<STDIM; i++)
		{
		free(GC->clover_array[r][i]);
		}
	free(GC->clover_array[r]);
	}
	free(GC->clover_array);
	}


#endif
