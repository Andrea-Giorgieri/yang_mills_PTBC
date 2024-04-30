#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"macro.h"

#include<complex.h>
#ifdef HASH_MODE
#include<openssl/md5.h>
#endif
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"
#include"u1.h"
#include"su2.h"
#include"sun.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

typedef struct Gauge_Conf {

	long update_index;

	GAUGE_GROUP **lattice;			// [volume] [STDIM]
	GAUGE_GROUP ***clover_array;	// [volume] [STDIM] [STDIM]
	
	// for the parallel tempering of volume defect
	double **C;		// C [volume] [STDIM], this factor changes the boundary condition for the link on the defect
	int conf_label;	// save the label of the configuration to keep track of the swaps

	// for the twisted boundary conditions
	double complex **Z;	// Z [volume] [STDIM*(STDIM-1)], this factor implements the twist for the links on the boundary
	
	// auxiliary conf for Metropolis tests and translations
	GAUGE_GROUP **lattice_copy;		// copy of lattice 
	double complex **Z_copy;		// copy of twist

	// for computing the polyakov loop correlator with multilevel
	TensProd ***ml_polycorr;	// [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]
	GAUGE_GROUP **loc_poly;		// [d_size[0]/d_ml_step[NLEVELS-1]] [space_vol] auxilliary vector to be used in the multilevel

	// for computing the polyakov loop correlator in the adjoint rep. with multilevel
	TensProdAdj ***ml_polycorradj;	// [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]

	// for the disconnected correlator for string width
	TensProd **ml_polyplaq;		// [NLEVELS] [only slice 0] [space_vol]
	double complex *loc_plaq;	// [only slice 0] [space_vol] auxilliary vector to be used in the multilevel

	// for the connected correlator for string width
	TensProd **ml_polyplaqconn;	// [NLEVELS] [only slice 0] [space_vol]
	GAUGE_GROUP *loc_plaqconn;	// [only slice 0][space_vol] auxilliary vector to be used in the multilevel
	
	// for multicanonical update: store running charge, copy of lattice to be cooled, rectangles for topcharge
	double stored_topcharge;
	GAUGE_GROUP **lattice_cold;			// aux lattice to be cooled
	GAUGE_GROUP **lattice_copy_cold;	// aux lattice to be cooled
} Gauge_Conf;
	
// to compute swap and multicanonic acceptances during parallel tempering evolution
// TO DO: refactor swaptrackfilep with this?
typedef struct Acc_Utils {
	long *num_accepted_swap;	// number of accepted swaps during parallel tempering
	long *num_swap;				// number of proposed swaps during parallel tempering
	double *metro_swap_prob;	// to store swap propabilities
	//FILE *swaptrackfilep;		// pointer to file to track replicas during parallel tempering
	
	long *num_accepted_metro_multicanonic;	// number of accepted multicanonic Metropolis updates
	long *num_metro_multicanonic;			// number of proposed multicanonic Metropolis updates
	FILE *multicanonic_acc_filep;			// pointer to file to write acceptance of multicanonic Metropolis updates
	}	Acc_Utils;

// auxiliary arrays and files for measures
typedef struct Meas_Utils {
	// for measurement of localobs
	double *meanplaq;
	double *clover_energy;
	double *charge;
	double *sum_q_timeslices;
	double *chi_prime;
	double **charge_prime;

	// for adaptive gradienf flow
	GAUGE_GROUP **lattice_aux[4];		// array of 4 lattices
	double local_max_dist[NTHREADS];
	
	// pointers to data files
	FILE *datafilep;
	FILE *chiprimefilep;
	FILE *topchar_tcorr_filep;
	}	Meas_Utils;


// in gauge_conf_def.c
void equal_gauge_conf(						Gauge_Conf *GC1,
											Gauge_Conf *GC2,
											GParam const * const param);

void accept_gauge_conf(						Gauge_Conf * const GC,
											GParam const * const param);

void restore_gauge_conf(					Gauge_Conf * const GC,
											GParam const * const param);

void accept_gauge_conf_rectangle(			Gauge_Conf * const GC,
											int const hierarc_level,
											Rect_Utils const * const rect_aux);

void restore_gauge_conf_rectangle(			Gauge_Conf * const GC,
											int const hierarc_level,
											Rect_Utils const * const rect_aux);

void init_gauge_conf_from_file_with_name(	Gauge_Conf *GC,
											GParam const * const param,
											char const * const filename);
											
void init_gauge_conf(						Gauge_Conf *GC,
											Geometry const * const geo,
											GParam const * const param);
											
void init_gauge_conf_step(					Gauge_Conf *GC,
											GParam const * const param,
											long step,
											int *stop);

void read_gauge_conf_step(					Gauge_Conf *GC,
											GParam const * const param,
											long step,
											int *stop);
											
void init_gauge_conf_replica(				Gauge_Conf **GC,
											Geometry const * const geo,
											GParam const * const param);
											
void init_bound_cond(						Gauge_Conf *GC,
											GParam const * const param,
											int const a);
											
void init_twist_cond_from_file_with_name(	Gauge_Conf *GC, GParam const * const param,
											char const * const filename);
											
void free_replica(							Gauge_Conf *GC,
											GParam const * const param);
											
void free_bound_cond(						Gauge_Conf *GC,
											GParam const * const param);
											
void free_twist_cond(						Gauge_Conf *GC,
											GParam const * const param);
											
void read_gauge_conf_from_file_with_name(	Gauge_Conf *GC,
											GParam const * const param, char const * const filename);
											
void read_twist_cond_from_file_with_name(	int *x_mu, int *x_nu, 
											GParam const * const param, char const * const filename);
											
void free_gauge_conf(						Gauge_Conf *GC,
											GParam const * const param);
											
void write_conf_on_file_with_name(			Gauge_Conf const * const GC,
											GParam const * const param,
											char const * const namefile);
											
void write_twist_on_file_with_name(		Gauge_Conf const * const GC,
										GParam const * const param,
										char const * const namefile);
										
void write_conf_on_file(				Gauge_Conf const * const GC,
										GParam const * const param);
										
void write_conf_on_file_back(			Gauge_Conf const * const GC,
										GParam const * const param);
										
void write_replica_on_file(				Gauge_Conf const * const GC,
										GParam const * const param);
										
void write_replica_on_file_back(		Gauge_Conf const * const GC,
										GParam const * const param);
										
void init_gauge_conf_from_gauge_conf(	Gauge_Conf *GC,
										Gauge_Conf const * const GC2,
										GParam const * const param);

void init_twist_cond_from_twist_cond(	double complex **Z, 
										double complex const * const * const Z2,
										GParam const * const param);
										
void compute_md5sum_conf(		char *res,		// the lenght is 2*MD5_DIGEST_LENGTH
								Gauge_Conf const * const GC,
								GParam const * const param);
								
void alloc_polycorr_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void free_polycorr_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void write_polycorr_on_file(	Gauge_Conf const * const GC,
								GParam const * const param,
								int iteration);
								
void read_polycorr_from_file(	Gauge_Conf const * const GC,
								GParam const * const param,
								int *iteration);
								
void compute_md5sum_polycorr(	char *res,		// the lenght is 2*MD5_DIGEST_LENGTH
								Gauge_Conf const * const GC,
								GParam const * const param);
								
void alloc_polycorradj(			Gauge_Conf *GC,
								GParam const * const param);
								
void free_polycorradj(			Gauge_Conf *GC,
								GParam const * const param);
								
void alloc_tube_disc_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void free_tube_disc_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void alloc_tube_conn_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void free_tube_conn_stuff(		Gauge_Conf *GC,
								GParam const * const param);
								
void write_tube_conn_stuff_on_file(	Gauge_Conf const * const GC,
									GParam const * const param,
									int iteration);
									
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
									GParam const * const param,
									int *iteration);
									
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param);

void alloc_clover_array(			Gauge_Conf *GC,
									GParam const * const param);
									
void end_clover_array(				Gauge_Conf *GC,
									GParam const * const param);
									

// in gauge_conf_meas.c
double plaquettep(					Gauge_Conf const * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i,
									int j);
					
double complex plaquettep_complex(	Gauge_Conf const * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i,
									int j);
									
void plaquettep_matrix(				Gauge_Conf const * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i,
									int j,
									GAUGE_GROUP *matrix);
						
void clover(			Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						long r,
						int i,
						int j,
						GAUGE_GROUP *M);
			
void plaquette(			Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						double *plaqs,
						double *plaqt);
				
void clover_disc_energy(Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						double *energy);
						
void polyakov(			Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						double *repoly,
						double *impoly);
				
void polyakov_adj(		Gauge_Conf const * const GC,
						Geometry const * const geo,
						GParam const * const param,
						double *repoly,
						double *impoly);
					
void polyakov_with_tracedef(			Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										double *repoly,
										double *impoly);

void check_correlation_decay_cooling(	Gauge_Conf const * const GC, 
										Geometry const * const geo, 
										GParam const * const param,
										double *ratio);

double sum_abs_topcharge_dens(			Gauge_Conf const * const GC,
										Geometry const * const geo, 
										GParam const * const param);

double loc_topcharge(		Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							long r);

double delta_loc_topcharge(		Gauge_Conf const * const GC1,
								Gauge_Conf const * const GC2,
								Geometry const * const geo,
								GParam const * const param,
								long r);
					 
double topcharge(			Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param);

double topcharge_rectangle(	Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							Rectangle const * const topcharge_rect);

double topcharge_prime(		Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param, int const dir);
				 
void topcharge_timeslices(	Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param, double *ris, int ncool, FILE*);
				 
void topcharge_timeslices_cooling(	Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux);
				 
void topcharge_timeslices_gradflow(	Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux);
				 
double topo_chi_prime(				Gauge_Conf const * const GC,
									Geometry const * const geo,
									GParam const * const param);
				 
void topo_obs_cooling(				Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux);
						
void topo_obs_gradflow(				Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Meas_Utils *meas_aux);
						
void topo_obs_clover_energy_gradflow(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
									
void topcharge_cooling(					Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
						
void topcharge_gradflow(				Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
						
void topcharge_clover_energy_gradflow(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
									
void loc_topcharge_corr(				Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int ncool,
										int dist,
										double *ris);
					
void perform_measures_aux(				Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int const meas_count,
										Meas_Utils *meas_aux);
							
void perform_measures_localobs(			Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
								
void perform_measures_localobs_notopo(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
								
void perform_measures_localobs_cooling(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
							   
void perform_measures_localobs_with_gradflow(					Gauge_Conf * const GC,
																Geometry const * const geo,
																GParam const * const param,
																Meas_Utils *meas_aux);
								
void perform_measures_localobs_with_adaptive_gradflow(			Gauge_Conf * const GC,
																Geometry const * const geo,
																GParam const * const param,
																Meas_Utils *meas_aux);
											
void perform_measures_localobs_with_adaptive_gradflow_debug(	Gauge_Conf * const GC,
																Geometry const * const geo,
																GParam const * const param,
																Meas_Utils *meas_aux,
																FILE *step_filep);
											
void perform_measures_localobs_with_adaptive_gradflow_debug2(	Gauge_Conf * const GC,
																Geometry const * const geo,
																GParam const * const param,
																Meas_Utils *meas_aux,
																FILE *step_filep);
											
void perform_measures_localobs_with_tracedef(	Gauge_Conf * const GC,
												Geometry const * const geo,
												GParam const * const param,
												Meas_Utils *meas_aux);
											 
void optimize_multihit_polycorr(				Gauge_Conf * const GC,
												Geometry const * const geo,
												GParam const * const param,
												FILE *datafilep);
								
void optimize_multilevel_polycorr(				Gauge_Conf * const GC,
												Geometry const * const geo,
												GParam const * const param,
												FILE *datafilep);
									
void perform_measures_polycorr(			Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
								
void optimize_multihit_polycorradj(		Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										FILE *datafilep);
									
void optimize_multilevel_polycorradj(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										FILE *datafilep);
									 
void perform_measures_polycorradj(		Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
									
void perform_measures_polycorr_long(	Gauge_Conf * const GC,
										GParam const * const param,
										Meas_Utils *meas_aux);
									
void optimize_multilevel_polycorr_long(	Gauge_Conf * const GC,
										GParam const * const param,
										FILE *datafilep);
										
void perform_measures_polycorr_long(	Gauge_Conf * const GC,
										GParam const * const param,
										Meas_Utils *meas_aux);
									
void perform_measures_tube_disc(		Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
								
void perform_measures_tube_conn(		Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Meas_Utils *meas_aux);
								
void perform_measures_tube_conn_long(	Gauge_Conf * const GC,
										GParam const * const param,
										Meas_Utils *meas_aux);
									 
void allocate_measures_arrays(			int const num_meas, GParam const * const param, double **meanplaq,
										double **clover_energy, double **charge, double **sum_q_timeslices,
										double **chi_prime, double ***charge_prime);
								
void allocate_measures_arrays_cooling(	int const num_meas, GParam const * const param,
										double **meanplaq, double **charge, double **chi_prime);

void open_data_file(			Meas_Utils *meas_aux,
								char const * const data,
								char const * const chiprime,
								char const * const topchar_tcorr,
								GParam const * const param);

void init_meas_utils(			Meas_Utils *meas_aux,
								GParam const * const param,
								int const replica_index);

void init_meas_utils_replica(	Meas_Utils **meas_aux,
								GParam const * const param);

void free_meas_utils(			Meas_Utils meas_aux,
								GParam const * const param,
								int const replica_index);

void free_meas_utils_replica(	Meas_Utils *meas_aux,
								GParam const * const param);

void print_measures_arrays(		int const num_meas, long const update_index, GParam const * const param,
								Meas_Utils *meas_aux);

// TO DO: remove after refactoring						
void free_measures_arrays(		int const num_meas, GParam const * const param, double *meanplaq, double *clover_energy,
								double *charge, double *sum_q_timeslices,
								double *chi_prime, double **charge_prime);


// in gauge_conf_multilevel.c
void multihit(								Gauge_Conf const * const GC,
											Geometry const * const geo,
											GParam const * const param,
											long r,
											int dir,
											int num_hit,
											GAUGE_GROUP *G);
				
void compute_local_poly(					Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param);
						
void update_for_multilevel(					Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int level);
							
void multilevel_polycorr(					Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt);
							
void multilevel_polycorradj(				Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt);
							
void multilevel_polycorr_long(				Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt,
											int iteration);
								
void compute_local_poly_and_plaq(			Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param);
								 
void multilevel_tube_disc(					Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt);
							
void compute_local_poly_plaq_and_plaqconn(	Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param);
											
void multilevel_tube_conn(					Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt);
							
void multilevel_tube_conn_long(				Gauge_Conf * const GC,
											Geometry const * const geo,
											GParam const * const param,
											int dt,
											int iteration);

// in gauge_conf_upd.c
void calcstaples_wilson(				Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const gparam,
										long r,
										int i,
										GAUGE_GROUP *M);

void calcstaples_tracedef(				Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										long r,
										int i,
										GAUGE_GROUP * M);

void calcstaples_with_topo(				Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										long r,
										int i,
										GAUGE_GROUP *M);

void compute_clovers(					Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int direction);

void compute_clovers_replica(			Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int dir);

void compute_clovers_replica_rect(		Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int dir,
										Rectangle const * const clover_rectangle);

void calcstaples_wilson_with_defect(	Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										long r,
										int i,
										GAUGE_GROUP *M);

void calcstaples_with_topo_with_defect(	Gauge_Conf const * const GC,
										Geometry const * const geo,
										GParam const * const param,
										long r,
										int i,
										GAUGE_GROUP *M);


void heatbath(						Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);

void overrelaxation(				Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);

int metropolis(						Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);

int metropolis_with_tracedef(		Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);

void heatbath_with_defect(			Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);

void overrelaxation_with_defect(	Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									long r,
									int i);


void update(						Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Acc_Utils *acc_counters);

void update_with_defect(			Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									Acc_Utils *acc_counters);

void update_rectangle_with_defect(	Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									int const hierarc_level,
									Rect_Utils const * const rect_aux,
									Acc_Utils *acc_counters);

void update_with_trace_def(			Gauge_Conf * const GC,
									Geometry const * const geo,
									GParam const * const param,
									double *acc_dt);


void hierarchical_update_rectangle_with_defect(		Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
													int const hierarc_level,
													Rect_Utils const * const rect_aux,
													Acc_Utils *acc_counters);

void parallel_tempering_with_hierarchical_update(	Gauge_Conf * const GC, Geometry const * const geo, GParam const * const param,
													Rect_Utils const * const rect_aux,
													Acc_Utils *acc_counters);


void cooling(							Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										int const n);

void hierarchical_cooling(				Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const * const param,
										Rectangle const * const cooling_rect);

void gradflow_RKstep(					Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const *const param,
										double dt,
										Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive(			Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const *const param,
										double *t,
										double *dt,
										int *accepted,
										Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive_debug(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const *const param,
										double *t,
										double *dt,
										int *accepted,
										double *total_error,
										Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive_debug2(	Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const *const param,
										double *t,
										double *dt,
										int *accepted,
										double *total_error,
										Meas_Utils *meas_aux);

void ape_smearing(						Gauge_Conf * const GC,
										Geometry const * const geo,
										GParam const *const param,
										double alpha,
										int n);
									
// in gauge_conf_paral_temp.c
void swap(					Gauge_Conf * const GC,
							Geometry const * const geo,
							GParam const * const param,
							Rectangle const * const swap_rectangle,
							Acc_Utils *acc_counters);

double delta_action_swap(	Gauge_Conf const * const GC,
							Geometry const * const geo,
							GParam const * const param,
							long const r,
							int const i,
							int const j,
							int const a,
							int const b);

void metropolis_single_swap(Gauge_Conf * const GC,
							int const a,
							int const b,
							double const p,
							Acc_Utils *acc_counters);

void conf_translation(		Gauge_Conf * const GC,
							Geometry const * const geo,
							GParam const * const param);	

void init_acc_utils(		Acc_Utils *acc_counters,
							GParam const * const param);

void free_acc_utils(		Acc_Utils *acc_counters,
							GParam const * const param);

void print_acceptances(		Acc_Utils const * const acc_counters,
							GParam const * const param);

void init_swap_track_file(	FILE **swaptrackfilep,
							GParam const * const param);

void print_conf_labels(		FILE *fp,
							Gauge_Conf const * const GC,
							GParam const * const param);

// in gauge_conf_multicanonic.c
void init_multicanonic_gauge_conf(	Gauge_Conf * const GC, 
									Geometry const * const geo, 
									GParam const * const param);

void init_multicanonic_acc_utils(	Acc_Utils *,
									GParam const * const);

void free_multicanonic_acc_utils(	Acc_Utils *);

void print_multicanonic_acceptance(	Gauge_Conf const * const,
									GParam const * const,
									Acc_Utils const * const);


double compute_topo_potential(			int const,
										double,
										GParam const * const);

double multicanonic_topcharge_cooling(	Gauge_Conf * const,
										Geometry const * const,
										GParam const * const);

double multicanonic_topcharge_cooling_rectangle(	Gauge_Conf * const,
													Geometry const * const,
													GParam const * const,
													int const,
													Rect_Utils const * const);

double topcharge_agf_multicanonic(		Gauge_Conf * const,
										Geometry const * const,
										GParam const * const,
										Meas_Utils *);


double metropolis_prob_multicanonic(int const,
									double const,
									double const,
									GParam const * const);

double delta_topo_potential_swap(	Gauge_Conf const * const,
									int const,
									int const,
									GParam const * const);


int multicanonic_metropolis_step_single_link(	Gauge_Conf * const,
												Geometry const * const,
												GParam const * const,
												long, int,
												GAUGE_GROUP);

int multicanonic_metropolis_step_all_links(		Gauge_Conf * const,
												Geometry const * const,
												GParam const * const);

int multicanonic_metropolis_step_rectangle(		Gauge_Conf * const,
												Geometry const * const,
												GParam const * const,
												int const,
												Rect_Utils const * const);

// old implementation of multicanonic

void compute_single_clover_insertion(	GAUGE_GROUP *,
										Gauge_Conf const * const,
										Geometry const * const,
										GParam const * const,
										long, int, int);

void compute_topostaple_alone(			Gauge_Conf const * const,
										Geometry const * const,
										GParam const * const,
										long, int,
										GAUGE_GROUP *);	
										
double delta_Q_upd(						Gauge_Conf const * const,
										Geometry const * const,
										GParam const * const,
										long, int,
										GAUGE_GROUP);

void multicanonic_parallel_tempering_with_hierarchical_update(	Gauge_Conf *,
																Geometry const * const,
																GParam const * const,
																Rectangle const * const,
																Rectangle const * const,
																Rectangle const * const,
																Acc_Utils *);

void multicanonic_agf_parallel_tempering_with_hierarchical_update(	Gauge_Conf *,
																	Geometry const * const,
																	GParam const * const,
																	Rectangle const * const,
																	Rectangle const * const,
																	Rectangle const * const,
																	Acc_Utils *);

void multicanonic_update_with_defect(				Gauge_Conf *,
													Geometry const * const,
													GParam const * const,
													Acc_Utils *);

void multicanonic_agf_update_with_defect(			Gauge_Conf *,
													Geometry const * const,
													GParam const * const,
													Acc_Utils *);

void multicanonic_update_rectangle_with_defect(		Gauge_Conf *,
													Geometry const * const,
													GParam const * const,
													Rectangle const * const,
													Rectangle const * const,
													Acc_Utils *);

void multicanonic_agf_update_rectangle_with_defect(	Gauge_Conf *,
													Geometry const * const,
													GParam const * const,
													Rectangle const * const,
													Rectangle const * const,
													Acc_Utils *);

void multicanonic_hierarchical_update_rectangle_with_defect(	Gauge_Conf *,
																Geometry const * const,
																GParam const * const,
																int const,
																Rectangle const * const,
																Rectangle const * const,
																Rectangle const * const,
																Acc_Utils *);

void multicanonic_agf_hierarchical_update_rectangle_with_defect(Gauge_Conf *,
																Geometry const * const,
																GParam const * const,
																int const,
																Rectangle const * const,
																Rectangle const * const,
																Rectangle const * const,
																Acc_Utils *);

int multicanonic_overrelaxation_with_defect(Gauge_Conf *,
											Geometry const * const,
											GParam const * const,
											long, int);
											
int multicanonic_heatbath_with_defect(		Gauge_Conf *,
											Geometry const * const,
											GParam const * const,
											long, int);

#endif
