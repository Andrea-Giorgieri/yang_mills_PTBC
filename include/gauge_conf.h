#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include "macro.h"

#include <complex.h>
#include <stdio.h>
#ifdef HASH_MODE
#include <openssl/md5.h>
#endif

#include "gauge_group.h"
#include "geometry.h"
#include "gparam.h"
#include "tens_prod.h"


typedef struct Gauge_Conf
	{
	long update_index;
	int replica_index; // index of the parallel tempering replica (defect coefficient), exchanged by the swaps
	int device;        // device (index in g_dev) holding this configuration, fixed for the whole run

	GAUGE_GROUP **lattice;       // [volume] [STDIM]
	GAUGE_GROUP ***clover_array; // [volume] [STDIM] [STDIM]

	// for the boundary conditions
	double complex **Z; // Z [volume] [STDIM*(STDIM-1)+1], this factor implements twisted and open (including the bulk identification) boundary conditions

	// for tracking random translations
	long translation[STDIM];

	// auxiliary conf for Metropolis tests and translations
	GAUGE_GROUP **lattice_copy; // copy of lattice
	double complex **Z_copy;    // copy of boundary conditions

	// for computing the polyakov loop correlator with multilevel
	TensProd ***ml_polycorr; // [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]
	GAUGE_GROUP **loc_poly;  // [d_size[0]/d_ml_step[NLEVELS-1]] [space_vol] auxiliary vector to be used in the multilevel

	// for the disconnected correlator for string width
	TensProd **ml_polyplaq;   // [NLEVELS] [only slice 0] [space_vol]
	double complex *loc_plaq; // [only slice 0] [space_vol] auxiliary vector to be used in the multilevel

	// for the connected correlator for string width
	TensProd **ml_polyplaqconn; // [NLEVELS] [only slice 0] [space_vol]
	GAUGE_GROUP *loc_plaqconn;  // [only slice 0][space_vol] auxiliary vector to be used in the multilevel

	// for multicanonical update: store running charge, copy of lattice to be cooled, rectangles for topcharge
	double stored_topcharge;
	GAUGE_GROUP **lattice_cold;      // aux lattice to be cooled
	GAUGE_GROUP **lattice_copy_cold; // aux lattice to be cooled

	// for symmetrical cooling update: store a random permutation of spacetime directions and parity of sites updated first
	int stdim_shuffle[STDIM];
	long parity_shuffle[2][2];

	#ifdef RNG_PHILOX_MODE
	// random number generators of the link updates, one stream per lattice site [volume]
	// (NULL for auxiliary configurations that are never updated with heatbath/overrelaxation/metropolis)
	rng_t *rng_site;
	#endif
	} Gauge_Conf;


// generator to be used for the update of the links of site r of configuration GCp:
// Philox: the stream of the site; dSFMT: the generator of the calling thread (as in the original code)
#ifdef RNG_PHILOX_MODE
#define RNG_SITE(GCp, r) (&((GCp)->rng_site[(r)]))
#else
#define RNG_SITE(GCp, r) rng_thread()
#endif


// boundary condition coefficient of the replica: d_pt_bound_cond_coeff[replica_index]
ACC_ROUTINE_SEQ
static inline double replica_defect_coeff(Gauge_Conf const *const GC, GParam const *const param)
	{
	return param->d_pt_bound_cond_coeff[GC->replica_index];
	}


// boundary condition factor of the link (r, i) for defect coefficient c: c on the defect, 1 elsewhere
ACC_ROUTINE_SEQ
static inline double link_defect_coeff(Geometry const *const geo, GParam const *const param, long const r, int const i, double const c)
	{
	return (i == param->d_defect_dir && geo->d_defect_mask[r] == 1) ? c : 1.0;
	}


// boundary condition factor K_ij(s) of the plaquette (s, i, j) for defect coefficient c, with spi = s+i and spj = s+j
ACC_ROUTINE_SEQ
static inline double plaquette_defect_coeff(Geometry const *const geo, GParam const *const param,
                                            long const s, long const spi, long const spj, int const i, int const j, double const c)
	{
	return link_defect_coeff(geo, param, s, i, c) * link_defect_coeff(geo, param, spi, j, c) * link_defect_coeff(geo, param, spj, i, c) * link_defect_coeff(geo, param, s, j, c);
	}


// number n of links of the plaquette (s, i, j) on the defect, with spi = s+i and spj = s+j: K_ij(s) = c^n
// (same links as plaquette_defect_coeff; n <= 2 because only one of the directions i != j is the defect direction)
ACC_ROUTINE_SEQ
static inline int plaquette_defect_links(Geometry const *const geo, GParam const *const param,
                                         long const s, long const spi, long const spj, int const i, int const j)
	{
	int const mu = param->d_defect_dir;
	return (i == mu && geo->d_defect_mask[s] == 1) + (j == mu && geo->d_defect_mask[spi] == 1) + (i == mu && geo->d_defect_mask[spj] == 1) + (j == mu && geo->d_defect_mask[s] == 1);
	}


// to compute swap and multicanonic acceptances during parallel tempering evolution
typedef struct Acc_Utils
	{
	// replica_to_gc[a] = position in GC of the configuration with replica_index a, inverse of GC[i].replica_index:
	// replica_to_gc[GC[i].replica_index] = i, GC[replica_to_gc[a]].replica_index = a
	int *replica_to_gc;

	long *num_accepted_swap; // number of accepted swaps during parallel tempering
	long *num_swap;          // number of proposed swaps during parallel tempering
	double *acc_swap;        // acceptance of swaps during parallel tempering
	double *metro_swap_prob; // to store swap probabilities
	FILE *swaptrackfilep;    // pointer to file to track replicas during parallel tempering

	long *num_accepted_metro_multicanonic; // number of accepted multicanonic Metropolis updates
	long *num_metro_multicanonic;          // number of proposed multicanonic Metropolis updates
	FILE *multicanonic_acc_filep;          // pointer to file to write acceptance of multicanonic Metropolis updates
	} Acc_Utils;


// auxiliary arrays and files for measures
typedef struct Meas_Utils
	{
	// for measurement of localobs
	double *meanplaq;
	double *clover_energy;
	double *charge;
	double **polyre;
	double **polyim;
	double *multipolyre;
	double *multipolyim;
	double *polyre_density;
	double *polyim_density;
	double *chi_prime;
	double **charge_prime;
	double *real_slices;
	double *imag_slices;
	double *scalar_density;
	double *action1;
	double *action2;
	double *action3;
	double *potential;

	// for adaptive gradient flow
	GAUGE_GROUP **lattice_aux[4]; // array of 4 lattices

	// pointers to data files
	FILE *datafilep;
	FILE *energydensityfilep;
	FILE *chargedensityfilep;
	FILE *polyakovdensityfilep[STDIM];
	FILE *chiprimefilep;
	FILE *e_slices_filep;
	FILE *q_slices_filep;
	} Meas_Utils;

// auxiliary arrays and flags for tuning of topo_potential
typedef struct Tune_Utils
	{
	double *grid_coarse;
	double *grid_fine;
	double *tuning_stp;
	double *grid_coarse_aux;
	int d_n_grid_aux;
	int *sweep_counter;
	int *sweep_check;
	int *sweep_flag;
	} Tune_Utils;


// in gauge_conf_def.c
void allocate_lattice_with_copy(Gauge_Conf *GC,
                                GParam const *const param);

void allocate_lattice_cold_with_copy(Gauge_Conf *GC,
                                     GParam const *const param);

void allocate_Z_with_copy(Gauge_Conf *GC,
                          GParam const *const param);

void allocate_clover_array(Gauge_Conf *GC,
                           GParam const *const param);

void equal_lattice(GAUGE_GROUP *const *const lattice1,
                   GAUGE_GROUP const *const *const lattice2,
                   GParam const *const param);

void equal_equal_lattice(GAUGE_GROUP *const *const lattice1,
                         GAUGE_GROUP *const *const lattice2,
                         GAUGE_GROUP const *const *const lattice3,
                         GParam const *const param);

double lattice_total_dist(GAUGE_GROUP const *const *const lattice1,
                          GAUGE_GROUP const *const *const lattice2,
                          GParam const *const param);

double lattice_max_dist(GAUGE_GROUP const *const *const lattice1,
                        GAUGE_GROUP const *const *const lattice2,
                        GParam const *const param);

void read_gauge_conf_from_file_with_name(Gauge_Conf *GC,
                                         GParam const *const param,
                                         char const *const filename);

void init_gauge_conf_from_file_with_name(Gauge_Conf *GC,
                                         GParam const *const param,
                                         char const *const filename);

void read_twist_cond_from_file_with_name(int *x_mu, int *x_nu, int *x_obc, int translation[STDIM],
                                         GParam const *const param,
                                         char const *const filename);

void initialize_Z_with_copy(Gauge_Conf *GC, GParam const *const param,
                            int x_mu, int x_nu, int x_obc, int const translation[STDIM]);

void init_twist_cond_from_file_with_name(Gauge_Conf *GC,
                                         GParam const *const param,
                                         char const *const filename);

void init_gauge_conf(Gauge_Conf *GC,
                     Geometry const *const geo,
                     GParam const *const param);

void free_gauge_conf(Gauge_Conf *GC,
                     GParam const *const param);

void init_rng_gauge_conf(Gauge_Conf *GC,
                         GParam const *const param,
                         unsigned int stream);

void free_rng_gauge_conf(Gauge_Conf *GC);

int init_gauge_conf_step(Gauge_Conf *GC,
                         GParam const *const param,
                         long step);

int read_gauge_conf_step(Gauge_Conf *GC,
                         GParam const *const param,
                         long step);

void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC,
                                     Gauge_Conf const *const GC2,
                                     GParam const *const param);

void init_gauge_conf_replica(Gauge_Conf **GC,
                             Geometry const *const geo,
                             GParam const *const param);

void free_gauge_conf_replica(Gauge_Conf *GC,
                             GParam const *const param);

void equal_gauge_conf(Gauge_Conf *GC1,
                      Gauge_Conf *GC2,
                      GParam const *const param);

void accept_gauge_conf(Gauge_Conf *const GC,
                       GParam const *const param);

void restore_gauge_conf(Gauge_Conf *const GC,
                        GParam const *const param);

void accept_gauge_conf_rectangle(Gauge_Conf *const GC,
                                 int const hierarc_level,
                                 Rect_Utils const *const rect_aux);

void restore_gauge_conf_rectangle(Gauge_Conf *const GC,
                                  int const hierarc_level,
                                  Rect_Utils const *const rect_aux);

void accept_or_restore_gauge_conf_device(Gauge_Conf *const GC,
                                         GParam const *const param,
                                         int const acc);

void accept_or_restore_gauge_conf_rectangle_device(Gauge_Conf *const GC,
                                                   int const hierarc_level,
                                                   Rect_Utils const *const rect_aux,
                                                   int const acc);

void write_conf_on_file_with_name(Gauge_Conf const *const GC,
                                  GParam const *const param,
                                  char const *const filename);

void write_twist_on_file_with_name(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   char const *const filename);

void write_conf_on_file(Gauge_Conf const *const GC,
                        GParam const *const param);

void write_conf_on_file_back(Gauge_Conf const *const GC,
                             GParam const *const param);

void write_replica_on_file(Gauge_Conf const *const GC,
                           GParam const *const param);

void write_replica_on_file_back(Gauge_Conf const *const GC,
                                GParam const *const param);

void allocate_polycorr_stuff(Gauge_Conf *GC,
                             GParam const *const param);

void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const *const param);

void allocate_tube_disc_stuff(Gauge_Conf *GC,
                              GParam const *const param);

void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const *const param);

void allocate_tube_conn_stuff(Gauge_Conf *GC,
                              GParam const *const param);

void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const *const param);

void allocate_multilevel_arrays(Gauge_Conf *GC,
                                GParam const *const param,
                                Multilevel_Obs const ml_obs);

void free_multilevel_arrays(Gauge_Conf *GC,
                            GParam const *const param,
                            Multilevel_Obs const ml_obs);

void write_polycorr_on_file(Gauge_Conf const *const GC,
                            GParam const *const param,
                            int const iteration);

void read_polycorr_stuff_from_file(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   int *iteration);

void write_tube_conn_stuff_on_file(Gauge_Conf const *const GC,
                                   GParam const *const param,
                                   int const iteration);

void read_tube_conn_stuff_from_file(Gauge_Conf const *const GC,
                                    GParam const *const param,
                                    int *iteration);

void write_multilevel_status_on_file(Gauge_Conf const *const GC,
                                     GParam const *const param,
                                     int const iteration,
                                     Multilevel_Obs const ml_obs);

void read_multilevel_status_from_file(Gauge_Conf const *const GC,
                                      GParam const *const param,
                                      int *iteration,
                                      Multilevel_Obs const ml_obs);

void compute_md5sum_conf(char *res,
                         Gauge_Conf const *const GC,
                         GParam const *const param);

void compute_md5sum_polycorr(char *res,
                             Gauge_Conf const *const GC,
                             GParam const *const param);

void compute_md5sum_tube_conn_stuff(char *res,
                                    Gauge_Conf const *const GC,
                                    GParam const *const param);


// in gauge_conf_meas.c
ACC_ROUTINE_SEQ
double plaquettep(Gauge_Conf const *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  long r,
                  int i,
                  int j);

double complex plaquettep_complex(Gauge_Conf const *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  long r,
                                  int i,
                                  int j);

ACC_ROUTINE_SEQ
void plaquettep_matrix(Gauge_Conf const *const GC,
                       Geometry const *const geo,
                       GParam const *const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix);

void plaquette(Gauge_Conf const *const GC,
               Geometry const *const geo,
               GParam const *const param,
               double *plaqs,
               double *plaqt);

ACC_ROUTINE_SEQ
void clover(Gauge_Conf const *const GC,
            Geometry const *const geo,
            GParam const *const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M);

ACC_ROUTINE_SEQ
double loc_clover_energy(Gauge_Conf const *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         long const r);

void energy_density(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    Meas_Utils *const meas_aux,
                    int const meas_count);

void clover_disc_energy(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        double *const energy);

void clover_energy_slices(Gauge_Conf const *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          int const mu,
                          double *slices,
                          int meas_count,
                          FILE *filep);

void polyakov_density(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param,
                      int mu,
                      Meas_Utils *const meas_aux,
                      int const meas_count);

void polyakov(Gauge_Conf const *const GC,
              Geometry const *const geo,
              GParam const *const param,
              int mu,
              int pwr,
              double *repoly,
              double *impoly);

void multipolyakov(Gauge_Conf const *const GC,
                   Geometry const *const geo,
                   GParam const *const param,
                   double *repoly,
                   double *impoly);

void polyakov_powers(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     int mu,
                     double *repoly_pwrs,
                     double *impoly_pwrs);

void polyakov_adj(Gauge_Conf const *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  double *repoly,
                  double *impoly);


ACC_ROUTINE_SEQ
double loc_topcharge(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     long r);

void charge_density(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    Meas_Utils *const meas_aux,
                    int const meas_count);

double topcharge(Gauge_Conf const *const GC,
                 Geometry const *const geo,
                 GParam const *const param);

void topcharge_slices(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param,
                      int const mu,
                      double *slices,
                      int const meas_count,
                      FILE *filep);

void topcharge_p_slices(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        int const mu,
                        int const nu,
                        Meas_Utils *const meas_aux,
                        int const meas_count);

double topcharge_prime(Gauge_Conf const *const GC,
                       Geometry const *const geo,
                       GParam const *const param, int const dir);

double topo_chi_prime(Gauge_Conf const *const GC,
                      Geometry const *const geo,
                      GParam const *const param);

void loc_topcharge_corr(Gauge_Conf *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        int const ncool,
                        int const dist,
                        double *res);

void action(Gauge_Conf const *const GC,
            Geometry const *const geo,
            GParam const *const param,
            double *S_wilson, double *S_theta, double *S_total, double *V_mc);

FILE *open_file_with_header_replica(char const *const name,
                                    char const *const header,
                                    int const replica_index,
                                    GParam const *const param,
                                    int const binary_flag);

void open_data_files(Meas_Utils *meas_aux,
                     int const replica_index,
                     GParam const *const param);

void close_data_files(Meas_Utils meas_aux,
                      int const replica_index,
                      GParam const *const param);

void init_meas_utils(Meas_Utils *meas_aux,
                     GParam const *const param,
                     int const replica_index);

void free_meas_utils(Meas_Utils meas_aux,
                     GParam const *const param,
                     int const replica_index);

void init_meas_utils_replica(Meas_Utils **meas_aux,
                             GParam const *const param);

int num_meas_utils_replica(GParam const *const param);

void free_meas_utils_replica(Meas_Utils *meas_aux,
                             GParam const *const param);

void perform_measures_localobs_hot(Gauge_Conf *const GC,
                                   Geometry const *const geo,
                                   GParam const *const param,
                                   Meas_Utils *meas_aux);

void perform_measures_aux(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          int const meas_count,
                          Meas_Utils *meas_aux);

void print_measures_aux(int const num_meas,
                        long const update_index,
                        GParam const *const param,
                        Meas_Utils const *const meas_aux,
                        Smoothing_Type const smoothing);

void perform_measures_localobs_cooling(Gauge_Conf *const GC,
                                       Geometry const *const geo,
                                       GParam const *const param,
                                       Meas_Utils *meas_aux, Cooling_Type type);

void perform_measures_localobs_gradflow(Gauge_Conf *const GC,
                                        Geometry const *const geo,
                                        GParam const *const param,
                                        Meas_Utils *meas_aux);

void perform_measures_localobs_adaptive_gradflow(Gauge_Conf *const GC,
                                                 Geometry const *const geo,
                                                 GParam const *const param,
                                                 Meas_Utils *meas_aux);

void perform_measures_localobs(Gauge_Conf *const GC,
                               Geometry const *const geo,
                               GParam const *const param,
                               Meas_Utils *meas_aux);

void perform_measures_localobs_with_adaptive_gradflow_debug1(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux);

void perform_measures_localobs_with_adaptive_gradflow_debug2(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux);

void perform_measures_localobs_with_adaptive_gradflow_debug3(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux);

void perform_measures_localobs_with_adaptive_gradflow_debug4(Gauge_Conf *const GC,
                                                             Geometry const *const geo,
                                                             GParam const *const param,
                                                             Meas_Utils *meas_aux);

void optimize_multihit_polycorr(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                FILE *datafilep);

void optimize_multilevel_polycorr(Gauge_Conf *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  FILE *datafilep);

void perform_measures_polycorr(Gauge_Conf *const GC,
                               Geometry const *const geo,
                               GParam const *const param,
                               Meas_Utils *meas_aux);

void optimize_multilevel_polycorr_long(Gauge_Conf *const GC,
                                       GParam const *const param,
                                       FILE *datafilep);

void perform_measures_polycorr_long(Gauge_Conf *const GC,
                                    GParam const *const param,
                                    Meas_Utils *meas_aux);

void perform_measures_tube_disc(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                Meas_Utils *meas_aux);

void perform_measures_tube_conn(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                Meas_Utils *meas_aux);

void perform_measures_tube_conn_long(Gauge_Conf *const GC,
                                     GParam const *const param,
                                     Meas_Utils *meas_aux);

void perform_multilevel_update_and_measures(Gauge_Conf *const GC,
                                            Geometry const *const geo,
                                            GParam const *const param,
                                            Meas_Utils *meas_aux,
                                            Multilevel_Obs const ml_obs);


// in gauge_conf_multilevel.c
void multihit(Gauge_Conf const *const GC,
              Geometry const *const geo,
              GParam const *const param,
              long const r,
              int const dir,
              int const num_hit,
              GAUGE_GROUP *G);

void compute_local_poly(Gauge_Conf *const GC,
                        Geometry const *const geo,
                        GParam const *const param);

void update_for_multilevel(Gauge_Conf *const GC,
                           Geometry const *const geo,
                           GParam const *const param,
                           int level);

void multilevel_polycorr(Gauge_Conf *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         int dt);

void multilevel_polycorr_long_zero(Gauge_Conf *const GC,
                                   Geometry const *const geo,
                                   GParam const *const param,
                                   int iteration);

void compute_local_poly_and_plaq(Gauge_Conf *const GC,
                                 Geometry const *const geo,
                                 GParam const *const param);

void multilevel_tube_disc(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          int dt);

void compute_local_poly_plaq_and_plaqconn(Gauge_Conf *const GC,
                                          Geometry const *const geo,
                                          GParam const *const param);

void multilevel_tube_conn(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          int dt);

void multilevel_tube_conn_long_zero(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    int iteration);

void perform_multilevel_long_update_zero(Gauge_Conf *GC,
                                         Geometry const *const geo,
                                         GParam const *const param,
                                         int const iteration,
                                         Multilevel_Obs const ml_obs);

// in gauge_conf_upd.c
ACC_ROUTINE_SEQ
void calcstaples_wilson(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const gparam,
                        long r,
                        int i,
                        GAUGE_GROUP *M);

void calcstaples_tracedef(Gauge_Conf const *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          long r,
                          int i,
                          GAUGE_GROUP *M);

ACC_ROUTINE_SEQ
void calcstaples_with_topo(Gauge_Conf const *const GC,
                           Geometry const *const geo,
                           GParam const *const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M);

void compute_clovers(Gauge_Conf const *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     int const dir);

void compute_clovers_replica(Gauge_Conf const *const GC,
                             Geometry const *const geo,
                             GParam const *const param,
                             int const dir);

void compute_clovers_replica_rect(Gauge_Conf const *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  int const dir,
                                  Rectangle const *const clover_rectangle);

ACC_ROUTINE_SEQ
void calcstaples_wilson_with_defect(Gauge_Conf const *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    long r,
                                    int i,
                                    GAUGE_GROUP *M);

ACC_ROUTINE_SEQ
void calcstaples_with_topo_with_defect(Gauge_Conf const *const GC,
                                       Geometry const *const geo,
                                       GParam const *const param,
                                       long r,
                                       int i,
                                       GAUGE_GROUP *M);


ACC_ROUTINE_SEQ
void heatbath(Gauge_Conf *const GC,
              Geometry const *const geo,
              GParam const *const param,
              long r,
              int i);

ACC_ROUTINE_SEQ
void overrelaxation(Gauge_Conf *const GC,
                    Geometry const *const geo,
                    GParam const *const param,
                    long r,
                    int i);

double delta_action(Gauge_Conf const *const GC,
                    Geometry const *const geo,
                    GParam const *const param);

int metropolis(Gauge_Conf *const GC,
               Geometry const *const geo,
               GParam const *const param,
               long r,
               int i);

int metropolis_with_tracedef(Gauge_Conf *const GC,
                             Geometry const *const geo,
                             GParam const *const param,
                             long r,
                             int i);

ACC_ROUTINE_SEQ
void heatbath_with_defect(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          long r,
                          int i);

ACC_ROUTINE_SEQ
void overrelaxation_with_defect(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param,
                                long r,
                                int i);


void update(Gauge_Conf *const GC,
            Geometry const *const geo,
            GParam const *const param,
            Acc_Utils *acc_counters);

void update_with_trace_def(Gauge_Conf *const GC,
                           Geometry const *const geo,
                           GParam const *const param,
                           double *acc_dt);

void update_with_defect(Gauge_Conf *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        Acc_Utils *acc_counters);

void update_rectangle_with_defect(Gauge_Conf *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param,
                                  int const hierarc_level,
                                  Rect_Utils const *const rect_aux,
                                  Acc_Utils *acc_counters);

void hierarchical_update_rectangle_with_defect(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                                               int const hierarc_level,
                                               Rect_Utils const *const rect_aux,
                                               Acc_Utils *acc_counters);

void parallel_tempering_with_hierarchical_update(Gauge_Conf *const GC, Geometry const *const geo, GParam const *const param,
                                                 Rect_Utils const *const rect_aux,
                                                 Acc_Utils *acc_counters);

void cooling_lex_dir_lexeo_site(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param);

void cooling_lex_site_lex_dir(Gauge_Conf *const GC,
                              Geometry const *const geo,
                              GParam const *const param);

void cooling_lex_dir_lex_site(Gauge_Conf *const GC,
                              Geometry const *const geo,
                              GParam const *const param);

void cooling_lexeo_site_lex_dir(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param);

void cooling_rnd_dir_rndeo_site(Gauge_Conf *const GC,
                                Geometry const *const geo,
                                GParam const *const param);

void cooling(Gauge_Conf *const GC,
             Geometry const *const geo,
             GParam const *const param,
             int const num_steps,
             Cooling_Type const type);

void cooling_with_defect(Gauge_Conf *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         int const n);

void hierarchical_cooling(Gauge_Conf *const GC,
                          Geometry const *const geo,
                          GParam const *const param,
                          Rectangle const *const cooling_rect);

void gradflow_RKstep(Gauge_Conf *const GC,
                     Geometry const *const geo,
                     GParam const *const param,
                     double const dt,
                     Meas_Utils *meas_aux);

double gradflow_RKstep_adaptive_aux(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double const dt,
                                    Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive(Gauge_Conf *const GC,
                              Geometry const *const geo,
                              GParam const *const param,
                              double *t,
                              double *dt,
                              int *accepted,
                              Meas_Utils *meas_aux);

double gradflow_RKstep_adaptive_no_advance(Gauge_Conf *const GC,
                                           Geometry const *const geo,
                                           GParam const *const param,
                                           double const *t,
                                           double const *dt,
                                           double *dt_new,
                                           int *accepted,
                                           Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive_check(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double *t,
                                    double *dt,
                                    int *accepted,
                                    Meas_Utils *meas_aux,
                                    Gauge_Conf *const GC_reset,
                                    Gauge_Conf *const conf_rk3dth);

void gradflow_RKstep_adaptive_debug(Gauge_Conf *const GC,
                                    Geometry const *const geo,
                                    GParam const *const param,
                                    double *t,
                                    double *dt,
                                    int *accepted,
                                    double *total_error,
                                    Meas_Utils *meas_aux);

void gradflow_RKstep_adaptive_debug2(Gauge_Conf *const GC,
                                     Geometry const *const geo,
                                     GParam const *const param,
                                     double *t,
                                     double *dt,
                                     int *accepted,
                                     double *total_error,
                                     Meas_Utils *meas_aux);

void ape_smearing(Gauge_Conf *const GC,
                  Geometry const *const geo,
                  GParam const *const param,
                  double const alpha,
                  int const n);


// in gauge_conf_paral_temp.c
void swap(Gauge_Conf *const GC,
          Geometry const *const geo,
          GParam const *const param,
          Rectangle const *const swap_rectangle,
          Acc_Utils *acc_counters);

void swap_plaquette_sums(Gauge_Conf const *const GC,
                         Geometry const *const geo,
                         GParam const *const param,
                         Rectangle const *const swap_rectangle,
                         double *const P1,
                         double *const P2);

ACC_ROUTINE_SEQ
double swap_action_site(Gauge_Conf const *const GC,
                        Geometry const *const geo,
                        GParam const *const param,
                        long const r,
                        int const i,
                        int const j,
                        double const c_a,
                        double const c_b);

void metropolis_single_swap(Gauge_Conf *const GC,
                            int const a,
                            int const b,
                            double const p,
                            Acc_Utils *acc_counters);

void conf_translation(Gauge_Conf *const GC,
                      Geometry const *const geo,
                      GParam const *const param);

void init_swap_track_file(FILE **swaptrackfilep,
                          GParam const *const param);

void init_acc_utils(Acc_Utils *acc_counters,
                    GParam const *const param);

void free_acc_utils(Acc_Utils *acc_counters,
                    GParam const *const param);

void print_acceptances(Acc_Utils const *const acc_counters,
                       GParam const *const param);

void print_conf_labels(Gauge_Conf const *GC,
                       GParam const *param,
                       Acc_Utils const *acc_counters);

// in gauge_conf_multicanonic.c
void init_multicanonic_gauge_conf(Gauge_Conf *const GC,
                                  Geometry const *const geo,
                                  GParam const *const param);

void init_multicanonic_acc_utils(Acc_Utils *,
                                 GParam const *const);

void free_multicanonic_acc_utils(Acc_Utils *);

void print_multicanonic_acceptance(Gauge_Conf const *const,
                                   GParam const *const,
                                   Acc_Utils const *const);

void init_tune_utils(Tune_Utils *const,
                     GParam const *const);

void free_tune_utils(Tune_Utils *);

void tune_topo_potential_unconstrained(Gauge_Conf const *const,
                                       GParam *const,
                                       Tune_Utils *const,
                                       Acc_Utils const *const);

void tune_topo_potential(Gauge_Conf const *const,
                         GParam *const,
                         Tune_Utils *const,
                         Acc_Utils *const);

int update_tuning_stp(Tune_Utils *const,
                      GParam const *const);

void print_tuning_stp(long,
                      Tune_Utils const *const,
                      GParam const *const);

ACC_ROUTINE_SEQ
double compute_topo_potential(int const,
                              double,
                              GParam const *const);

void compute_topo_potential_debug(char *filename,
                                  double const x_min, double const x_max, double const x_stp,
                                  GParam const *const param);

double multicanonic_topcharge_cooling(Gauge_Conf *const,
                                      Geometry const *const,
                                      GParam const *const);

double multicanonic_topcharge_cooling_rectangle(Gauge_Conf *const,
                                                Geometry const *const,
                                                GParam const *const,
                                                int const,
                                                Rect_Utils const *const);

double topcharge_agf_multicanonic(Gauge_Conf *const,
                                  Geometry const *const,
                                  GParam const *const,
                                  Meas_Utils *);


ACC_ROUTINE_SEQ
double metropolis_prob_multicanonic(int const,
                                    double const,
                                    double const,
                                    GParam const *const);

ACC_ROUTINE_SEQ
double delta_topo_potential_swap(Gauge_Conf const *const,
                                 Gauge_Conf const *const,
                                 GParam const *const);


ACC_ROUTINE_SEQ
int multicanonic_metropolis_step_single_link(Gauge_Conf *const,
                                             Geometry const *const,
                                             GParam const *const,
                                             long, int,
                                             GAUGE_GROUP const *const);

int multicanonic_metropolis_step_all_links(Gauge_Conf *const,
                                           Geometry const *const,
                                           GParam const *const);

int multicanonic_metropolis_step_rectangle(Gauge_Conf *const,
                                           Geometry const *const,
                                           GParam const *const,
                                           int const,
                                           Rect_Utils const *const);

void multicanonic_metropolis_replica(Gauge_Conf *const GC,
                                     Geometry const *const geo,
                                     GParam const *const param,
                                     int *const acc);

void multicanonic_metropolis_rectangle_replica(Gauge_Conf *const GC,
                                               Geometry const *const geo,
                                               GParam const *const param,
                                               int const hierarc_level,
                                               Rect_Utils const *const rect_aux,
                                               int *const acc);

// old implementation of multicanonic

ACC_ROUTINE_SEQ
void compute_single_clover_insertion(GAUGE_GROUP *,
                                     Gauge_Conf const *const,
                                     Geometry const *const,
                                     GParam const *const,
                                     long, int, int);

ACC_ROUTINE_SEQ
void compute_topostaple_alone(Gauge_Conf const *const,
                              Geometry const *const,
                              GParam const *const,
                              long, int,
                              GAUGE_GROUP *);

ACC_ROUTINE_SEQ
double delta_Q_upd(Gauge_Conf const *const,
                   Geometry const *const,
                   GParam const *const,
                   long, int,
                   GAUGE_GROUP const *const);

void multicanonic_parallel_tempering_with_hierarchical_update(Gauge_Conf *,
                                                              Geometry const *const,
                                                              GParam const *const,
                                                              Rectangle const *const,
                                                              Rectangle const *const,
                                                              Rectangle const *const,
                                                              Acc_Utils *);

void multicanonic_agf_parallel_tempering_with_hierarchical_update(Gauge_Conf *,
                                                                  Geometry const *const,
                                                                  GParam const *const,
                                                                  Rectangle const *const,
                                                                  Rectangle const *const,
                                                                  Rectangle const *const,
                                                                  Acc_Utils *);

void multicanonic_update_with_defect(Gauge_Conf *,
                                     Geometry const *const,
                                     GParam const *const,
                                     Acc_Utils *);

void multicanonic_agf_update_with_defect(Gauge_Conf *,
                                         Geometry const *const,
                                         GParam const *const,
                                         Acc_Utils *);

void multicanonic_update_rectangle_with_defect(Gauge_Conf *,
                                               Geometry const *const,
                                               GParam const *const,
                                               Rectangle const *const,
                                               Rectangle const *const,
                                               Acc_Utils *);

void multicanonic_agf_update_rectangle_with_defect(Gauge_Conf *,
                                                   Geometry const *const,
                                                   GParam const *const,
                                                   Rectangle const *const,
                                                   Rectangle const *const,
                                                   Acc_Utils *);

void multicanonic_hierarchical_update_rectangle_with_defect(Gauge_Conf *,
                                                            Geometry const *const,
                                                            GParam const *const,
                                                            int const,
                                                            Rectangle const *const,
                                                            Rectangle const *const,
                                                            Rectangle const *const,
                                                            Acc_Utils *);

void multicanonic_agf_hierarchical_update_rectangle_with_defect(Gauge_Conf *,
                                                                Geometry const *const,
                                                                GParam const *const,
                                                                int const,
                                                                Rectangle const *const,
                                                                Rectangle const *const,
                                                                Rectangle const *const,
                                                                Acc_Utils *);

ACC_ROUTINE_SEQ
int multicanonic_overrelaxation_with_defect(Gauge_Conf *,
                                            Geometry const *const,
                                            GParam const *const,
                                            long, int);

ACC_ROUTINE_SEQ
int multicanonic_heatbath_with_defect(Gauge_Conf *,
                                      Geometry const *const,
                                      GParam const *const,
                                      long, int);

#endif
