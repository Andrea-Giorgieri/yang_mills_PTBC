#ifndef GPARAM_H
#define GPARAM_H

#include"macro.h"
#include"timing.h"

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>

// cooling implementations
typedef enum
	{
	LEX_SITE_LEX_DIR,
	LEX_DIR_LEX_SITE,
	LEX_DIR_LEXEO_SITE,
	LEXEO_SITE_LEX_DIR,
	RND_DIR_RNDEO_SITE
	} Cooling_Type;

typedef struct GParam
	{
	// lattice dimensions
	int d_size[STDIM];

	// simulation parameters
	double d_beta;
	double d_h[NCOLOR]; // parameters for the trace deformation
	double d_theta;

	// parallel tempering parameters
	int d_defect_dir;              // defect boundary
	int d_L_defect[STDIM - 1];     // defect sizes
	int d_N_replica_pt;            // numbers of replica used in parallel tempering
	double *d_pt_bound_cond_coeff; // boundary conditions coefficients

	// twist parameters and open boundary conditions
	int d_k_twist[STDIM * (STDIM - 1) / 2]; // twist parameter for each plane
	int d_obc_dir;                          // direction of obc (-1 if pbc)
	int d_obc_default_pos;                  // default starting position of obc
	int d_obc_bulk;                         // size of the bulk along d_obc_dir

	// hierarchical update (parallel tempering)
	int d_N_hierarc_levels; // number of hierarchical levels
	int *d_L_rect;          // d_L_rect is a vector of length d_N_hierarc_levels, d_L_rect[i] is the extension of the rectangle at the i-th hierarchical level
	int *d_N_sweep_rect;    // d_N_sweep_rect is vector of length d_N_hierarch_levels, d_N_sweep_rect[i] is the number of sweep of the rectangle at the i-th hierarchical level

	// simulation details
	int d_sample;
	int d_thermal;
	int d_overrelax;
	int d_measevery;

	// time limit in hours
	double d_walltime;

	// initialization & saving
	int d_start;
	int d_saveconf_back_every;
	int d_saveconf_analysis_every;

	// for metropolis
	double d_epsilon_metro;

	// for observables to measure
	int d_plaquette_meas;
	int d_clover_energy_meas;
	int d_energy_density_meas;
	int d_charge_meas;
	int d_charge_density_meas;
	int d_polyakov_meas;
	int d_polyakov_powers_meas;
	int d_polyakov_density_meas;
	int d_chi_prime_meas;
	int d_charge_prime_meas;
	int d_action_meas;
	int d_energy_slices_meas;
	int d_charge_slices_meas;
	int d_charge_p_slices_meas;

	int d_multipolyakov_order;
	int *d_multipolyakov_dirs;
	int d_meas_effective_charge;

	// for cooling in measures
	char d_cooling_type_str[STD_STRING_LENGTH];
	Cooling_Type d_cooling_type;
	int d_coolsteps;
	int d_coolrepeat;

	// for gradient-flow evolution
	double d_gfstep;
	int d_ngfsteps;
	int d_gf_meas_each;
	int d_gf_num_meas;

	// for adaptive-step gradient-flow evolution
	double d_agf_length;
	double d_agf_meas_each;
	double d_agf_step;
	double d_agf_delta;
	double d_agf_time_bin;
	int d_agf_num_meas;

	// for multilevel
	int d_multihit;
	int d_ml_step[NLEVELS];
	int d_ml_upd[NLEVELS];
	int d_ml_level0_repeat;
	int d_dist_poly;
	int d_trasv_dist;
	int d_plaq_dir[2];
	int d_ml_num_hit;
	int d_ml_num_slices[NLEVELS];

	// output file names
	char d_conf_file[STD_STRING_LENGTH];             // save gauge configuration
	char d_twist_file[STD_STRING_LENGTH];            // save twist configuration
	char d_data_file[STD_STRING_LENGTH];             // print measures of simple observables
	char d_chiprime_file[STD_STRING_LENGTH];         // print chi prime measures
	char d_energy_slices_file[STD_STRING_LENGTH];    // print energy slices measures
	char d_charge_slices_file[STD_STRING_LENGTH];    // print topological charge time correlator measures
	char d_energydensity_file[STD_STRING_LENGTH];    // print energy density measures
	char d_chargedensity_file[STD_STRING_LENGTH];    // print charge density measures
	char d_polyakovdensity_file[STD_STRING_LENGTH];  // print polyakov density measures
	char d_log_file[STD_STRING_LENGTH];              // print program details
	char d_ml_file[STD_STRING_LENGTH];               //
	char d_swap_acc_file[STD_STRING_LENGTH];         // print swap Metropolis acceptance
	char d_swap_tracking_file[STD_STRING_LENGTH];    // print swap tracks
	char d_multicanonic_acc_file[STD_STRING_LENGTH]; // print multicanonic Metropolis acceptance

	// random seed
	unsigned int d_randseed;

	// derived constants
	int d_max_size;                // max lattice size
	int d_min_size;                // min lattice size
	int d_even_size[STDIM];        // sizes of the largest even sublattice
	long d_volume;                 // total volume
	long d_even_volume;            // total volume of larges even sublattice
	double d_inv_vol;              // 1 / total volume
	long d_space_vol[STDIM];       // volume without given component
	double d_inv_space_vol[STDIM]; // 1 / volume without given component
	long d_volume_defect;          // volume of the defect (only for parallel tempering)
	int d_n_grid;                  // total grid points (only for multicanonic)
	int d_n_planes;                // number of planes (only for twisted boundary conditions)
	long d_n_even;                 // number of even lattice sites in the largest even sublattice
	long d_n_border;               // number of lattice sites outside the largest even sublattice

	// for multicanonic
	char d_topo_potential_file[STD_STRING_LENGTH];
	double d_grid_step;
	double d_grid_max;
	double **d_grid;        // d_grid [a][x] = V_a(x) is the topo potential for replica a
	int d_topo_cooling;     // cooling strat for topcharge before evaluating V_a: 0 = none, 1 = cooling
	int d_topo_coolsteps;   // cooling steps of the charge before evaluating V_a
	double d_topo_agf_time; // adaptive gradflow time before evaluating V_a. TODO: unused, debug only, remove?
	double d_topo_alpha;    // alpha parameter for alpha-rounding of cooled charge

	// for tuning
	int d_topo_tuning_even;       // force V_a(x) to be even during tuning (0 = False, 1 = True)
	int d_topo_tuning_save_every; // save tuned V_a(x) every this number of steps (0 = Never)
	double d_topo_tuning_thr;     // threshold below which tuning of topo potential is completed
	double d_topo_tuning_stp;     // initial variation of topo potential during tuning

	// for debugging and testing
	int d_test_flag;
	} GParam;

// functions to impose conditions on params
int param_any_ui(unsigned int val, char *msg);

int param_any_int(int val, char *msg);

int param_bool_int(int val, char *msg);

int param_positive_int(int val, char *msg);

int param_nonnegative_int(int val, char *msg);

int param_any_double(double val, char *msg);

int param_positive_double(double val, char *msg);

int param_nonnegative_double(double val, char *msg);

int param_any_string(char *val, char *msg);

// functions to check that a required param was found
void check_required_string(char *val, char *name, int required);

// functions to set values of params
void set_ui_param(FILE *fp, unsigned int *ptr, char const *const name, int (*condition)(unsigned int, char *));

void set_int_param(FILE *fp, int *ptr, char const *const name, int (*condition)(int, char *));

void set_double_param(FILE *fp, double *ptr, char const *const name, int (*condition)(double, char *));

void set_string_param(FILE *fp, char *ptr, char const *const name, int (*condition)(char *, char *));

void set_defaults(GParam *const param);

// read and write parameters
void readinput(char const *const in_file, GParam *const param);

void remove_white_line_and_comments(FILE *input);

void read_topo_potential(GParam *param);

void write_topo_potential(GParam const *const param, char *filename);

void init_derived_constants(GParam *const param);

void edit_params_meas_only(GParam *const param);

void free_hierarc_params(GParam *const param);

// print simulation parameters aux
void print_configuration_parameters(FILE *fp);

void print_pt_parameters(FILE *fp, GParam const *const param);

void print_multicanonic_parameters(FILE *fp, GParam const *const param);

void print_multicanonic_tuning_parameters(FILE *fp, GParam const *const param);

void print_simul_parameters(FILE *fp, GParam const *const param);

void print_smoothing_parameters(FILE *fp, GParam const *const param);

void print_multilevel_parameters(FILE *fp, GParam const *const param);

void print_metro_parameters(FILE *fp, GParam const *const param, double acc);

// print simulation parameters
void print_parameters_local(GParam const *const param, Time_Utils const *const timers);

void print_parameters_local_agf(GParam const *const param, Time_Utils const *const timers);

void print_parameters_local_pt(GParam const *const param, Time_Utils const *const timers);

void print_parameters_local_pt_gf(GParam const *const param, Time_Utils const *const timers);

void print_parameters_local_pt_agf(GParam const *const param, Time_Utils const *const timers);

void print_parameters_debug_agf_vs_gf(GParam const *const param, time_t time_start, time_t time_end, time_t agf_time, time_t dagf_time, time_t gf_time);

void print_parameters_debug_agf_vs_delta(GParam const *const param, time_t time_mc, time_t time_agf0, time_t time_agf1, time_t time_agf2, time_t time_agf3);

void print_parameters_local_pt_multicanonic(GParam const *const param, Time_Utils const *const timers);

void print_parameters_polycorr_long(GParam *param, Time_Utils const *const timers);

void print_parameters_polycorr(GParam *param, Time_Utils const *const timers);

void print_parameters_t0(GParam *param, Time_Utils const *const timers);

void print_parameters_gf(GParam *param, Time_Utils const *const timers);

void print_parameters_agf(GParam *param, Time_Utils const *const timers);

void print_parameters_tracedef(GParam const *const param, Time_Utils const *const timers, double acc);

void print_parameters_tube_disc(GParam *param, Time_Utils const *const timers);

void print_parameters_tube_conn(GParam *param, Time_Utils const *const timers);

void print_parameters_tube_conn_long(GParam *param, Time_Utils const *const timers);

void print_parameters_tuning_pt_mc(GParam const *const param, Time_Utils const *const timers, int count);

// print template input aux
void print_template_volume_parameters(FILE *fp);

void print_template_simul_parameters(FILE *fp);

void print_template_pt_parameters(FILE *fp);

void print_template_twist_parameters(FILE *fp);

void print_template_adaptive_gradflow_parameters(FILE *fp);

void print_template_gradflow_parameters(FILE *fp);

void print_template_cooling_parameters(FILE *fp);

void print_template_metro_parameters(FILE *fp);

void print_template_multicanonic_parameters(FILE *fp);

void print_template_multicanonic_tuning_parameters(FILE *fp);

void print_template_multilevel_parameters(FILE *fp);

void print_template_output_parameters(FILE *fp);

// print program details
void print_authors(int parallel_tempering, int twisted_bc);

void print_compilation_details(void);

#endif
