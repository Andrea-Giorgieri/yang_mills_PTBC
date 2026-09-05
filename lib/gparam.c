#ifndef GPARAM_C
#define GPARAM_C

#include "../include/macro.h"
#include "../include/gparam.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/endianness.h"
#include "../include/memalign.h"
#include "../include/timing.h"


// functions to impose conditions on params

static inline char const *param_any_ui(unsigned int const val)
	{
	(void) val;
	return NULL;
	}


static inline char const *param_any_int(int const val)
	{
	(void) val;
	return NULL;
	}


static inline char const *param_bool_int(int const val)
	{
	if((val == 0) || (val == 1)) return NULL;
	return "must be either 0 or 1";
	}


static inline char const *param_positive_int(int const val)
	{
	if(val > 0) return NULL;
	return "must be positive";
	}


static inline char const *param_nonnegative_int(int const val)
	{
	if(val >= 0) return NULL;
	return "must be non-negative";
	}


static inline char const *param_any_double(double const val)
	{
	(void) val;
	return NULL;
	}


static inline char const *param_positive_double(double const val)
	{
	if(val > 0) return NULL;
	return "must be positive";
	}


static inline char const *param_nonnegative_double(double const val)
	{
	if(val >= 0) return NULL;
	return "must be non-negative";
	}


static inline char const *param_any_string(char const *val)
	{
	(void) val;
	return NULL;
	}


static inline void check_required_string(char const *val, char const *name, int const required)
	{
	if(required)
		{
		REQUIRE(val != NULL && val[0] != '\0', "parameter '%s' is required", name);
		}
	}


// functions to set values of params

static inline void set_ui_param(FILE *fp, unsigned int *ptr, char const *const name, char const *(*condition)(unsigned int))
	{
	unsigned int temp;
	REQUIRE(fscanf(fp, "%u", &temp) == 1, "error reading parameter '%s' from input file", name);
	char const *msg = condition(temp);
	REQUIRE(msg == NULL, "invalid parameter '%s': %s", name, msg);
	*ptr = temp;
	}


static inline void set_int_param(FILE *fp, int *ptr, char const *const name, char const *(*condition)(int))
	{
	int temp;
	REQUIRE(fscanf(fp, "%d", &temp) == 1, "error reading parameter '%s' from input file", name);
	char const *msg = condition(temp);
	REQUIRE(msg == NULL, "invalid parameter '%s': %s", name, msg);
	*ptr = temp;
	}


static inline void set_double_param(FILE *fp, double *ptr, char const *const name, char const *(*condition)(double))
	{
	double temp;
	REQUIRE(fscanf(fp, "%lf", &temp) == 1, "error reading parameter '%s' from input file", name);
	char const *msg = condition(temp);
	REQUIRE(msg == NULL, "invalid parameter '%s': %s", name, msg);
	*ptr = temp;
	}


static inline void set_string_param(FILE *fp, char *ptr, char const *const name, char const *(*condition)(char const *))
	{
	char temp[STD_STRING_LENGTH];
	REQUIRE(fscanf(fp, "%s", temp) == 1, "error reading parameter '%s' from input file", name);
	char const *msg = condition(temp);
	REQUIRE(msg == NULL, "invalid parameter '%s': %s", name, msg);
	strcpy(ptr, temp);
	}


static int set_parameter(FILE *input, ParamDef const *param)
	{
	switch(param->type)
		{
		case PARAM_INT:
			set_int_param(input, (int *) param->ptr, param->name, param->condition.integer);
			break;

		case PARAM_UINT:
			set_ui_param(input, (unsigned int *) param->ptr, param->name, param->condition.ui);
			break;

		case PARAM_DOUBLE:
			set_double_param(input, (double *) param->ptr, param->name, param->condition.dbl);
			break;

		case PARAM_STRING:
			set_string_param(input, (char *) param->ptr, param->name, param->condition.string);
			break;

		case PARAM_INT_ARRAY:
			for(int i = 0; i < param->count; i++)
				{
				set_int_param(input, &((int *) param->ptr)[i], param->name, param->condition.integer);
				}
			break;

		case PARAM_DOUBLE_ARRAY:
			for(int i = 0; i < param->count; i++)
				{
				set_double_param(input, &((double *) param->ptr)[i], param->name, param->condition.dbl);
				}
			break;

		default:
			REQUIRE(0, "unknown type for parameter '%s'", param->name);
		}
	return 1;
	}


static int find_and_set_parameter(FILE *input, char const *name, ParamDef const *param_defs, size_t const num_params)
	{
	for(size_t i = 0; i < num_params; i++)
		{
		if(strcmp(name, param_defs[i].name) == 0)
			return set_parameter(input, &param_defs[i]);
		}
	return 0;
	}


static void read_n_replica_pt(FILE *input, GParam *param)
	{
	set_int_param(input, &param->d_N_replica_pt, "N_replica_pt", &param_positive_int);
	allocate_array_double(&param->d_pt_bound_cond_coeff, param->d_N_replica_pt, __FILE__, __LINE__);
	for(int i = 0; i < param->d_N_replica_pt; ++i)
		{
		set_double_param(input, &param->d_pt_bound_cond_coeff[i], "N_replica_pt", &param_any_double);
		}
	}


static void read_hierarc_upd(FILE *input, GParam *param)
	{
	set_int_param(input, &param->d_N_hierarc_levels, "hierarc_upd", &param_nonnegative_int);
	if(param->d_N_hierarc_levels == 0)
		return;
	allocate_array_int(&param->d_L_rect, param->d_N_hierarc_levels, __FILE__, __LINE__);
	allocate_array_int(&param->d_N_sweep_rect, param->d_N_hierarc_levels, __FILE__, __LINE__);
	for(int i = 0; i < param->d_N_hierarc_levels; ++i)
		{
		set_int_param(input, &param->d_L_rect[i], "hierarc_upd", &param_nonnegative_int);
		}
	for(int i = 0; i < param->d_N_hierarc_levels; ++i)
		{
		set_int_param(input, &param->d_N_sweep_rect[i], "hierarc_upd", &param_nonnegative_int);
		}
	}


static void read_multipolyakov_order(FILE *input, GParam *param)
	{
	set_int_param(input, &param->d_multipolyakov_order, "multipolyakov_order", &param_nonnegative_int);
	if(param->d_multipolyakov_order == 0)
		return;
	allocate_array_int(&param->d_multipolyakov_dirs, param->d_multipolyakov_order, __FILE__, __LINE__);
	for(int i = 0; i < param->d_multipolyakov_order; ++i)
		{
		set_int_param(input, &param->d_multipolyakov_dirs[i], "multipolyakov_order", &param_nonnegative_int);
		}
	}


void set_defaults(GParam *const param)
	{
	int i;

	// multilevel
	strcpy(param->d_ml_obs_str, "NONE");
	for(i = 0; i < NLEVELS; i++)
		{
		param->d_ml_step[i] = NLEVELS - i;
		param->d_ml_upd[i] = 0;
		}
	param->d_multihit = 0;
	param->d_ml_level0_repeat = 0;
	param->d_dist_poly = 0;
	param->d_transv_dist = 0;
	param->d_plaq_dir[0] = 0;
	param->d_plaq_dir[1] = 1;

	// trace deformation and theta term
	for(i = 0; i < NCOLOR / 2; i++) param->d_h[i] = 0.0;
	param->d_theta = 0.0;
	param->d_theta_profile_dir = -1;

	// twist and open boundary conditions
	for(i = 0; i < STDIM * (STDIM - 1) / 2; i++) param->d_k_twist[i] = 0;
	param->d_obc_dir = -1;
	param->d_obc_bulk = -1;

	// parallel tempering
	for(i = 0; i < STDIM - 1; i++) param->d_L_defect[i] = 0;
	param->d_defect_dir = 0;
	param->d_N_replica_pt = 1;
	allocate_array_double(&(param->d_pt_bound_cond_coeff), 1, __FILE__, __LINE__);
	param->d_pt_bound_cond_coeff[0] = 1.0;

	// gradient flow
	param->d_ngfsteps = 0;
	param->d_gf_meas_each = 1;
	param->d_gfstep = 1;

	// adaptive gradient flow
	param->d_agf_length = 0;
	param->d_agf_meas_each = 1;
	param->d_agf_step = 0.1;
	param->d_agf_delta = 0.00001;
	param->d_agf_time_bin = 0;

	// cooling
	strcpy(param->d_cooling_type_str, "LEX_DIR_LEXEO_SITE");
	param->d_coolrepeat = 0;
	param->d_coolsteps = 1;

	// multicanonical with cold topcharge
	param->d_topo_cooling = 0;
	param->d_topo_coolsteps = 0;
	param->d_topo_alpha = 0;
	param->d_topo_tuning_thr = 0.1;
	param->d_topo_tuning_stp = 0.1;
	param->d_topo_tuning_save_every = 0;
	param->d_topo_tuning_even = 1;

	// walltime
	param->d_walltime = 24.0;

	// measurements
	param->d_plaquette_meas = 1;
	param->d_clover_energy_meas = 0;
	param->d_energy_density_meas = 0;
	param->d_charge_meas = 0;
	param->d_charge_density_meas = 0;
	param->d_chi_prime_meas = 0;
	param->d_charge_prime_meas = 0;
	param->d_polyakov_meas = 0;
	param->d_polyakov_powers_meas = 0;
	param->d_polyakov_density_meas = 0;
	param->d_energy_slices_meas = 0;
	param->d_charge_slices_meas = 0;
	param->d_charge_p_slices_meas = 0;
	param->d_action_meas = 0;

	param->d_multipolyakov_order = 0;
	param->d_meas_effective_charge = 0;

	// filenames
	strcpy(param->d_conf_file, "");
	strcpy(param->d_twist_file, "");
	strcpy(param->d_data_file, "");
	strcpy(param->d_energydensity_file, "");
	strcpy(param->d_chargedensity_file, "");
	strcpy(param->d_polyakovdensity_file, "");
	strcpy(param->d_chiprime_file, "");
	strcpy(param->d_energy_slices_file, "");
	strcpy(param->d_charge_slices_file, "");
	strcpy(param->d_log_file, "");
	strcpy(param->d_swap_acc_file, "");
	strcpy(param->d_swap_tracking_file, "");
	strcpy(param->d_theta_profile_file, "");
	strcpy(param->d_multicanonic_acc_file, "");
	strcpy(param->d_topo_potential_file, "");
	strcpy(param->d_ml_file, "");
	}


// remove white lines and comments starting with # from input file
void remove_white_line_and_comments(FILE *input)
	{
	int temp_i;

	// skip empty line
	{
		do temp_i = getc(input); while(temp_i == '\n' || temp_i == ' ');
	}

	// skip comment, from \043 = ascii oct for # to first newline or EOF
	if(temp_i == '\043') { do temp_i = getc(input); while(temp_i != '\n' && temp_i != EOF); }

	// return if EOF reached or nothing else to remove, after pushing back last char
	if(temp_i == EOF) return;
	ungetc(temp_i, input);
	if(temp_i != '\n' && temp_i != ' ' && temp_i != '\043') return;

	// recursive call if another white line or comment
	remove_white_line_and_comments(input);
	}


// read params from input file
void readinput(char const *const in_file, GParam *const param)
	{
	char str[STD_STRING_LENGTH];
	int i, err;

	#define INT_PARAM(name, member, condition) { name, PARAM_INT, &(param->member), 1, {.integer = condition} }
	#define UINT_PARAM(name, member, condition) { name, PARAM_UINT, &(param->member), 1, {.ui = condition} }
	#define DOUBLE_PARAM(name, member, condition) { name, PARAM_DOUBLE, &(param->member), 1, {.dbl = condition} }
	#define STRING_PARAM(name, member, condition) { name, PARAM_STRING, param->member, 1, {.string = condition} }
	#define INT_ARRAY_PARAM(name, member, count, condition) { name, PARAM_INT_ARRAY, param->member, count, {.integer = condition} }
	#define DOUBLE_ARRAY_PARAM(name, member, count, condition) { name, PARAM_DOUBLE_ARRAY, param->member, count, {.dbl = condition} }

	ParamDef const params[] =
			{
			// lattice
			INT_ARRAY_PARAM("size", d_size, STDIM, param_positive_int),
			DOUBLE_PARAM("beta", d_beta, param_positive_double),
			DOUBLE_ARRAY_PARAM("htracedef", d_h, NCOLOR / 2, param_any_double),
			DOUBLE_PARAM("theta", d_theta, param_any_double),
			INT_PARAM("theta_profile_dir", d_theta_profile_dir, param_any_int),
			STRING_PARAM("theta_profile_file", d_theta_profile_file, param_any_string),
			INT_ARRAY_PARAM("k_twist", d_k_twist, STDIM * (STDIM - 1) / 2, param_any_int),
			INT_PARAM("obc_dir", d_obc_dir, param_any_int),
			INT_PARAM("obc_bulk", d_obc_bulk, param_nonnegative_int),

			// parallel tempering
			INT_PARAM("defect_dir", d_defect_dir, param_nonnegative_int),
			INT_ARRAY_PARAM("defect_size", d_L_defect, STDIM - 1, param_any_int),

			// simulation
			INT_PARAM("sample", d_sample, param_nonnegative_int),
			INT_PARAM("thermal", d_thermal, param_nonnegative_int),
			INT_PARAM("overrelax", d_overrelax, param_nonnegative_int),
			INT_PARAM("measevery", d_measevery, param_nonnegative_int),
			INT_PARAM("start", d_start, param_nonnegative_int),
			INT_PARAM("saveconf_back_every", d_saveconf_back_every, param_nonnegative_int),
			INT_PARAM("saveconf_analysis_every", d_saveconf_analysis_every, param_nonnegative_int),
			UINT_PARAM("randseed", d_randseed, param_any_ui),
			DOUBLE_PARAM("walltime", d_walltime, param_any_double),

			// measure flags
			INT_PARAM("plaquette_meas", d_plaquette_meas, param_bool_int),
			INT_PARAM("clover_energy_meas", d_clover_energy_meas, param_bool_int),
			INT_PARAM("energy_density_meas", d_energy_density_meas, param_bool_int),
			INT_PARAM("charge_meas", d_charge_meas, param_bool_int),
			INT_PARAM("charge_density_meas", d_charge_density_meas, param_bool_int),
			INT_PARAM("polyakov_meas", d_polyakov_meas, param_bool_int),
			INT_PARAM("polyakov_powers_meas", d_polyakov_powers_meas, param_bool_int),
			INT_PARAM("polyakov_density_meas", d_polyakov_density_meas, param_bool_int),
			INT_PARAM("chi_prime_meas", d_chi_prime_meas, param_bool_int),
			INT_PARAM("charge_prime_meas", d_charge_prime_meas, param_bool_int),
			INT_PARAM("action_meas", d_action_meas, param_bool_int),
			INT_PARAM("energy_slices_meas", d_energy_slices_meas, param_bool_int),
			INT_PARAM("charge_slices_meas", d_charge_slices_meas, param_bool_int),
			INT_PARAM("charge_p_slices_meas", d_charge_p_slices_meas, param_bool_int),
			INT_PARAM("meas_effective_charge", d_meas_effective_charge, param_bool_int),

			// cooling
			STRING_PARAM("cooling_type", d_cooling_type_str, param_any_string),
			INT_PARAM("coolsteps", d_coolsteps, param_nonnegative_int),
			INT_PARAM("coolrepeat", d_coolrepeat, param_nonnegative_int),

			// gradient flow
			DOUBLE_PARAM("gfstep", d_gfstep, param_nonnegative_double),
			INT_PARAM("num_gfsteps", d_ngfsteps, param_nonnegative_int),
			INT_PARAM("gf_meas_each", d_gf_meas_each, param_positive_int),

			// adaptive gradient flow
			DOUBLE_PARAM("agf_length", d_agf_length, param_nonnegative_double),
			DOUBLE_PARAM("agf_meas_each", d_agf_meas_each, param_nonnegative_double),
			DOUBLE_PARAM("agf_step", d_agf_step, param_positive_double),
			DOUBLE_PARAM("agf_delta", d_agf_delta, param_positive_double),
			DOUBLE_PARAM("agf_time_bin", d_agf_time_bin, param_nonnegative_double),

			// multilevel
			INT_PARAM("multihit", d_multihit, param_any_int),
			STRING_PARAM("ml_obs", d_ml_obs_str, param_any_string),
			INT_ARRAY_PARAM("ml_step", d_ml_step, NLEVELS, param_any_int),
			INT_ARRAY_PARAM("ml_upd", d_ml_upd, NLEVELS, param_any_int),
			INT_PARAM("ml_level0_repeat", d_ml_level0_repeat, param_any_int),
			STRING_PARAM("ml_file", d_ml_file, param_any_string),
			INT_PARAM("dist_poly", d_dist_poly, param_any_int),
			INT_PARAM("transv_dist", d_transv_dist, param_any_int),
			INT_ARRAY_PARAM("plaq_dir", d_plaq_dir, 2, param_any_int),

			// configuration filenames
			STRING_PARAM("conf_file", d_conf_file, param_any_string),
			STRING_PARAM("twist_file", d_twist_file, param_any_string),

			// data filenames
			STRING_PARAM("data_file", d_data_file, param_any_string),
			STRING_PARAM("energy_density_file", d_energydensity_file, param_any_string),
			STRING_PARAM("charge_density_file", d_chargedensity_file, param_any_string),
			STRING_PARAM("polyakov_density_file", d_polyakovdensity_file, param_any_string),
			STRING_PARAM("chiprime_data_file", d_chiprime_file, param_any_string),
			STRING_PARAM("energy_slices_file", d_energy_slices_file, param_any_string),
			STRING_PARAM("charge_slices_file", d_charge_slices_file, param_any_string),

			// log, acceptance and tracking
			STRING_PARAM("log_file", d_log_file, param_any_string),
			STRING_PARAM("swap_acc_file", d_swap_acc_file, param_any_string),
			STRING_PARAM("swap_track_file", d_swap_tracking_file, param_any_string),
			STRING_PARAM("multicanonic_acc_file", d_multicanonic_acc_file, param_any_string),

			// multicanonical
			STRING_PARAM("topo_potential_file", d_topo_potential_file, param_any_string),
			DOUBLE_PARAM("grid_step", d_grid_step, param_positive_double),
			DOUBLE_PARAM("grid_max", d_grid_max, param_positive_double),
			INT_PARAM("topo_cooling", d_topo_cooling, param_any_int),
			INT_PARAM("topo_coolsteps", d_topo_coolsteps, param_nonnegative_int),
			DOUBLE_PARAM("topo_alpha", d_topo_alpha, param_nonnegative_double),

			// multicanonical tuning
			DOUBLE_PARAM("topo_tuning_thr", d_topo_tuning_thr, param_nonnegative_double),
			DOUBLE_PARAM("topo_tuning_stp", d_topo_tuning_stp, param_positive_double),
			INT_PARAM("topo_tuning_save_every", d_topo_tuning_save_every, param_nonnegative_int),
			INT_PARAM("topo_tuning_even", d_topo_tuning_even, param_bool_int),

			// other
			DOUBLE_PARAM("epsilon_metro", d_epsilon_metro, param_nonnegative_double),
			INT_PARAM("test_flag", d_test_flag, param_any_int)
			};

	size_t const num_params = sizeof(params) / sizeof(params[0]);

	#undef INT_PARAM
	#undef UINT_PARAM
	#undef DOUBLE_PARAM
	#undef STRING_PARAM
	#undef INT_ARRAY_PARAM
	#undef DOUBLE_ARRAY_PARAM

	// set default values
	set_defaults(param);

	// open the input file
	FILE *input = fopen(in_file, "r");
	REQUIRE(input != NULL, "failed to open input file %s", in_file);

	// slide the input file
	for(;;)
		{
		// check for EOF after skipping white line or comments
		remove_white_line_and_comments(input);
		i = getc(input);
		if(i == EOF) break;
		ungetc(i, input);

		// read the name of a parameter in str
		err = fscanf(input, "%s", str);
		REQUIRE(err == 1, "error reading the name of a parameter from input file");

		// ordinary parameters
		if(find_and_set_parameter(input, str, params, num_params))
			continue;

		// special parameters that require custom reading functions
		if(strcmp(str, "N_replica_pt") == 0)
			{
			read_n_replica_pt(input, param);
			continue;
			}

		if(strcmp(str, "hierarc_upd") == 0)
			{
			read_hierarc_upd(input, param);
			continue;
			}

		if(strcmp(str, "multipolyakov_order") == 0)
			{
			read_multipolyakov_order(input, param);
			continue;
			}

		// unknown parameter
		REQUIRE(0, "unrecognized parameter '%s' in input file %s", str, in_file);
		}

	// close the input file
	fclose(input);

	// convert walltime from hours to seconds
	param->d_walltime *= 3600;

	// Further checks

	// multilevel
	#if NLEVELS > 1
	check_required_string(param->d_ml_file, "ml_file", 1);

	REQUIRE(param->d_size[0] % param->d_ml_step[0] == 0 &&
	        param->d_size[0] >= param->d_ml_step[0],
	        "size[0] must be divisible by ml_step[0] and >= ml_step[0]");

	for(i = 1; i < NLEVELS; i++)
		{
		REQUIRE(param->d_ml_step[i - 1] % param->d_ml_step[i] == 0 &&
		        param->d_ml_step[i - 1] > param->d_ml_step[i],
		        "ml_step[%d] must divide ml_step[%d] and be smaller",
		        i, i - 1);
		}

	REQUIRE(param->d_ml_step[NLEVELS - 1] > 1, "ml_step[%d] must be > 1", NLEVELS - 1);
	#endif

	// Along odd sides L_mu, x_mu = 0 and x_mu = L_mu-1 are neighbors but even. This prevents naive even-odd parallelization of updates.
	// This is dealt with only when sweeping the full lattice, not when sweeping the spatial lattice as required by trace deformation.
	#ifndef OPENMP_MODE
	for(i = 0; i < NCOLOR / 2; i++)
		{
		if(param->d_h[i] != 0.0)
			{
			for(int j = 1; j < STDIM; j++)
				{
				REQUIRE(param->d_size[j] % 2 == 0, "when using OpenMP and trace deformation, spatial lattice sizes must be even");
				}
			break;
			}
		}
	#endif

	// sizes
	for(i = 0; i < STDIM; i++)
		{
		REQUIRE(param->d_size[i] > 1, "all lattice sizes must be larger than 1");
		}

	// open boundary conditions
	REQUIRE(param->d_obc_dir >= -1 && param->d_obc_dir < STDIM,
	        "direction of open boundary conditions must be -1 (PBC) or in [0, %d)",
	        STDIM);
	if(param->d_obc_dir != -1 && param->d_obc_bulk == -1)
		{
		param->d_obc_bulk = param->d_size[param->d_obc_dir];
		}

	// parallel tempering
	REQUIRE(param->d_defect_dir >= 0 && param->d_defect_dir < STDIM,
	        "defect_dir must be in [0, %d)",
	        STDIM);
	for(i = 0; i < STDIM - 1; i++)
		{
		int j = (i < param->d_defect_dir) ? i : i + 1;
		REQUIRE(param->d_L_defect[i] <= param->d_size[j],
		        "defect length %d (%d) exceeds lattice size %d (%d)",
		        i, param->d_L_defect[i], j, param->d_size[j]);
		}

	// check on gradflow parameters
	REQUIRE(param->d_agf_meas_each > param->d_agf_time_bin,
	        "agf_meas_each must be greater than agf_time_bin");

	REQUIRE(param->d_agf_meas_each > param->d_agf_step,
	        "agf_meas_each must be greater than agf_step");

	REQUIRE(param->d_agf_meas_each == 0 ||
	        param->d_agf_meas_each > MIN_VALUE,
	        "if nonzero, agf_meas_each must be > MIN_VALUE in macro.h");

	// check on topological observables and theta term
	#if !((STDIM == 4 && NCOLOR > 1) || (STDIM == 2 && NCOLOR == 1))

	err = 0;
	err += param->d_charge_meas;
	err += param->d_charge_prime_meas;
	err += param->d_chi_prime_meas;
	err += param->d_charge_slices_meas;
	err += param->d_charge_p_slices_meas;
	REQUIRE(err == 0, "topological observables not allowed with STDIM=%d and NCOLOR=%d", STDIM, NCOLOR);

	#ifdef MULTICANONICAL_MODE
	REQUIRE(0, "multicanonical mode not supported with STDIM=%d and NCOLOR=%d", STDIM, NCOLOR);
	#endif

	#endif

	#ifdef THETA_MODE
	REQUIRE(STDIM == 4, "theta term can only be used in 4 dimensions");
	REQUIRE(param->d_theta_profile_dir >= -1 && param->d_theta_profile_dir < STDIM,
	        "direction of theta-term profile must be -1 (constant profile) or in [0, %d)",
	        STDIM);
	#endif

	for(i = 0; i < param->d_multipolyakov_order; i++)
		{
		REQUIRE(param->d_multipolyakov_dirs[i] >= 0 &&
		        param->d_multipolyakov_dirs[i] < STDIM,
		        "multipolyakov_dirs[%d] must be in [0, %d)",
		        i, STDIM);
		}

	// check on filenames
	check_required_string(param->d_conf_file, "conf_file", 1);
	check_required_string(param->d_twist_file, "twist_file", 1);
	check_required_string(param->d_log_file, "log_file", 1);
	check_required_string(param->d_data_file, "data_file", 1);
	check_required_string(param->d_energydensity_file, "energy_density_file", param->d_energy_density_meas);
	check_required_string(param->d_chargedensity_file, "charge_density_file", param->d_charge_density_meas);
	check_required_string(param->d_polyakovdensity_file, "polyakov_density_file", param->d_polyakov_density_meas);
	check_required_string(param->d_chiprime_file, "chiprime_data_file", param->d_chi_prime_meas);
	check_required_string(param->d_energy_slices_file, "energy_slices_file", param->d_energy_slices_meas);
	check_required_string(param->d_charge_slices_file, "charge_slices_file", param->d_charge_slices_meas);
	check_required_string(param->d_charge_slices_file, "charge_slices_file", param->d_charge_p_slices_meas);
	check_required_string(param->d_swap_acc_file, "swap_acc_file", param->d_N_replica_pt > 1);
	check_required_string(param->d_swap_tracking_file, "swap_track_file", param->d_N_replica_pt > 1);

	init_derived_constants(param);

	#ifdef THETA_MODE
	if(param->d_theta_profile_dir != -1)
		{
		check_required_string(param->d_theta_profile_file, "theta_profile_file", 1);
		read_theta_profile(param);
		}
	#endif

	#ifdef MULTICANONICAL_MODE
	check_required_string(param->d_multicanonic_acc_file, "multicanonic_acc_file", 1);
	check_required_string(param->d_topo_potential_file, "topo_potential_file", 1);
	read_topo_potential(param);
	#endif
	}


static inline Cooling_Type cooling_type_from_string(char const *str)
	{
	if(strcmp(str, "LEX_SITE_LEX_DIR") == 0)
		return LEX_SITE_LEX_DIR;
	if(strcmp(str, "LEX_DIR_LEX_SITE") == 0)
		return LEX_DIR_LEX_SITE;
	if(strcmp(str, "LEX_DIR_LEXEO_SITE") == 0)
		return LEX_DIR_LEXEO_SITE;
	if(strcmp(str, "LEXEO_SITE_LEX_DIR") == 0)
		return LEXEO_SITE_LEX_DIR;
	if(strcmp(str, "RND_DIR_RNDEO_SITE") == 0)
		return RND_DIR_RNDEO_SITE;
	REQUIRE(0, "unknown cooling type \"%s\"", str);
	}

static inline Multilevel_Obs ml_obs_from_string(char const *str)
	{
	if(strcmp(str, "NONE") == 0)
		return NONE;
	if(strcmp(str, "POLYCORR") == 0)
		return POLYCORR;
	if(strcmp(str, "POLYCORR_LONG") == 0)
		return POLYCORR_LONG;
	if(strcmp(str, "TUBE_DISC") == 0)
		return TUBE_DISC;
	if(strcmp(str, "TUBE_CONN") == 0)
		return TUBE_CONN;
	if(strcmp(str, "TUBE_CONN_LONG") == 0)
		return TUBE_CONN_LONG;
	REQUIRE(0, "unknown multilevel observable \"%s\"", str);
	}


// read normalized profile of theta term from file
void read_theta_profile(GParam *const param)
	{
	FILE *fp = fopen(param->d_theta_profile_file, "r");
	REQUIRE(fp != NULL, "failed to open theta-term profile file %s", param->d_theta_profile_file);

	allocate_array_double(&(param->d_theta_profile), param->d_theta_profile_size, __FILE__, __LINE__);

	// read theta profile from file
	for(int i = 0; i < param->d_theta_profile_size; i++)
		{
		int const err = fscanf(fp, "%lf", &(param->d_theta_profile[i]));
		REQUIRE(err == 1, "can't read theta profile at index %d", i);
		}
	double extra;
	int const err = fscanf(fp, "%lf", &extra);
	REQUIRE(err == EOF, "theta profile contains more than %d values", param->d_theta_profile_size);

	fclose(fp);
	}


// write normalized profile of theta term to file with name
void write_theta_profile(GParam const *const param, char const *const filename)
	{
	FILE *fp = fopen(filename, "w");
	REQUIRE(fp != NULL, "failed to open theta-term profile file %s", filename);

	// write theta profile to file
	for(int i = 0; i < param->d_theta_profile_size; i++)
		{
		fprintf(fp, "% 12.6e ", param->d_theta_profile[i]);
		}
	fprintf(fp, "\n");

	fclose(fp);
	}


// read topo potential from file
void read_topo_potential(GParam *const param)
	{
	FILE *fp = fopen(param->d_topo_potential_file, "r");
	REQUIRE(fp != NULL, "failed to open topological potential file %s", param->d_topo_potential_file);

	allocate_array_double_pointer(&(param->d_grid), param->d_N_replica_pt, __FILE__, __LINE__);
	for(int a = 0; a < param->d_N_replica_pt; a++)
		{
		allocate_array_double(&(param->d_grid[a]), param->d_n_grid, __FILE__, __LINE__);
		}

	// read x and V_a(x) from topo_potential file
	for(int i = 0; i < param->d_n_grid; i++)
		{
		double x, V;

		// read x
		int err = fscanf(fp, "%lf", &x);
		REQUIRE(err == 1, "can't read x at row %d", i);
		int j = (int) floor((x + param->d_grid_max + (param->d_grid_step / 2.0)) / param->d_grid_step);
		REQUIRE(i == j, "grid mismatch: found %d for x=%lf, expected %d", j, x, i);

		// read V_a(x)
		for(int a = 0; a < param->d_N_replica_pt; a++)
			{
			err = fscanf(fp, "%lf", &V);
			REQUIRE(err == 1, "can't read V_%d at row %d", a, i);
			param->d_grid[a][i] = V;
			}
		}

	fclose(fp);
	}


// write topo potential to file with name
void write_topo_potential(GParam const *const param, char const *const filename)
	{
	FILE *fp = fopen(filename, "w");
	REQUIRE(fp != NULL, "failed to open topological potential file %s", filename);

	// write x and V_a(x)
	for(int i = 0; i < param->d_n_grid; i++)
		{
		// write x
		double x = i * param->d_grid_step - param->d_grid_max;
		fprintf(fp, "% 12.6e ", x);

		// write V_a(x)
		for(int a = 0; a < param->d_N_replica_pt; a++)
			{
			fprintf(fp, "% 12.6e ", param->d_grid[a][i]);
			}
		fprintf(fp, "\n");
		}
	fprintf(fp, "\n");
	fclose(fp);
	}


void init_derived_constants(GParam *const param)
	{
	int i;

	// derived constants
	param->d_max_size = param->d_size[0];
	param->d_min_size = param->d_size[0];
	param->d_volume = 1;
	param->d_even_volume = 1;
	for(i = 0; i < STDIM; i++) param->d_space_vol[i] = 1;
	for(i = 0; i < STDIM; i++)
		{
		param->d_even_size[i] = param->d_size[i] - (param->d_size[i] % 2);
		if(param->d_size[i] > param->d_max_size) param->d_max_size = param->d_size[i];
		if(param->d_size[i] < param->d_min_size) param->d_min_size = param->d_size[i];
		(param->d_volume) *= (param->d_size[i]);
		(param->d_even_volume) *= (param->d_even_size[i]);
		for(int j = 0; j < STDIM; j++) if(j != i) (param->d_space_vol[j]) *= (param->d_size[i]);
		}

	param->d_n_even = param->d_even_volume / 2;
	param->d_n_border = param->d_volume - param->d_even_volume;
	param->d_inv_vol = 1.0 / ((double) param->d_volume);
	for(i = 0; i < STDIM; i++) param->d_inv_space_vol[i] = 1.0 / ((double) param->d_space_vol[i]);

	// volume of the defect
	param->d_volume_defect = 1;
	for(i = 0; i < STDIM - 1; i++)
		{
		param->d_volume_defect *= param->d_L_defect[i];
		}

	// number of planes (twisted boundary conditions only)
	param->d_n_planes = STDIM * (STDIM - 1);

	// default open boundary position
	param->d_obc_default_pos = 0;
	if(param->d_obc_dir != -1)
		param->d_obc_default_pos = param->d_size[param->d_obc_dir] - 1;

	// number of slices for theta-term profile
	param->d_theta_profile_size = param->d_size[param->d_theta_profile_dir];

	// number of grid points for multicanonical potential
	param->d_n_grid = (int) ((2.0 * param->d_grid_max / param->d_grid_step) + 1.0);

	// multilevel observable
	param->d_ml_obs = ml_obs_from_string(param->d_ml_obs_str);

	// number of multilevel hits (assuming Polyakov loops are separated along the "1" direction)
	if(param->d_dist_poly > 1 && param->d_size[1] - param->d_dist_poly > 1)
		param->d_ml_num_hit = param->d_multihit;
	else
		param->d_ml_num_hit = 0;

	// number of multilevel slices
	for(i = 0; i < NLEVELS; i++)
		param->d_ml_num_slices[i] = param->d_size[0] / param->d_ml_step[i];

	// number of measurements during gradient-flow evolution
	if(param->d_gf_meas_each > 0)
		param->d_gf_num_meas = (int) (param->d_ngfsteps / param->d_gf_meas_each);
	else
		param->d_gf_num_meas = 0;

	if(param->d_agf_meas_each > 0)
		param->d_agf_num_meas = (int) ((param->d_agf_length + MIN_VALUE) / param->d_agf_meas_each);
	else
		param->d_agf_num_meas = 0;

	// cooling type
	param->d_cooling_type = cooling_type_from_string(param->d_cooling_type_str);
	}


// edit parameters for measure-only mode
void edit_params_meas_only(GParam *const param)
	{
	// start from saved conf.
	param->d_start = 2;

	// not to overwrite files of runs with online measurements
	char suffix[STD_STRING_LENGTH];
	sprintf(suffix, "%s", "");
	if(param->d_agf_num_meas > 0)
		strcat(suffix, "_agf");
	if(param->d_gf_num_meas > 0)
		strcat(suffix, "_gf");
	if(param->d_coolrepeat > 0)
		{
		char aux[STD_STRING_LENGTH + 6];
		sprintf(aux, "_cool_%s", param->d_cooling_type_str);
		strcat(suffix, aux);
		}
	if(strcmp(suffix, "") == 0)
		strcat(suffix, "_hot");

	char *files[] = {
			param->d_data_file,
			param->d_energydensity_file,
			param->d_energydensity_file,
			param->d_polyakovdensity_file,
			param->d_chiprime_file,
			param->d_energy_slices_file,
			param->d_charge_slices_file,
			param->d_log_file
			};

	for(size_t i = 0; i < sizeof(files) / sizeof(files[0]); i++)
		strcat(files[i], suffix);
	}


// free allocated memory for hierarc update parameters
void free_hierarc_params(GParam *const param)
	{
	if(param->d_N_hierarc_levels == 0)
		{
		(void) param; // to avoid compiler warning about unused variable
		}
	else
		{
		free(param->d_L_rect);
		free(param->d_N_sweep_rect);
		}
	}

// print simulation parameters aux

void print_configuration_parameters(FILE *fp)
	{
	#ifdef OPENMP_MODE
	fprintf(fp, "Using OpenMP with %d threads\n\n", NTHREADS);
	#endif

	if(endian() == 0) fprintf(fp, "Little endian machine\n\n");
	else fprintf(fp, "Big endian machine\n\n");
	}


void print_pt_parameters(FILE *fp, GParam const *const param)
	{
	int i;
	fprintf(fp, "Using Parallel Tempering with\n");
	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
	fprintf(fp, "defect sizes: ");
	for(i = 0; i < STDIM - 1; i++) fprintf(fp, "%d ", param->d_L_defect[i]);
	fprintf(fp, "\n");
	fprintf(fp, "number of replicas: %d\n", param->d_N_replica_pt);
	fprintf(fp, "boundary conditions: ");
	for(i = 0; i < param->d_N_replica_pt; i++) fprintf(fp, "%lf ", param->d_pt_bound_cond_coeff[i]);
	fprintf(fp, "\n");
	fprintf(fp, "hierarchical levels: %d\n", param->d_N_hierarc_levels);
	if(param->d_N_hierarc_levels > 0)
		{
		fprintf(fp, "extentions of hierarchical rectangles: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) fprintf(fp, "%d ", param->d_L_rect[i]);
		fprintf(fp, "\n");
		fprintf(fp, "sweeps per hierarchical level: ");
		for(i = 0; i < param->d_N_hierarc_levels; i++) fprintf(fp, "%d ", param->d_N_sweep_rect[i]);
		fprintf(fp, "\n");
		}
	fprintf(fp, "\n");
	}


void print_multicanonic_parameters(FILE *fp, GParam const *const param)
	{
	fprintf(fp, "Using Multicanonical method with\n");
	fprintf(fp, "Multicanonic topo-potential read from file %s\nPotential defined on a grid with step=%.10lf and max=%.10lf\n", param->d_topo_potential_file, param->d_grid_step, param->d_grid_max);
	fprintf(fp, "topo_cooling:     %d\n", param->d_topo_cooling);
	fprintf(fp, "topo_coolsteps:   %d\n", param->d_topo_coolsteps);
	fprintf(fp, "topo_alpha:       %lf\n", param->d_topo_alpha);
	fprintf(fp, "\n");
	}


void print_multicanonic_tuning_parameters(FILE *fp, GParam const *const param)
	{
	fprintf(fp, "Tuning Multicanonical method with\n");
	fprintf(fp, "topo_tuning_thr:         %lf\n", param->d_topo_tuning_thr);
	fprintf(fp, "topo_tuning_stp:         %lf\n", param->d_topo_tuning_stp);
	fprintf(fp, "topo_tuning_save_every:  %d\n", param->d_topo_tuning_save_every);
	fprintf(fp, "topo_tuning_even:        %d\n", param->d_topo_tuning_even);
	fprintf(fp, "\n");
	}


void print_simul_parameters(FILE *fp, GParam const *const param)
	{
	int i;
	fprintf(fp, "colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimension: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(i = 1; i < STDIM; i++) fprintf(fp, "x%d", param->d_size[i]);
	fprintf(fp, "\n\n");

	if(param->d_obc_dir != -1)
		{
		fprintf(fp, "obc_dir:  %d\n", param->d_obc_dir);
		fprintf(fp, "obc_bulk: %d\n", param->d_obc_bulk);
		}
	fprintf(fp, "twist parameters: ");
	for(i = 0; i < STDIM * (STDIM - 1) / 2; i++) fprintf(fp, "%d ", param->d_k_twist[i]);
	fprintf(fp, "\n");
	fprintf(fp, "beta:  %.10lf\n", param->d_beta);
	#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
	#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "plaquette_meas:        %d\n", param->d_plaquette_meas);
	fprintf(fp, "clover_energy_meas:    %d\n", param->d_clover_energy_meas);
	fprintf(fp, "energy_density_meas:   %d\n", param->d_energy_density_meas);
	fprintf(fp, "charge_meas:           %d\n", param->d_charge_meas);
	fprintf(fp, "charge_density_meas:   %d\n", param->d_charge_density_meas);
	fprintf(fp, "polyakov_meas:         %d\n", param->d_polyakov_meas);
	fprintf(fp, "polyakov_powers_meas:  %d\n", param->d_polyakov_powers_meas);
	fprintf(fp, "polyakov_density_meas: %d\n", param->d_polyakov_density_meas);
	fprintf(fp, "chi_prime_meas:        %d\n", param->d_chi_prime_meas);
	fprintf(fp, "energy_slices_meas:    %d\n", param->d_energy_slices_meas);
	fprintf(fp, "charge_slices_meas:    %d\n", param->d_charge_slices_meas);
	fprintf(fp, "charge_p_slices_meas:  %d\n", param->d_charge_p_slices_meas);
	fprintf(fp, "\n");

	fprintf(fp, "multipolyakov_order:   %d    ", param->d_multipolyakov_order);
	for(i = 0; i < param->d_multipolyakov_order; i++) fprintf(fp, "%d ", param->d_multipolyakov_dirs[i]);
	fprintf(fp, "\n");
	fprintf(fp, "meas_effective_charge: %d\n", param->d_meas_effective_charge);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");
	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");
	}


void print_smoothing_parameters(FILE *fp, GParam const *const param)
	{
	if(param->d_agf_num_meas > 0)
		{
		fprintf(fp, "Using adaptive gradient flow with\n");
		fprintf(fp, "agf_length     %lf\n", param->d_agf_length);
		fprintf(fp, "agf_step:      %lf\n", param->d_agf_step);
		fprintf(fp, "agf_meas_each  %lf\n", param->d_agf_meas_each);
		fprintf(fp, "agf_delta      %e\n", param->d_agf_delta);
		fprintf(fp, "\n");
		}
	if(param->d_gf_num_meas > 0)
		{
		fprintf(fp, "Using fixed-step gradient flow with\n");
		fprintf(fp, "gfstep:        %lf\n", param->d_gfstep);
		fprintf(fp, "num_gfsteps    %d\n", param->d_ngfsteps);
		fprintf(fp, "gf_meas_each   %d\n", param->d_gf_meas_each);
		fprintf(fp, "\n");
		}
	if(param->d_coolrepeat > 0)
		{
		fprintf(fp, "Using cooling with\n");
		fprintf(fp, "cooling_type:  %s\n", param->d_cooling_type_str);
		fprintf(fp, "coolrepeat:    %d\n", param->d_coolrepeat);
		fprintf(fp, "coolsteps:     %d\n", param->d_coolsteps);
		fprintf(fp, "\n");
		}
	}


void print_multilevel_parameters(FILE *fp, GParam const *const param)
	{
	int i;
	fprintf(fp, "Using Multilevel algorithm for %s with\n", param->d_ml_obs_str);
	fprintf(fp, "multihit:  %d\n", param->d_multihit);
	fprintf(fp, "levels:    %d\n", NLEVELS);
	fprintf(fp, "steps:     ");
	for(i = 0; i < NLEVELS; i++)
		{
		fprintf(fp, "%d ", param->d_ml_step[i]);
		}
	fprintf(fp, "\n");
	fprintf(fp, "updates:   ");
	for(i = 0; i < NLEVELS; i++)
		{
		fprintf(fp, "%d ", param->d_ml_upd[i]);
		}
	fprintf(fp, "\n\n");
	}


void print_metro_parameters(FILE *fp, GParam const *const param, double acc)
	{
	fprintf(fp, "epsilon_metro:         %.10lf\n", param->d_epsilon_metro);
	fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
	fprintf(fp, "\n");
	}

// print simulation parameters

void print_parameters_local_agf(GParam const *const param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+---------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_agf |\n");
	fprintf(fp, "+---------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_local_pt_multicanonic(GParam const *const param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+---------------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_pt_multicanonic |\n");
	fprintf(fp, "+---------------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	print_multicanonic_parameters(fp, param);
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_local_pt_agf(GParam const *const param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_pt_agf |\n");
	fprintf(fp, "+------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_debug_agf_vs_gf(GParam const *const param, time_t time_start, time_t time_end, time_t agf_time, time_t dagf_time, time_t gf_time)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+----------------------------------------+\n");
	fprintf(fp, "| Simulation details for debug_agf_vs_gf |\n");
	fprintf(fp, "+----------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	double diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time:              %.3lf seconds\n", diff_sec);
	fprintf(fp, "Adaptive gradflow time:       %d seconds\n", (int) agf_time);
	fprintf(fp, "Debug adaptive gradflow time: %d seconds\n", (int) dagf_time);
	fprintf(fp, "Gradflow time:                %d seconds\n", (int) gf_time);
	fprintf(fp, "\n");

	fclose(fp);
	}


void print_parameters_debug_agf_vs_delta(GParam const *const param, time_t time_mc, time_t time_agf0, time_t time_agf1, time_t time_agf2, time_t time_agf3)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+----------------------------------------+\n");
	fprintf(fp, "| Simulation details for debug_agf_vs_gf |\n");
	fprintf(fp, "+----------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	fprintf(fp, "Simulation time:              %d seconds\n", (int) time_mc);
	fprintf(fp, "Adaptive gradflow time:\n");
	fprintf(fp, "    delta = %e:       %d seconds\n", param->d_agf_delta / 1.0000, (int) time_agf0);
	fprintf(fp, "    delta = %e:       %d seconds\n", param->d_agf_delta / 10.000, (int) time_agf1);
	fprintf(fp, "    delta = %e:       %d seconds\n", param->d_agf_delta / 100.00, (int) time_agf2);
	fprintf(fp, "    delta = %e:       %d seconds\n", param->d_agf_delta / 1000.0, (int) time_agf3);
	fprintf(fp, "\n");

	fclose(fp);
	}


void print_parameters_conf_analysis(GParam *param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+-------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_conf_analysis |\n");
	fprintf(fp, "+-------------------------------------------------+\n\n");

	print_configuration_parameters(fp);

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for(int i = 1; i < STDIM; i++)
		{
		fprintf(fp, "x%d", param->d_size[i]);
		}
	fprintf(fp, "\n\n");
	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	print_smoothing_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_tracedef(GParam const *const param, Time_Utils const *const timers, double acc)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+--------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tracedef |\n");
	fprintf(fp, "+--------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	fprintf(fp, "htracedef: ");
	for(int i = 0; i < NCOLOR / 2; i++)
		{
		fprintf(fp, "%lf ", param->d_h[i]);
		}
	fprintf(fp, "\n\n");
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_metro_parameters(fp, param, acc);
	print_smoothing_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_polycorr_long(GParam *param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+-------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_polycorr_long |\n");
	fprintf(fp, "+-------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_multilevel_parameters(fp, param);

	fprintf(fp, "level0_repeat: %d\n", param->d_ml_level0_repeat);
	fprintf(fp, "\n");
	fprintf(fp, "dist_poly:     %d\n", param->d_dist_poly);
	fprintf(fp, "\n");

	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_multilevel(GParam *param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+----------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_multilevel |\n");
	fprintf(fp, "+----------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_multilevel_parameters(fp, param);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_multilevel_long(GParam *param, Time_Utils const *const timers)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+---------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_multilevel_long |\n");
	fprintf(fp, "+---------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	#ifdef MULTICANONICAL_MODE
	print_multicanonic_parameters(fp, param);
	#endif
	print_simul_parameters(fp, param);
	print_multilevel_parameters(fp, param);
	fprintf(fp, "level0_repeat:    %d\n", param->d_ml_level0_repeat);
	fprintf(fp, "\n");

	print_time_utils(fp, timers);

	fclose(fp);
	}


void print_parameters_tuning_pt_mc(GParam const *const param, Time_Utils const *const timers, int count)
	{
	FILE *fp = fopen(param->d_log_file, "w");
	REQUIRE(fp != NULL, "failed to open log file %s", param->d_log_file);

	fprintf(fp, "+------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_tuning_pt_mc |\n");
	fprintf(fp, "+------------------------------------------------+\n\n");

	print_configuration_parameters(fp);
	print_pt_parameters(fp, param);
	print_multicanonic_parameters(fp, param);
	print_multicanonic_tuning_parameters(fp, param);
	print_simul_parameters(fp, param);
	print_smoothing_parameters(fp, param);

	fprintf(fp, "Tuning steps: %d\n", count);
	print_time_utils(fp, timers);

	fclose(fp);
	}

// print template input aux

void print_template_volume_parameters(FILE *fp)
	{
	fprintf(fp, "size 12 4 4 12  # lattice sizes along space-time directions 0 (time), 1, ...\n");
	fprintf(fp, "\n");
	}


void print_template_simul_parameters(FILE *fp)
	{
	fprintf(fp, "# Simulations parameters\n");
	fprintf(fp, "beta   6.4881\n");
	#ifdef THETA_MODE
	fprintf(fp, "theta  0.5\n");
	#endif
	fprintf(fp, "\n");
	fprintf(fp, "sample     10 # total number of Monte Carlo steps\n");
	fprintf(fp, "thermal     0 # number of thermalization steps before measurements\n");
	fprintf(fp, "overrelax   5 # number of overrelaxation sweeps for heatbath sweep (the set of all these sweeps counts as 1 MC step)\n");
	fprintf(fp, "measevery   1 # number of MC steps between measurements\n");
	fprintf(fp, "\n");

	fprintf(fp, "start                    3  # 0=all links to identity  1=random  2=from saved configuration 3=ordered with twisted bc\n");
	fprintf(fp, "saveconf_back_every      5  # if 0, do not save; otherwise, save backup configurations every this many MC steps\n");
	fprintf(fp, "saveconf_analysis_every  5  # if 0, do not save; otherwise, save configurations for analysis every this many MC steps\n");
	fprintf(fp, "\n");

	fprintf(fp, "randseed 0    # RNG seed (0=time)\n");
	fprintf(fp, "walltime 24   # execution wall time in hours\n");
	fprintf(fp, "\n");

	fprintf(fp, "# Observables to measure\n");
	fprintf(fp, "plaquette_meas        0  # 1=YES, 0=NO\n");
	fprintf(fp, "clover_energy_meas    1  # 1=YES, 0=NO\n");
	fprintf(fp, "energy_density_meas   0  # 1=YES, 0=NO\n");
	fprintf(fp, "charge_meas           1  # 1=YES, 0=NO\n");
	fprintf(fp, "charge_density_meas   0  # 1=YES, 0=NO\n");
	fprintf(fp, "polyakov_meas         0  # 1=YES, 0=NO\n");
	fprintf(fp, "polyakov_density_meas 0  # 1=YES, 0=NO\n");
	fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
	fprintf(fp, "energy_slices_meas    0  # 1=YES, 0=NO\n");
	fprintf(fp, "charge_slices_meas    0  # 1=YES, 0=NO\n");
	fprintf(fp, "\n");
	fprintf(fp, "multipolyakov_order   0  # n  mu_1 mu_2 ... mu_n\n");
	fprintf(fp, "meas_effective_charge 0  # if 1, eventual x-dependence of theta-term moved in definition of Q\n");
	fprintf(fp, "\n");
	}


void print_template_pt_parameters(FILE *fp)
	{
	fprintf(fp, "# Parallel tempering parameters\n");
	fprintf(fp, "defect_dir    0             # space-time direction orthogonal to the topological defect\n");
	fprintf(fp, "defect_size   2 2 2         # sizes of the defect along its dimensions (skipping defect_dir)\n");
	fprintf(fp, "N_replica_pt  2    1.0 0.0  # number of parallel tempering replicas, followed by corresponding boundary conditions coefficients\n");
	fprintf(fp, "\n");
	fprintf(fp, "# Hierarchical update parameters\n");
	fprintf(fp, "# number of hierarchical levels, extensions of rectangles, sweeps per rectangle\n");
	fprintf(fp, "hierarc_upd 2    2 1    1 1\n");
	fprintf(fp, "\n");
	}


void print_template_twist_parameters(FILE *fp)
	{
	fprintf(fp, "# Twist parameters\n");
	fprintf(fp, "k_twist 0 0 0 1 0 0    # twist parameter on the plane (0,1), (0,2), ..., (0,STDIM-1), (1, 2), ...\n");
	fprintf(fp, "\n");
	}


void print_template_adaptive_gradflow_parameters(FILE *fp)
	{
	fprintf(fp, "# For gradient flow evolution with adaptive-step integrator\n");
	fprintf(fp, "agf_length       10    # total integration time\n");
	fprintf(fp, "agf_step       0.01    # initial integration step\n");
	fprintf(fp, "agf_meas_each     1    # time interval between measurements\n");
	fprintf(fp, "agf_delta     0.001    # error threshold on gauge links for adapting the integration step\n");
	fprintf(fp, "agf_time_bin      0    # error threshold on time at which measure are taken\n");
	fprintf(fp, "\n");
	}


void print_template_gradflow_parameters(FILE *fp)
	{
	fprintf(fp, "# For gradient flow evolution with fixed-step integrator\n");
	fprintf(fp, "gfstep      0.01    # integration step\n");
	fprintf(fp, "num_gfsteps  100    # total number of integration steps\n");
	fprintf(fp, "gf_meas_each   5    # number of integration steps between measurements\n");
	fprintf(fp, "\n");
	}


void print_template_cooling_parameters(FILE *fp)
	{
	fprintf(fp, "# For cooling\n");
	fprintf(fp, "coolsteps             3  # number of cooling sweeps between measurements\n");
	fprintf(fp, "coolrepeat            5  # number of measurements during cooling\n");
	fprintf(fp, "\n");
	}


void print_template_metro_parameters(FILE *fp)
	{
	fprintf(fp, "epsilon_metro    0.25  # distance from the identity of the random matrix for metropolis\n");
	fprintf(fp, "\n");
	}


void print_template_multicanonic_parameters(FILE *fp)
	{
	fprintf(fp, "# Multicanonic parameters\n");
	fprintf(fp, "grid_step                0.05                  # charge steps at which topo_potential is defined in topo_potential_file\n");
	fprintf(fp, "grid_max                 3.0                   # abs value of charge at which topo_potential saturates in topo_potential_file\n");
	fprintf(fp, "topo_cooling             0                     # cooling strat before evaluating the topo potential: 0 = none, 1 = cooling\n");
	fprintf(fp, "topo_coolsteps           5                     # cooling steps before evaluating the topo potential (if topo_cooling = 1)\n");
	fprintf(fp, "topo_alpha               1.0                   # used for alpha-rounding if >0, no alpha-rounding if =0\n");
	fprintf(fp, "topo_potential_file      topo_potential        # file to read the topo_potential from\n");
	fprintf(fp, "multicanonic_acc_file    multicanonic_acc.dat  # file to save acceptances of Metropolis tests with topo_potential\n");
	fprintf(fp, "\n");
	}


void print_template_multicanonic_tuning_parameters(FILE *fp)
	{
	fprintf(fp, "# Multicanonic tuning parameters\n");
	fprintf(fp, "topo_tuning_thr           0.05 # tuning ends if topo potential changes less than this threshold\n");
	fprintf(fp, "topo_tuning_stp           0.1  # maximum variation of topo_potential at each point every step \n");
	fprintf(fp, "topo_tuning_save_every    1    # save topo_potential at each point every step \n");
	fprintf(fp, "topo_tuning_even          1    # force topo_potential to be even during tuning (0 = False, 1 = True) \n");
	fprintf(fp, "\n");
	}


void print_template_multilevel_parameters(FILE *fp)
	{
	fprintf(fp, "# Multilevel parameters\n");
	fprintf(fp, "multihit         10  # number of multihit steps\n");
	fprintf(fp, "ml_step           2  # timeslices for multilevel (from largest to smallest)\n");
	fprintf(fp, "ml_upd           10  # number of updates for various levels\n");
	fprintf(fp, "ml_file      ml.dat  # multilevel output file\n");
	fprintf(fp, "ml_obs     POLYCORR  # observable to be measured during multilevel (see Multilevel_Obs in gparam.h)\n");
	fprintf(fp, "\n");
	fprintf(fp, "dist_poly      2 # distance between the polyakov loop\n");
	fprintf(fp, "transv_dist    2 # transverse distance from the polyakov correlator\n");
	fprintf(fp, "plaq_dir     1 0 # plaquette orientation for flux tube\n");
	fprintf(fp, "\n");
	}


void print_template_output_parameters(FILE *fp)
	{
	fprintf(fp, "# Output files\n");
	fprintf(fp, "conf_file             conf.dat\n");
	fprintf(fp, "twist_file            twist.dat\n");
	fprintf(fp, "data_file             dati.dat\n");
	fprintf(fp, "energy_density_file   energy_density.dat\n");
	fprintf(fp, "charge_density_file   charge_density.dat\n");
	fprintf(fp, "polyakov_density_file polyakov_density.dat\n");
	fprintf(fp, "chiprime_data_file    chi_prime_cool.dat\n");
	fprintf(fp, "energy_slices_file    energy_slices.dat\n");
	fprintf(fp, "charge_slices_file    charge_slices.dat\n");
	fprintf(fp, "log_file              log.dat\n");
	fprintf(fp, "swap_acc_file         swap_acc.dat\n");
	fprintf(fp, "swap_track_file       swap_track.dat\n");
	fprintf(fp, "\n");
	}

// print program details

void print_authors(int parallel_tempering, int twisted_bc)
	{
	printf("\n");
	if(parallel_tempering == 1)
		printf("SU(N) Hasenbusch Parallel Tempering implemented by Claudio Bonanno (claudiobonanno93@gmail.com)\n");
	if(twisted_bc == 1)
		printf("Twisted Boundary Conditions implemented by Andrea Giorgieri (andrea.giorgieri.pi@gmail.com)\n");
	if(parallel_tempering == 1 || twisted_bc == 1)
		printf("within yang-mills package\n\n");

	printf("Details about yang-mills package:\n");
	printf("\tPackage %s version: %s\n", PACKAGE_NAME, PACKAGE_VERSION);
	printf("\tAuthor: Claudio Bonati %s\n\n", PACKAGE_BUGREPORT);
	}


void print_compilation_details(void)
	{
	printf("Compilation details:\n");
	printf("\tN_c (number of colors): %d\n", NCOLOR);
	printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
	printf("\tNum_levels (number of levels): %d\n", NLEVELS);
	printf("\n");
	printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
	printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

	#ifdef DEBUG
	printf("\n\tDEBUG mode\n");
	#endif

	#ifdef OPENMP_MODE
	printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
	#endif

	#ifdef THETA_MODE
	printf("\n\tusing imaginary theta\n");
	#endif

	#ifdef MULTICANONICAL_MODE
	printf("\n\tusing multicanonical algorithm\n");
	#endif

	#ifdef OPT_MULTIHIT
	printf("\tcompiled for multihit optimization\n");
	#endif

	#ifdef OPT_MULTILEVEL
	printf("\tcompiled for multilevel optimization\n");
	#endif

	printf("\n");

	#ifdef __INTEL_COMPILER
	printf("\tcompiled with icc\n");
	#elif defined(__clang__)
	printf("\tcompiled with clang\n");
	#elif defined( __GNUC__ )
	printf("\tcompiled with gcc version: %d.%d.%d\n",
	       __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	#endif
	printf("\n");
	}


#endif
