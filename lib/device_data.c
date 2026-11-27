#ifndef DEVICE_DATA_C
#define DEVICE_DATA_C

#include "../include/macro.h"
#include "../include/device_data.h"

#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef OPENACC_MODE
#include <openacc.h>
#endif

#include "../include/device.h"
#include "../include/gauge_conf.h"
#include "../include/gauge_group.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/memalign.h"
#include "../include/random.h"


#ifdef OPENACC_MODE

// ---------------------------------------------------------------------------------------------------------------------
// generic helpers, acting on the current device; "member" is the address of a pointer member of a structure that is
// already present on the device: after the copy the member is attached (its device copy points to device memory)
// ---------------------------------------------------------------------------------------------------------------------

// copy a table of n host pointers into the block h_target to the device, translated into pointers into the device copy
// of the block (which must be already present); the table must be copied only once per device
static void table_to_device(void **h_table, long const n, void *h_target)
	{
	if(g_dev.on_host) return; // shared memory (-acc=multicore): the host table is already valid, nothing to translate
	void **d_table = (void **) acc_copyin(h_table, (size_t) n * sizeof(void *));
	uintptr_t const h_base = (uintptr_t) h_target;
	uintptr_t const d_base = (uintptr_t) acc_deviceptr(h_target);
	#pragma acc parallel loop gang vector deviceptr(d_table)
	for(long i = 0; i < n; i++)
		d_table[i] = (void *) ((uintptr_t) d_table[i] - h_base + d_base);
	}


// array [bytes]: copied (copy = 1) or only allocated (copy = 0)
static void enter_1d(void **member, size_t const bytes, int const copy)
	{
	if(copy) acc_copyin(*member, bytes);
	else acc_create(*member, bytes);
	acc_attach(member);
	}

static void exit_1d(void **member, size_t const bytes)
	{
	acc_detach(member);
	acc_delete(*member, bytes);
	}


// array [n1][...] with block (*member)[0] of block_bytes bytes (allocate_array_*_2d)
static void enter_2d(void ***member, long const n1, size_t const block_bytes, int const copy)
	{
	void **table = *member;
	if(copy) acc_copyin(table[0], block_bytes);
	else acc_create(table[0], block_bytes);
	table_to_device(table, n1, table[0]);
	acc_attach((void **) member);
	}

static void exit_2d(void ***member, long const n1, size_t const block_bytes)
	{
	void **table = *member;
	acc_detach((void **) member);
	acc_delete(table, (size_t) n1 * sizeof(void *));
	acc_delete(table[0], block_bytes);
	}


// array [n_outer][...][...] with n_rows rows (*member)[0][k] and block (*member)[0][0] of block_bytes bytes
// (allocate_array_*_3d, geo->d_mutimespace)
static void enter_3d(void ****member, long const n_outer, long const n_rows, size_t const block_bytes, int const copy)
	{
	void ***outer = *member;
	void **rows = outer[0];
	if(copy) acc_copyin(rows[0], block_bytes);
	else acc_create(rows[0], block_bytes);
	table_to_device(rows, n_rows, rows[0]);
	table_to_device((void **) outer, n_outer, rows);
	acc_attach((void **) member);
	}

static void exit_3d(void ****member, long const n_outer, long const n_rows, size_t const block_bytes)
	{
	void ***outer = *member;
	void **rows = outer[0];
	void *block = rows[0];
	acc_detach((void **) member);
	acc_delete(outer, (size_t) n_outer * sizeof(void *));
	acc_delete(rows, (size_t) n_rows * sizeof(void *));
	acc_delete(block, block_bytes);
	}


// array of n rectangles (with their sites)
static void enter_rect_array(Rectangle **member, long const n)
	{
	acc_copyin(*member, (size_t) n * sizeof(Rectangle));
	acc_attach((void **) member);
	for(long i = 0; i < n; i++)
		enter_1d((void **) &((*member)[i].rect_sites), (size_t) (*member)[i].d_vol_rect * sizeof(long), 1);
	}

static void exit_rect_array(Rectangle **member, long const n)
	{
	for(long i = 0; i < n; i++)
		exit_1d((void **) &((*member)[i].rect_sites), (size_t) (*member)[i].d_vol_rect * sizeof(long));
	acc_detach((void **) member);
	acc_delete(*member, (size_t) n * sizeof(Rectangle));
	}


// update of the block of a 2d array
static void update_2d(void **table, size_t const bytes, int const to_device)
	{
	if(to_device) acc_update_device(table[0], bytes);
	else acc_update_self(table[0], bytes);
	}


// ---------------------------------------------------------------------------------------------------------------------
// sizes
// ---------------------------------------------------------------------------------------------------------------------

static size_t lattice_bytes(GParam const *const param)
	{
	return (size_t) param->d_volume * STDIM * sizeof(GAUGE_GROUP);
	}

static size_t Z_bytes(GParam const *const param)
	{
	return (size_t) param->d_volume * (size_t) (param->d_n_planes + 1) * sizeof(double complex);
	}

static long mutimespace_rows(GParam const *const param)
	{
	long n_rows = 0;
	for(int i = 0; i < STDIM; i++)
		n_rows += param->d_size[i];
	return n_rows;
	}


// ---------------------------------------------------------------------------------------------------------------------
// shared data (current device)
// ---------------------------------------------------------------------------------------------------------------------

static void enter_shared(Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int const n_meas)
	{
	long const V = param->d_volume;
	int const Nr = param->d_N_replica_pt;

	// parameters
	acc_copyin(param, sizeof(GParam));
	enter_1d((void **) &(param->d_pt_bound_cond_coeff), (size_t) Nr * sizeof(double), 1);
	#ifdef MULTICANONICAL_MODE
	enter_2d((void ***) &(param->d_grid), Nr, (size_t) Nr * (size_t) param->d_n_grid * sizeof(double), 1);
	#endif

	// geometry
	acc_copyin(geo, sizeof(Geometry));
	enter_1d((void **) &(geo->d_sip_to_si), (size_t) V * sizeof(long), 1);
	enter_2d((void ***) &(geo->d_nnp), V, (size_t) V * STDIM * sizeof(long), 1);
	enter_2d((void ***) &(geo->d_nnm), V, (size_t) V * STDIM * sizeof(long), 1);
	enter_2d((void ***) &(geo->d_mucomp), STDIM, (size_t) STDIM * (size_t) V * sizeof(int), 1);
	enter_2d((void ***) &(geo->d_muorth), STDIM, (size_t) STDIM * (size_t) V * sizeof(long), 1);
	enter_3d((void ****) &(geo->d_mutimespace), STDIM, mutimespace_rows(param), (size_t) STDIM * (size_t) V * sizeof(long), 1);
	enter_1d((void **) &(geo->d_defect_mask), (size_t) V * sizeof(int), 1);

	// global lookup tables
	#if STDIM == 4
	#pragma acc update device(g_signed_ord_orth_dir, g_indep_perm_dir)
	#endif

	// rectangles (NULL for the programs without parallel tempering)
	int const L = param->d_N_hierarc_levels;
	if(rect_aux != NULL)
		{
		acc_copyin(rect_aux, sizeof(Rect_Utils));
		enter_1d((void **) &(rect_aux->swap_rect.rect_sites), (size_t) rect_aux->swap_rect.d_vol_rect * sizeof(long), 1);
		}
	if(rect_aux != NULL && L > 0)
		{
		enter_rect_array(&(rect_aux->update_rect), L);
		#ifdef THETA_MODE
		enter_rect_array(&(rect_aux->clover_rect), L);
		#endif
		#ifdef MULTICANONICAL_MODE
		enter_rect_array(&(rect_aux->topcharge_rect), L);
		if(rect_aux->cooling_rect != NULL)
			{
			acc_copyin(rect_aux->cooling_rect, (size_t) L * sizeof(Rectangle *));
			acc_attach((void **) &(rect_aux->cooling_rect));
			for(int i = 0; i < L; i++)
				enter_rect_array(&(rect_aux->cooling_rect[i]), param->d_topo_coolsteps);
			}
		#endif
		}

	// auxiliary lattices of the measurements: allocated by meas_aux[0], shared by the other elements
	acc_copyin(meas_aux, (size_t) n_meas * sizeof(Meas_Utils));
	if(meas_aux[0].lattice_aux[0] != NULL)
		{
		for(int k = 0; k < 4; k++)
			{
			enter_2d((void ***) &(meas_aux[0].lattice_aux[k]), V, lattice_bytes(param), 0);
			for(int m = 1; m < n_meas; m++)
				acc_attach((void **) &(meas_aux[m].lattice_aux[k]));
			}
		}
	}


static void exit_shared(Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int const n_meas)
	{
	long const V = param->d_volume;
	int const Nr = param->d_N_replica_pt;

	if(meas_aux[0].lattice_aux[0] != NULL)
		{
		for(int k = 0; k < 4; k++)
			{
			for(int m = 1; m < n_meas; m++)
				acc_detach((void **) &(meas_aux[m].lattice_aux[k]));
			exit_2d((void ***) &(meas_aux[0].lattice_aux[k]), V, lattice_bytes(param));
			}
		}
	acc_delete(meas_aux, (size_t) n_meas * sizeof(Meas_Utils));

	int const L = param->d_N_hierarc_levels;
	if(rect_aux != NULL && L > 0)
		{
		#ifdef MULTICANONICAL_MODE
		if(rect_aux->cooling_rect != NULL)
			{
			for(int i = 0; i < L; i++)
				exit_rect_array(&(rect_aux->cooling_rect[i]), param->d_topo_coolsteps);
			acc_detach((void **) &(rect_aux->cooling_rect));
			acc_delete(rect_aux->cooling_rect, (size_t) L * sizeof(Rectangle *));
			}
		exit_rect_array(&(rect_aux->topcharge_rect), L);
		#endif
		#ifdef THETA_MODE
		exit_rect_array(&(rect_aux->clover_rect), L);
		#endif
		exit_rect_array(&(rect_aux->update_rect), L);
		}
	if(rect_aux != NULL)
		{
		exit_1d((void **) &(rect_aux->swap_rect.rect_sites), (size_t) rect_aux->swap_rect.d_vol_rect * sizeof(long));
		acc_delete(rect_aux, sizeof(Rect_Utils));
		}

	exit_1d((void **) &(geo->d_defect_mask), (size_t) V * sizeof(int));
	exit_3d((void ****) &(geo->d_mutimespace), STDIM, mutimespace_rows(param), (size_t) STDIM * (size_t) V * sizeof(long));
	exit_2d((void ***) &(geo->d_muorth), STDIM, (size_t) STDIM * (size_t) V * sizeof(long));
	exit_2d((void ***) &(geo->d_mucomp), STDIM, (size_t) STDIM * (size_t) V * sizeof(int));
	exit_2d((void ***) &(geo->d_nnm), V, (size_t) V * STDIM * sizeof(long));
	exit_2d((void ***) &(geo->d_nnp), V, (size_t) V * STDIM * sizeof(long));
	exit_1d((void **) &(geo->d_sip_to_si), (size_t) V * sizeof(long));
	acc_delete(geo, sizeof(Geometry));

	#ifdef MULTICANONICAL_MODE
	exit_2d((void ***) &(param->d_grid), Nr, (size_t) Nr * (size_t) param->d_n_grid * sizeof(double));
	#endif
	exit_1d((void **) &(param->d_pt_bound_cond_coeff), (size_t) Nr * sizeof(double));
	acc_delete(param, sizeof(GParam));
	}


// ---------------------------------------------------------------------------------------------------------------------
// configurations (current device = device of the slots)
// ---------------------------------------------------------------------------------------------------------------------

static void enter_gauge_conf_members(Gauge_Conf *GC, GParam const *const param)
	{
	long const V = param->d_volume;
	enter_2d((void ***) &(GC->lattice), V, lattice_bytes(param), 1);
	enter_2d((void ***) &(GC->lattice_copy), V, lattice_bytes(param), 1);
	enter_2d((void ***) &(GC->Z), V, Z_bytes(param), 1);
	enter_2d((void ***) &(GC->Z_copy), V, Z_bytes(param), 1);
	#ifdef MULTICANONICAL_MODE
	enter_2d((void ***) &(GC->lattice_cold), V, lattice_bytes(param), 1);
	enter_2d((void ***) &(GC->lattice_copy_cold), V, lattice_bytes(param), 1);
	#endif
	#ifdef THETA_MODE
	// clovers are recomputed from the lattice before every use: only allocated
	enter_3d((void ****) &(GC->clover_array), V, V * STDIM, (size_t) V * STDIM * STDIM * sizeof(GAUGE_GROUP), 0);
	#endif
	#ifdef RNG_PHILOX_MODE
	if(GC->rng_site != NULL) enter_1d((void **) &(GC->rng_site), (size_t) V * sizeof(rng_t), 1); // NULL: never updated
	#endif
	}


static void exit_gauge_conf_members(Gauge_Conf *GC, GParam const *const param)
	{
	long const V = param->d_volume;
	#ifdef RNG_PHILOX_MODE
	if(GC->rng_site != NULL) exit_1d((void **) &(GC->rng_site), (size_t) V * sizeof(rng_t));
	#endif
	#ifdef THETA_MODE
	exit_3d((void ****) &(GC->clover_array), V, V * STDIM, (size_t) V * STDIM * STDIM * sizeof(GAUGE_GROUP));
	#endif
	#ifdef MULTICANONICAL_MODE
	exit_2d((void ***) &(GC->lattice_copy_cold), V, lattice_bytes(param));
	exit_2d((void ***) &(GC->lattice_cold), V, lattice_bytes(param));
	#endif
	exit_2d((void ***) &(GC->Z_copy), V, Z_bytes(param));
	exit_2d((void ***) &(GC->Z), V, Z_bytes(param));
	exit_2d((void ***) &(GC->lattice_copy), V, lattice_bytes(param));
	exit_2d((void ***) &(GC->lattice), V, lattice_bytes(param));
	}


// update of the fields of GC on the current device
static void update_gauge_conf(Gauge_Conf *GC, GParam const *const param, int const fields, int const to_device)
	{
	if(fields & GC_LATTICE) update_2d((void **) GC->lattice, lattice_bytes(param), to_device);
	if(fields & GC_COPY)
		{
		update_2d((void **) GC->lattice_copy, lattice_bytes(param), to_device);
		update_2d((void **) GC->Z_copy, Z_bytes(param), to_device);
		}
	if(fields & GC_Z) update_2d((void **) GC->Z, Z_bytes(param), to_device);
	#ifdef MULTICANONICAL_MODE
	if(fields & GC_COLD)
		{
		update_2d((void **) GC->lattice_cold, lattice_bytes(param), to_device);
		update_2d((void **) GC->lattice_copy_cold, lattice_bytes(param), to_device);
		}
	#endif
	#ifdef RNG_PHILOX_MODE
	if((fields & GC_RNG) && GC->rng_site != NULL)
		{
		size_t const bytes = (size_t) param->d_volume * sizeof(rng_t);
		if(to_device) acc_update_device(GC->rng_site, bytes);
		else acc_update_self(GC->rng_site, bytes);
		}
	#endif
	}

#endif // OPENACC_MODE


// ---------------------------------------------------------------------------------------------------------------------
// views (every build)
// ---------------------------------------------------------------------------------------------------------------------

static Gauge_Conf const *s_gc_base = NULL; // array of configurations given to init_device_data
static long s_n_slots = 0;                 // its number of elements
static Gauge_Conf *s_view[NUM_VIEWS];      // [kind][slot]
static Meas_Utils *s_meas = NULL;          // meas_aux given to init_device_data (auxiliary lattices on the devices)
static int s_force_host = 0;               // force_host_compute


static void init_views(Gauge_Conf const *GC, GParam const *const param)
	{
	(void) param;
	s_gc_base = GC;
	s_n_slots = device_num_slots();
	for(int k = 0; k < NUM_VIEWS; k++)
		{
		s_view[k] = (Gauge_Conf *) calloc((size_t) s_n_slots, sizeof(Gauge_Conf));
		REQUIRE(s_view[k] != NULL, "allocation of the views failed");
		for(long g = 0; g < s_n_slots; g++)
			s_view[k][g].device = GC[g].device;
		}

	#ifdef OPENACC_MODE
	// shallow copies: the pointers are assigned on the device by gauge_conf_view
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		long const first = device_slot_begin(d);
		long const n = device_slot_end(d) - first;
		for(int k = 0; k < NUM_VIEWS; k++)
			acc_copyin(&(s_view[k][first]), (size_t) n * sizeof(Gauge_Conf));
		}
	#endif
	}


static void free_views(void)
	{
	#ifdef OPENACC_MODE
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		long const first = device_slot_begin(d);
		long const n = device_slot_end(d) - first;
		for(int k = 0; k < NUM_VIEWS; k++)
			acc_delete(&(s_view[k][first]), (size_t) n * sizeof(Gauge_Conf));
		}
	#endif

	for(int k = 0; k < NUM_VIEWS; k++)
		{
		free(s_view[k]);
		s_view[k] = NULL;
		}
	s_gc_base = NULL;
	s_n_slots = 0;
	}


Gauge_Conf *gauge_conf_view(Gauge_Conf const *GC, View_Kind const kind)
	{
	long const g = (long) (GC - s_gc_base);
	REQUIRE(s_gc_base != NULL && g >= 0 && g < s_n_slots, "configuration without views (init_device_data not called?)");
	Gauge_Conf *v = &(s_view[kind][g]);

	switch(kind)
		{
		case VIEW_COLD:
			v->lattice = GC->lattice_cold;
			v->Z = GC->Z;
			v->lattice_copy = GC->lattice_copy_cold;
			v->Z_copy = GC->Z_copy;
			break;
		case VIEW_COPY:
			v->lattice = GC->lattice_copy;
			v->Z = GC->Z_copy;
			break;
		case VIEW_COLD_COPY:
			v->lattice = GC->lattice_copy_cold;
			v->Z = GC->Z_copy;
			break;
		case VIEW_FLOW:
			REQUIRE(s_meas != NULL && s_meas->lattice_aux[0] != NULL, "no auxiliary lattices of the measurements for the gradient flow");
			v->lattice = s_meas->lattice_aux[0];
			v->Z = GC->Z;
			break;
		default:
			REQUIRE(0, "unknown view kind %d", (int) kind);
		}

	#ifdef OPENACC_MODE
	// the same assignment on the device, with the device pointers of GC (and of the auxiliary lattices of the device)
	select_device(GC->device);
	int const k = (int) kind;
	Meas_Utils const *const ma = s_meas;
	if(k == VIEW_FLOW)
		{
		#pragma acc serial present(GC[0:1], v[0:1], ma[0:1]) async(QUEUE_DEVICE)
			{
			v->lattice = ma->lattice_aux[0];
			v->Z = GC->Z;
			}
		return v;
		}
	#pragma acc serial present(GC[0:1], v[0:1]) async(QUEUE_DEVICE)
		{
		if(k == VIEW_COLD)
			{
			v->lattice = GC->lattice_cold;
			v->Z = GC->Z;
			v->lattice_copy = GC->lattice_copy_cold;
			v->Z_copy = GC->Z_copy;
			}
		else if(k == VIEW_COPY)
			{
			v->lattice = GC->lattice_copy;
			v->Z = GC->Z_copy;
			}
		else
			{
			v->lattice = GC->lattice_copy_cold;
			v->Z = GC->Z_copy;
			}
		}
	#endif

	return v;
	}


// ---------------------------------------------------------------------------------------------------------------------
// public functions
// ---------------------------------------------------------------------------------------------------------------------

void select_device(int const d)
	{
	#ifdef OPENACC_MODE
	if(g_dev.on_host)
		acc_set_device_num(0, acc_device_host); // for simulated splitting among devices
	else
		acc_set_device_num(g_dev.id[d], acc_device_nvidia);
	#else
	(void) d;
	#endif
	}


void init_device_data(Gauge_Conf *GC, Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int const n_meas)
	{
	#ifdef OPENACC_MODE
	REQUIRE(g_dev.n_used > 0, "init_device_utils must be called before init_device_data");
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		enter_shared(geo, param, rect_aux, meas_aux, n_meas);

		long const first = device_slot_begin(d);
		long const n = device_slot_end(d) - first;
		acc_copyin(&(GC[first]), (size_t) n * sizeof(Gauge_Conf));
		for(long g = first; g < first + n; g++)
			enter_gauge_conf_members(&(GC[g]), param);
		}
	acc_wait_all();
	#else
	(void) geo;
	(void) rect_aux;
	(void) n_meas;
	#endif

	s_meas = meas_aux;
	init_views(GC, param);
	}


int gc_on_device(Gauge_Conf const *GC)
	{
	#ifdef OPENACC_MODE
	if(s_force_host || s_gc_base == NULL) return 0;
	if(GC->device < 0 || GC->device >= g_dev.n_used) return 0; // configuration outside the slots (device not set)
	select_device(GC->device);
	return acc_is_present((void *) GC, sizeof(Gauge_Conf));
	#else
	(void) GC;
	return 0;
	#endif
	}


void wait_device_queue(int const dev)
	{
	#ifdef OPENACC_MODE
	if(dev) acc_wait(QUEUE_DEVICE);
	#else
	(void) dev;
	#endif
	}


int force_host_compute(int const on)
	{
	int const previous = s_force_host;
	s_force_host = on;
	return previous;
	}


Meas_Utils *device_meas_aux(void)
	{
	return s_meas;
	}


void wait_devices(void)
	{
	#ifdef OPENACC_MODE
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		acc_wait_all();
		}
	#endif
	}


void free_device_data(Gauge_Conf *GC, Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int const n_meas)
	{
	wait_devices();
	free_views();
	s_meas = NULL;

	#ifdef OPENACC_MODE
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		acc_wait_all();

		long const first = device_slot_begin(d);
		long const n = device_slot_end(d) - first;
		for(long g = first; g < first + n; g++)
			exit_gauge_conf_members(&(GC[g]), param);
		acc_delete(&(GC[first]), (size_t) n * sizeof(Gauge_Conf));

		exit_shared(geo, param, rect_aux, meas_aux, n_meas);
		}
	#else
	(void) GC;
	(void) geo;
	(void) param;
	(void) rect_aux;
	(void) meas_aux;
	(void) n_meas;
	#endif
	}


void copy_gauge_conf_to_host(Gauge_Conf *GC, GParam const *param, int const fields)
	{
	#ifdef OPENACC_MODE
	select_device(GC->device);
	acc_wait_all(); // the copy is synchronous: wait for the kernels in the asynchronous queues
	update_gauge_conf(GC, param, fields, 0);
	#else
	(void) GC;
	(void) param;
	(void) fields;
	#endif
	}


void copy_gauge_conf_to_device(Gauge_Conf *GC, GParam const *param, int const fields)
	{
	#ifdef OPENACC_MODE
	select_device(GC->device);
	acc_wait_all();
	update_gauge_conf(GC, param, fields, 1);
	#else
	(void) GC;
	(void) param;
	(void) fields;
	#endif
	}


void copy_gauge_conf_replica_to_host(Gauge_Conf *GC, GParam const *param, int const fields)
	{
	for(int g = 0; g < param->d_N_replica_pt; g++)
		copy_gauge_conf_to_host(&(GC[g]), param, fields);
	}


void copy_gauge_conf_replica_to_device(Gauge_Conf *GC, GParam const *param, int const fields)
	{
	for(int g = 0; g < param->d_N_replica_pt; g++)
		copy_gauge_conf_to_device(&(GC[g]), param, fields);
	}


void copy_replica_index_to_device(Gauge_Conf *GC, GParam const *param)
	{
	#ifdef OPENACC_MODE
	// all the compute regions run in QUEUE_DEVICE, so the copy is ordered after the kernels still reading the previous
	// value and before the next ones; the host changes replica_index again only in the next swap, after waiting for
	// the queues of all the devices (swap sums)
	for(int g = 0; g < param->d_N_replica_pt; g++)
		{
		select_device(GC[g].device);
		acc_update_device_async(&(GC[g].replica_index), sizeof(int), QUEUE_DEVICE);
		}
	#else
	(void) GC;
	(void) param;
	#endif
	}


double *alloc_double_on_devices(long const n)
	{
	double *a;
	allocate_array_double(&a, n, __FILE__, __LINE__);
	#ifdef OPENACC_MODE
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		acc_create(a, (size_t) n * sizeof(double));
		}
	#endif
	return a;
	}


void free_double_on_devices(double *a, long const n)
	{
	if(a == NULL) return;
	#ifdef OPENACC_MODE
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		if(acc_is_present(a, (size_t) n * sizeof(double)))
			{
			acc_wait_all();
			acc_delete(a, (size_t) n * sizeof(double));
			}
		}
	#else
	(void) n;
	#endif
	free(a);
	}


void copy_topo_grid_to_device(GParam *param)
	{
	#if defined(OPENACC_MODE) && defined(MULTICANONICAL_MODE)
	FOR_EACH_DEVICE(d)
		{
		select_device(d);
		acc_wait_all();
		update_2d((void **) param->d_grid, (size_t) param->d_N_replica_pt * (size_t) param->d_n_grid * sizeof(double), 1);
		}
	#else
	(void) param;
	#endif
	}


void swap_lattice_cold(Gauge_Conf *GC)
	{
	GAUGE_GROUP **aux = GC->lattice_cold;
	GC->lattice_cold = GC->lattice_copy_cold;
	GC->lattice_copy_cold = aux;

	#ifdef OPENACC_MODE
	// same exchange on the device copy of GC (not present during the initialization and for configurations not on
	// the devices)
	if(g_dev.n_used > 0)
		{
		select_device(GC->device);
		if(acc_is_present(GC, sizeof(Gauge_Conf)))
			{
			// in the queue of the device: ordered with the kernels using the cold lattices
			#pragma acc serial present(GC[0:1]) async(QUEUE_DEVICE)
				{
				GAUGE_GROUP **tmp = GC->lattice_cold;
				GC->lattice_cold = GC->lattice_copy_cold;
				GC->lattice_copy_cold = tmp;
				}
			}
		}
	#endif
	}

#endif
