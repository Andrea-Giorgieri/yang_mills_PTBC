#ifndef DEVICE_DATA_H
#define DEVICE_DATA_H

#include "macro.h"

#include "device.h"
#include "gauge_conf.h"
#include "geometry.h"
#include "gparam.h"

// Data on the devices (OpenACC): every used device holds a copy of the shared data (param, geo, rect_aux, auxiliary
// lattices of the measurements) and the configurations of its own slots GC[first : first + n_slots].
//
// Pointer tables (contiguous arrays of memalign.h) are copied to the device and translated into device pointers,
// then attached to the member of the structure. Never update a whole structure (GC[g], param, geo, ...) or a pointer
// table: host pointers would overwrite the device ones. Only blocks and scalar members are updated, with the
// functions below.
//
// Without OPENACC_MODE all the functions do nothing (swap_lattice_cold only swaps on the host, the views live on the
// host).

// fields of a configuration to be copied between host and device (bit mask)
enum
	{
	GC_LATTICE = 1 << 0, // lattice
	GC_COPY = 1 << 1,    // lattice_copy, Z_copy
	GC_Z = 1 << 2,       // Z
	GC_COLD = 1 << 3,    // lattice_cold, lattice_copy_cold (multicanonical)
	GC_RNG = 1 << 4,     // rng_site (Philox)
	GC_ALL = GC_LATTICE | GC_COPY | GC_Z | GC_COLD | GC_RNG,
	GC_ALL_BUT_RNG = GC_LATTICE | GC_COPY | GC_Z | GC_COLD // the generators are updated only on the device (heatbath)
	};

// asynchronous queues: QUEUE_DEVICE for the kernels over all the slots of a device, queue_slot(g) for the kernels of a
// single slot g (reductions); the copy functions wait for all the queues of the device before copying
#define QUEUE_DEVICE 1

static inline int queue_slot(long const g)
	{
	return 2 + (int) g;
	}

// make d (index in g_dev) the current device
void select_device(int d);

// 1 if the compute regions acting on GC run on the device holding it (OPENACC_MODE, init_device_data called, GC present
// on its device and the host not forced), and that device becomes the current one; 0 if they run on the host.
// The functions written for both use it in the clause if() of their compute regions.
int gc_on_device(Gauge_Conf const *GC);

// if dev (result of gc_on_device), wait for the queue QUEUE_DEVICE of the current device (results of reductions)
void wait_device_queue(int dev);

// 1: the compute regions of the functions using gc_on_device run on the host (data of the host), 0: default;
// returns the previous value
int force_host_compute(int on);

// element of meas_aux given to init_device_data (its auxiliary lattices are on every device), NULL before
Meas_Utils *device_meas_aux(void);

// wait for all the queues of all the used devices
void wait_devices(void);

// copy the shared data to every used device and the configurations GC[0 : device_num_slots()] to the device of their
// slot (after init_device_utils, the initialization of the configurations and the other initializations);
// meas_aux[0 : n_meas] share the auxiliary lattices of meas_aux[0] (n_meas = num_meas_utils_replica(param) for the array
// of init_meas_utils_replica); rect_aux = NULL for the programs without parallel tempering; configurations with
// rng_site = NULL (never updated) have no generators on the device
void init_device_data(Gauge_Conf *GC, Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int n_meas);

// release the device copies (before freeing the host data)
void free_device_data(Gauge_Conf *GC, Geometry *geo, GParam *param, Rect_Utils *rect_aux, Meas_Utils *meas_aux, int n_meas);

// copy the fields of a configuration from/to the device holding it
void copy_gauge_conf_to_host(Gauge_Conf *GC, GParam const *param, int fields);
void copy_gauge_conf_to_device(Gauge_Conf *GC, GParam const *param, int fields);

// the same for all the configurations GC[0 : N_replica_pt]
void copy_gauge_conf_replica_to_host(Gauge_Conf *GC, GParam const *param, int fields);
void copy_gauge_conf_replica_to_device(Gauge_Conf *GC, GParam const *param, int fields);

// copy replica_index of all the configurations to their devices (after the swaps)
void copy_replica_index_to_device(Gauge_Conf *GC, GParam const *param);

// copy the topological potential d_grid to every used device (after tune_topo_potential)
void copy_topo_grid_to_device(GParam *param);

// exchange lattice_cold and lattice_copy_cold, on the host and on the device holding GC (if present there)
void swap_lattice_cold(Gauge_Conf *GC);

// views: persistent auxiliary configurations, one per slot and kind, whose lattice and Z point to other fields of the
// configuration of the slot (allocated by init_device_data in every build, on the devices with OPENACC_MODE)
typedef enum
	{
	VIEW_COLD,      // lattice = lattice_cold, Z = Z, lattice_copy = lattice_copy_cold, Z_copy = Z_copy
	VIEW_COPY,      // lattice = lattice_copy, Z = Z_copy
	VIEW_COLD_COPY, // lattice = lattice_copy_cold, Z = Z_copy
	VIEW_FLOW,      // lattice = device_meas_aux()->lattice_aux[0] (on the device of the slot), Z = Z (gradient flow)
	NUM_VIEWS
	} View_Kind;

// view of the given kind of GC (an element of the array given to init_device_data); the pointers are assigned on the
// host and, in the queue QUEUE_DEVICE, on the device holding GC: call it again after swap_lattice_cold
Gauge_Conf *gauge_conf_view(Gauge_Conf const *GC, View_Kind kind);

#endif
