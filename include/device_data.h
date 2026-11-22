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
// Without OPENACC_MODE all the functions do nothing (swap_lattice_cold only swaps on the host).

// fields of a configuration to be copied between host and device (bit mask)
enum
	{
	GC_LATTICE = 1 << 0, // lattice
	GC_COPY = 1 << 1,    // lattice_copy, Z_copy
	GC_Z = 1 << 2,       // Z
	GC_COLD = 1 << 3,    // lattice_cold, lattice_copy_cold (multicanonical)
	GC_RNG = 1 << 4,     // rng_site (Philox)
	GC_ALL = GC_LATTICE | GC_COPY | GC_Z | GC_COLD | GC_RNG
	};

// make d (index in g_dev) the current device
void select_device(int d);

// copy the shared data to every used device and the configurations to the device of their slot
// (after init_gauge_conf_replica and the other initializations); meas_aux[0 : n_meas] share the auxiliary lattices of
// meas_aux[0] (n_meas = num_meas_utils_replica(param) for the array of init_meas_utils_replica)
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

#endif
