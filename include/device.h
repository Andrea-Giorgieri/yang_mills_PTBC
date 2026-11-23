#ifndef DEVICE_H
#define DEVICE_H

#include "macro.h"

#include <stddef.h>
#include <stdio.h>

#include "gparam.h"

// Distribution of the configurations (parallel tempering replicas, or the single configuration of the other programs)
// among the devices (GPUs).
//
// The configurations GC[0..n_conf-1] (slots) are split in contiguous blocks: device d holds the slots
// first[d], ..., first[d] + n_slots[d] - 1 for the whole run (configurations never move, the swaps exchange the
// replica indices). The minimal number of devices able to host all the replicas is used, then the slots are
// distributed in a balanced way (e.g. 10 replicas, at most 4 per device -> 3 devices with 4, 3, 3 slots).
//
// CPU (no OpenACC): one device with all the slots, unless num_gpus > 1 in the input file, which simulates the
// splitting among num_gpus devices (to test the per-device loops).

#define DEVICE_MEMORY_FRACTION 0.9                         // fraction of the free memory of a device that can be used
#define DEVICE_MEMORY_OVERHEAD ((size_t) 512 * 1024 * 1024) // memory reserved for the runtime on every device (bytes)

typedef struct Device_Utils
	{
	int n_used;               // number of devices used
	int n_available;          // number of devices available (at most MAX_DEVICES)
	int id[MAX_DEVICES];      // device number (for acc_set_device_num) of the d-th device used
	int first[MAX_DEVICES];   // first slot of the d-th device used
	int n_slots[MAX_DEVICES]; // number of slots of the d-th device used
	int n_conf;               // total number of slots (configurations)
	size_t bytes_per_replica; // estimated memory of a replica
	size_t bytes_shared;      // estimated memory of the data needed on every device used
	} Device_Utils;

extern Device_Utils g_dev; // devices used in this run, set once by init_device_utils


// estimated memory (bytes) of the data of one replica (lattices, Z, clovers, random number generators)
size_t bytes_per_replica(GParam const *const param);

// estimated memory (bytes) of the data needed on every device used (geometry, auxiliary lattices for measurements)
size_t bytes_shared(GParam const *const param);

// distribute n_replica slots among the devices with free memory free_mem[0..n_available-1]:
// num_forced = 0 -> minimal number of devices, num_forced > 0 -> exactly num_forced devices
// return 0 on success, 1 if the replicas do not fit in the memory, 2 if num_forced is not valid
int assign_devices(Device_Utils *const dev, int const n_replica, int const n_available, size_t const *const free_mem,
                   size_t const per_replica, size_t const shared, int const num_forced);

// set g_dev for this run with n_conf configurations (N_replica_pt for the parallel tempering programs, 1 for the
// programs with a single configuration); queries the GPUs with OpenACC, simulated devices on CPU
void init_device_utils(GParam const *const param, int n_conf);

void print_device_utils(FILE *fp, Device_Utils const *const dev);


// loop over the devices used; the slots of device d are device_slot_begin(d) <= i < device_slot_end(d)
#define FOR_EACH_DEVICE(d) for(int d = 0; d < g_dev.n_used; d++)

static inline long device_slot_begin(int const d)
	{
	return (long) g_dev.first[d];
	}

static inline long device_slot_end(int const d)
	{
	return (long) (g_dev.first[d] + g_dev.n_slots[d]);
	}

// total number of slots (configurations on the devices)
static inline long device_num_slots(void)
	{
	return (long) g_dev.n_conf;
	}

#endif
