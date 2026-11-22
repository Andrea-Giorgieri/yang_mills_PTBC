#ifndef DEVICE_C
#define DEVICE_C

#include "../include/macro.h"
#include "../include/device.h"

#include <complex.h>
#include <stdio.h>
#ifdef OPENACC_MODE
#include <openacc.h>
#endif

#include "../include/gauge_group.h"
#include "../include/gparam.h"
#include "../include/random.h"

Device_Utils g_dev; // devices used in this run


// estimated memory of an array [n1][n2] allocated with allocate_array_*_2d: block + table of pointers
static size_t bytes_2d(long const n1, long const n2, size_t const elem_size)
	{
	return (size_t) n1 * (size_t) n2 * elem_size + (size_t) n1 * sizeof(void *);
	}


size_t bytes_per_replica(GParam const *const param)
	{
	long const V = param->d_volume;
	size_t bytes = 0;

	bytes += 2 * bytes_2d(V, STDIM, sizeof(GAUGE_GROUP));                     // lattice, lattice_copy
	bytes += 2 * bytes_2d(V, param->d_n_planes + 1, sizeof(double complex));  // Z, Z_copy
	#ifdef MULTICANONICAL_MODE
	bytes += 2 * bytes_2d(V, STDIM, sizeof(GAUGE_GROUP));  // lattice_cold, lattice_copy_cold
	#endif
	#ifdef THETA_MODE
	bytes += bytes_2d(V * STDIM, STDIM, sizeof(GAUGE_GROUP)) + (size_t) V * sizeof(void *);  // clover_array
	#endif
	#ifdef RNG_PHILOX_MODE
	bytes += (size_t) V * sizeof(rng_t);  // rng_site
	#endif

	return bytes;
	}


size_t bytes_shared(GParam const *const param)
	{
	long const V = param->d_volume;
	size_t bytes = 0;

	// geometry
	bytes += (size_t) V * sizeof(long);                   // d_sip_to_si
	bytes += 2 * bytes_2d(V, STDIM, sizeof(long));        // d_nnp, d_nnm
	bytes += bytes_2d(STDIM, V, sizeof(int));             // d_mucomp
	bytes += bytes_2d(STDIM, V, sizeof(long));            // d_muorth
	bytes += (size_t) STDIM * (size_t) V * sizeof(long);  // d_mutimespace (block)
	bytes += (size_t) V * sizeof(int);                    // d_defect_mask

	// auxiliary lattices for smoothing measurements, one set per device (the physical replica can be on any device)
	if(param->d_agf_num_meas > 0 || param->d_gf_num_meas > 0 || param->d_coolrepeat > 0)
		bytes += 4 * bytes_2d(V, STDIM, sizeof(GAUGE_GROUP));

	return bytes;
	}


int assign_devices(Device_Utils *const dev, int const n_replica, int const n_available, size_t const *const free_mem,
                   size_t const per_replica, size_t const shared, int const num_forced)
	{
	dev->n_used = 0;
	dev->n_available = (n_available < MAX_DEVICES) ? n_available : MAX_DEVICES;
	dev->bytes_per_replica = per_replica;
	dev->bytes_shared = shared;

	if(num_forced < 0 || num_forced > dev->n_available || num_forced > n_replica) return 2;

	// maximum number of replicas that fit in each device
	long max_slots[MAX_DEVICES];
	for(int d = 0; d < dev->n_available; d++)
		{
		double const usable = DEVICE_MEMORY_FRACTION * (double) free_mem[d] - (double) DEVICE_MEMORY_OVERHEAD - (double) shared;
		max_slots[d] = (usable > 0.0) ? (long) (usable / (double) per_replica) : 0;
		}

	// smallest number of devices (or the forced one) such that the balanced distribution fits
	int const n_min = (num_forced > 0) ? num_forced : 1;
	int const n_max = (num_forced > 0) ? num_forced : dev->n_available;
	for(int n = n_min; n <= n_max && n <= n_replica; n++)
		{
		int fits = 1;
		for(int d = 0; d < n; d++)
			{
			long const slots = n_replica / n + ((d < n_replica % n) ? 1 : 0); // first devices take the remainder
			if(slots > max_slots[d]) fits = 0;
			}
		if(fits)
			{
			dev->n_used = n;
			int first = 0;
			for(int d = 0; d < n; d++)
				{
				dev->id[d] = d;
				dev->first[d] = first;
				dev->n_slots[d] = n_replica / n + ((d < n_replica % n) ? 1 : 0);
				first += dev->n_slots[d];
				}
			return 0;
			}
		}
	return 1;
	}


void init_device_utils(GParam const *const param)
	{
	size_t free_mem[MAX_DEVICES];
	int n_available;

	#ifdef OPENACC_MODE
	n_available = acc_get_num_devices(acc_device_nvidia);
	REQUIRE(n_available > 0, "no NVIDIA GPU available");
	if(n_available > MAX_DEVICES) n_available = MAX_DEVICES;
	for(int d = 0; d < n_available; d++)
		free_mem[d] = acc_get_property(d, acc_device_nvidia, acc_property_free_memory);
	#else
	// CPU: memory is not limited by the devices, num_gpus > 1 only simulates the splitting
	n_available = MAX_DEVICES;
	for(int d = 0; d < n_available; d++)
		free_mem[d] = (size_t) -1;
	#endif

	int const err = assign_devices(&g_dev, param->d_N_replica_pt, n_available, free_mem,
	                               bytes_per_replica(param), bytes_shared(param), param->d_num_gpus);
	REQUIRE(err != 2, "invalid num_gpus = %d: must be in [0, %d] and not larger than the number of replicas %d",
	        param->d_num_gpus, (n_available < MAX_DEVICES) ? n_available : MAX_DEVICES, param->d_N_replica_pt);
	REQUIRE(err != 1, "the %d replicas (%.3f GB each, plus %.3f GB per device) do not fit in the memory of %d device(s)",
	        param->d_N_replica_pt, (double) bytes_per_replica(param) / 1e9, (double) bytes_shared(param) / 1e9, n_available);
	}


void print_device_utils(FILE *fp, Device_Utils const *const dev)
	{
	#ifdef OPENACC_MODE
	fprintf(fp, "Devices (GPUs) used: %d of %d available\n", dev->n_used, dev->n_available);
	#else
	fprintf(fp, "Devices (simulated on CPU) used: %d\n", dev->n_used);
	#endif
	fprintf(fp, "estimated memory per replica: %.3f GB, per device (shared data): %.3f GB\n",
	        (double) dev->bytes_per_replica / 1e9, (double) dev->bytes_shared / 1e9);
	for(int d = 0; d < dev->n_used; d++)
		{
		fprintf(fp, "  device %d: slots %d - %d (%d replicas)\n",
		        dev->id[d], dev->first[d], dev->first[d] + dev->n_slots[d] - 1, dev->n_slots[d]);
		}
	fprintf(fp, "\n");
	}

#endif
