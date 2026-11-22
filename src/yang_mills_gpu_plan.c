#ifndef YM_GPU_PLAN_C
#define YM_GPU_PLAN_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/device.h"
#include "../include/gparam.h"

// Estimate the memory needed by the replicas and the number of GPUs required, without using any GPU.
// The memory of every GPU is given on the command line.

static double const default_gpu_memory_gb = 64.0;


void real_main(char const *in_file, double const gpu_memory_gb)
	{
	GParam param;
	readinput(in_file, &param);

	size_t free_mem[MAX_DEVICES];
	for(int d = 0; d < MAX_DEVICES; d++)
		free_mem[d] = (size_t) (gpu_memory_gb * 1e9);

	size_t const per_replica = bytes_per_replica(&param);
	size_t const shared = bytes_shared(&param);
	double const usable = DEVICE_MEMORY_FRACTION * gpu_memory_gb * 1e9 - (double) DEVICE_MEMORY_OVERHEAD - (double) shared;
	long const max_per_device = (usable > 0.0) ? (long) (usable / (double) per_replica) : 0;

	printf("lattice: %d", param.d_size[0]);
	for(int i = 1; i < STDIM; i++) printf("x%d", param.d_size[i]);
	printf(", replicas: %d\n", param.d_N_replica_pt);
	printf("memory per GPU: %.1f GB (usable fraction %.2f, runtime overhead %.2f GB)\n",
	       gpu_memory_gb, DEVICE_MEMORY_FRACTION, (double) DEVICE_MEMORY_OVERHEAD / 1e9);
	printf("estimated memory per replica: %.3f GB\n", (double) per_replica / 1e9);
	printf("estimated memory per GPU (shared data): %.3f GB\n", (double) shared / 1e9);
	printf("maximum replicas per GPU: %ld\n\n", max_per_device);

	Device_Utils dev;
	int const err = assign_devices(&dev, param.d_N_replica_pt, MAX_DEVICES, free_mem, per_replica, shared, param.d_num_gpus);
	if(err == 0)
		{
		print_device_utils(stdout, &dev);
		printf("=> Necessary GPUs: %d\n", dev.n_used);
		}
	else if(err == 1)
		{
		printf("the replicas do not fit in %d GPUs of %.1f GB%s\n", MAX_DEVICES, gpu_memory_gb,
		       (param.d_num_gpus > 0) ? " with the requested num_gpus" : "");
		}
	else
		{
		printf("invalid num_gpus = %d (must be in [0, %d] and not larger than the number of replicas)\n", param.d_num_gpus, MAX_DEVICES);
		}

	free_gparam(&param);
	}


int main(int argc, char **argv)
	{
	if(argc != 2 && argc != 3)
		{
		printf("Usage: %s input_file [memory_per_gpu_GB (default %.0f)]\n\n", argv[0], default_gpu_memory_gb);
		printf("Estimate the memory needed by the parallel tempering replicas and the number of GPUs required.\n");
		printf("The input file is the one of yang_mills_local_pt_agf or yang_mills_tuning_pt_mc.\n");
		printf("At most %d GPUs per node (Num_devices at configure time).\n", MAX_DEVICES);
		return EXIT_SUCCESS;
		}

	REQUIRE(strlen(argv[1]) < STD_STRING_LENGTH, "input filename too long, increase STD_STRING_LENGTH in macro.h");
	double const gpu_memory_gb = (argc == 3) ? atof(argv[2]) : default_gpu_memory_gb;
	REQUIRE(gpu_memory_gb > 0.0, "invalid memory per GPU: %s", argv[2]);

	real_main(argv[1], gpu_memory_gb);

	return EXIT_SUCCESS;
	}

#endif
