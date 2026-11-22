#ifndef DEBUG_DEVICE_C
#define DEBUG_DEVICE_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>

#include "../include/device.h"

// tests of the distribution of the replicas among the devices with fictitious memories


static size_t const GB = (size_t) 1000 * 1000 * 1000;


// check the result of assign_devices against the expected status and slots per device (expected_slots ends with 0)
static int check_case(char const *name, int n_replica, int n_available, size_t mem_gb, size_t per_replica, size_t shared,
                      int num_forced, int expected_err, int const *expected_slots)
	{
	size_t free_mem[MAX_DEVICES];
	for(int d = 0; d < MAX_DEVICES; d++) free_mem[d] = mem_gb * GB;

	Device_Utils dev;
	int const err = assign_devices(&dev, n_replica, n_available, free_mem, per_replica, shared, num_forced);

	int ok = (err == expected_err);
	if(ok && err == 0)
		{
		int n = 0;
		while(expected_slots[n] != 0) n++;
		ok = (dev.n_used == n);
		int first = 0;
		for(int d = 0; ok && d < n; d++)
			{
			ok = (dev.n_slots[d] == expected_slots[d] && dev.first[d] == first && dev.id[d] == d);
			first += dev.n_slots[d];
			}
		ok = ok && (first == n_replica);
		}

	printf("%-52s err=%d n_used=%d slots=", name, err, (err == 0) ? dev.n_used : 0);
	for(int d = 0; err == 0 && d < dev.n_used; d++) printf("%d ", dev.n_slots[d]);
	printf("  %s\n", ok ? "OK" : "FAILED");
	return ok;
	}


int main(void)
	{
	int ok = 1;
	printf("MAX_DEVICES = %d\n\n", MAX_DEVICES);

	// usable memory of a 64 GB device: 0.9 * 64 - 0.537 = 57.06 GB minus the shared data
	{ int const s[] = {10, 0};       ok = check_case("10 replicas of 5 GB, 1 GB shared: one device", 10, 4, 64, 5 * GB, 1 * GB, 0, 0, s) && ok; }
	{ int const s[] = {4, 3, 3, 0};  ok = check_case("10 replicas, at most 4 per device: 4,3,3", 10, 4, 64, 13 * GB, 1 * GB, 0, 0, s) && ok; }
	{ int const s[] = {3, 3, 3, 3, 0}; ok = check_case("12 replicas, at most 3 per device: 3,3,3,3", 12, 4, 64, 18 * GB, 1 * GB, 0, 0, s) && ok; }
	{ int const s[] = {0};           ok = check_case("13 replicas, at most 3 per device: do not fit", 13, 4, 64, 18 * GB, 1 * GB, 0, 1, s) && ok; }
	{ int const s[] = {0};           ok = check_case("shared data larger than the device: do not fit", 2, 4, 64, 1 * GB, 60 * GB, 0, 1, s) && ok; }
	{ int const s[] = {0};           ok = check_case("only 2 devices available, 3 needed: do not fit", 10, 2, 64, 13 * GB, 1 * GB, 0, 1, s) && ok; }
	{ int const s[] = {3, 2, 0};     ok = check_case("num_gpus = 2 forced (1 would be enough): 3,2", 5, 4, 64, 1 * GB, 1 * GB, 2, 0, s) && ok; }
	{ int const s[] = {0};           ok = check_case("num_gpus = 2 forced but 3 needed: do not fit", 10, 4, 64, 13 * GB, 1 * GB, 2, 1, s) && ok; }
	{ int const s[] = {0};           ok = check_case("num_gpus larger than the devices: invalid", 10, 4, 64, 1 * GB, 1 * GB, MAX_DEVICES + 1, 2, s) && ok; }
	{ int const s[] = {0};           ok = check_case("num_gpus larger than the replicas: invalid", 2, 4, 64, 1 * GB, 1 * GB, 3, 2, s) && ok; }
	{ int const s[] = {1, 0};        ok = check_case("1 replica", 1, 4, 64, 1 * GB, 1 * GB, 0, 0, s) && ok; }
	{ int const s[] = {1, 1, 1, 1, 0}; ok = check_case("4 replicas, 1 per device: 1,1,1,1", 4, 4, 64, 40 * GB, 1 * GB, 0, 0, s) && ok; }

	printf("\n%s\n", ok ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
	}

#endif
