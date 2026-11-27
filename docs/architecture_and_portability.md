# Architecture and portability

This page describes how the same source runs on CPUs (OpenMP) and on NVIDIA GPUs (OpenACC), how data are laid out and
moved, how the lattice is parallelized (including odd lattice sizes) and how replicas are distributed among GPUs.

## Two exclusive parallel modes, one source

| Mode | Selection (compile time) | Target |
|---|---|---|
| Serial | neither option | one CPU core |
| OpenMP | `--enable-use-openmp`, `Num_threads=n` | shared-memory CPU node, n threads (fixed at compile time) |
| OpenACC | `--enable-use-openacc`, `OPENACC_CFLAGS` | NVIDIA GPUs, up to `Num_devices` per process; also the host CPU with `-acc=host`/`-acc=multicore` (for testing) |

OpenMP and OpenACC cannot be combined (configure stops). In OpenACC mode one host thread drives all the GPUs.

Every parallel loop carries both directives, each active only in its mode. Typical pattern
([lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c), `update_sites_parallel`):

```c
// dev = gc_on_device(GC), computed by the caller
#ifdef OPENACC_MODE
#pragma acc parallel loop gang vector vector_length(ACC_VL_UPD) if(dev) present(GC[0:1], geo[0:1], param[0:1]) async(QUEUE_DEVICE)
#endif
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS)
#endif
for(long s = begin; s < end; s++)
	{
	if(over) overrelaxation(GC, geo, param, geo->d_sip_to_si[s], dir);
	else heatbath(GC, geo, param, geo->d_sip_to_si[s], dir);
	}
```

Every OpenACC compute loop has an explicit `gang vector` schedule (nested loops in the called routines are then
sequential for every compiler and N; with nvc configure also adds `-acc=noautopar`). The vector length of each family
of kernels is a macro of [include/accel.h](../include/accel.h), which can be overridden with `-D`: `ACC_VL_UPD` = 32
(updates, cooling, clovers), `ACC_VL_FLOW` = 128 (gradient flow), `ACC_VL_MEAS` = 128 (measurements and the
reductions that drive decisions).

- Functions called inside compute regions are marked `ACC_ROUTINE_SEQ` (`#pragma acc routine seq` in OpenACC mode,
  nothing otherwise), defined in [include/accel.h](../include/accel.h). This covers the group algebra
  ([include/gauge_group.h](../include/gauge_group.h)), staples, single-link updates, clovers, the local observables
  and the random number generator.
- `gc_on_device(GC)` ([lib/device_data.c](../lib/device_data.c)) returns 1 when the configuration has a copy on its
  GPU, and makes that GPU current. The `if(dev)` clause then runs the same loop on the GPU or on the host. On CPU builds
  it always returns 0. `force_host_compute(1)` forces the host path (used by tests and by the sequential cooling
  orderings).
- Reductions (observables, charge differences) use `reduction(+:...)` in both modes; the host waits for the device
  queue only when it needs the scalar result (`wait_device_queue`). Exception: on the GPUs the swap sums are computed
  in blocks, one kernel per device, and the block sums are added on the host in a fixed order
  ([swap_translation_gpu.md](swap_translation_gpu.md), in Italian).
- Kernels are launched asynchronously in the queue `QUEUE_DEVICE` of each GPU; the host synchronizes only for
  swaps, multicanonical tests, measurements and I/O.

## Data layout

| Structure | Content | File |
|---|---|---|
| `Gauge_Conf` | one configuration: `lattice[site][μ]`, `lattice_copy` (state before the step, for restore and smoothing), `Z[site][plane]` and `Z_copy` (boundary factors and θ profile), `clover_array` (θ term), `lattice_cold`, `lattice_copy_cold` (multicanonical cooled copies), `rng_site` (Philox), `replica_index`, `device`, `translation`, `update_index`, `stored_topcharge` | [include/gauge_conf.h](../include/gauge_conf.h) |
| `Geometry` | neighbour tables `d_nnp`, `d_nnm`, update order `d_sip_to_si`, coordinate tables, defect mask | [include/geometry.h](../include/geometry.h) |
| `GParam` | input parameters and derived constants | [include/gparam.h](../include/gparam.h) |
| `Rect_Utils` | site lists of the PTBC rectangles | [include/geometry.h](../include/geometry.h) |

- Link matrices are stored as structures: U(1) one complex number, SU(2) four real numbers
  (a_0 + i a·σ), SU(N) N×N complex numbers in row-major order (`comp[m(i,j)]`), aligned to 32 bytes
  (`DOUBLE_ALIGN`). The alignment also rounds `sizeof` up to a multiple of 32 bytes: U(1) 32 B (16 of data), SU(2)
  32 B, SU(3) 160 B (144), SU(5) 416 B (400); a 4D SU(3) site row is 640 B. A different storage (AoSoA with blocks of
  `Link_block` sites) is planned: [link_layout.md](link_layout.md) (in Italian).
- Links are stored site by site, all directions of a site together (`lattice[r][μ]`).
- Every 2D array (`lattice`, `Z`, neighbour tables, …) is one contiguous memory block plus a table of row pointers
  (`allocate_array_*_2d` in [lib/memalign.c](../lib/memalign.c)). This makes host–device copies single block
  transfers.
- The site index used in memory is the even/odd lexicographic index (`cart_to_lexeo` in
  [lib/geometry.c](../lib/geometry.c)): lexicographic order with direction 0 running fastest, even sites first.

## Device memory management (OpenACC)

Implemented in [lib/device_data.c](../lib/device_data.c) and [include/device_data.h](../include/device_data.h). Without
OpenACC all these functions do nothing, so the programs are the same in every build.

- `init_device_data` copies to **every** GPU in use the shared data (`param`, `geo`, rectangles, auxiliary lattices of
  the measurements) and, to the GPU of each slot, the configurations of that GPU. It is called once after
  initialization; `free_device_data` releases everything at the end.
- Pointer tables are copied once and translated into device pointers on the device (`table_to_device`), then attached
  to the structure members (`acc_attach`). Whole structures are never updated after that (it would overwrite device
  pointers with host ones): only data blocks and scalar members are copied.
- Explicit copies use bit masks of fields: `GC_LATTICE`, `GC_COPY`, `GC_Z`, `GC_COLD`, `GC_RNG`
  (`copy_gauge_conf_to_host`, `copy_gauge_conf_to_device` and the `_replica_` variants). The programs copy the
  configurations to the host only to save them or for the observables computed on the host.
- "Views" are persistent auxiliary `Gauge_Conf` headers whose pointers alias other fields of a slot (cooled copy,
  saved copy, flow lattice), so that the same kernels can act on them on the device (`gauge_conf_view`).
- The Philox generator states live on the device and are never copied back during the run.

What runs where in OpenACC builds:

| On the GPU | On the host |
|---|---|
| heatbath, overrelaxation (full lattice, rectangles, border sites), accept/restore, cooling of the multicanonical charge, swap sums, translations, plaquette, clover energy, topological charge, Polyakov loops, gradient flow (fixed and adaptive), cooling orderings `LEX_DIR_LEXEO_SITE`, `LEXEO_SITE_LEX_DIR`, `RND_DIR_RNDEO_SITE` | Metropolis tests of swaps and multicanonical steps, random choices of translations and cooling orders; observables that write densities or slices, `charge_prime`, `chi_prime`, `multipolyakov`, `polyakov_powers`, `action` (the configuration is copied to the host before them); sequential cooling orderings; I/O |

Programs with data on the GPU: `yang_mills_local_pt_agf`, `yang_mills_tuning_pt_mc`, `yang_mills_local_agf`,
`yang_mills_conf_analysis`. `yang_mills_multilevel*` and `yang_mills_tracedef` are not built with OpenACC.

## Several GPUs: distribution of the replicas

Implemented in [lib/device.c](../lib/device.c) and [include/device.h](../include/device.h).

- The configurations (PTBC replicas, or the single configuration of the other programs) are **slots** `GC[0 … n−1]`.
  Each slot lives on one GPU for the whole run; swaps exchange replica indices, not data, so configurations never move
  between GPUs.
- The slots are split in contiguous, balanced blocks. By default the code uses the **smallest** number of GPUs whose
  memory can hold the replicas (e.g. 10 replicas with at most 4 per GPU → 3 GPUs with 4, 3, 3). `num_gpus k` (runtime)
  forces exactly k GPUs.
- Memory estimate: `bytes_per_replica` (lattice and copies, Z, clovers, cold lattices, generators) and `bytes_shared`
  (geometry, rectangles, potential, auxiliary lattices for smoothing). Usable memory per GPU = 90% of the free memory
  minus 512 MiB and the shared data (`DEVICE_MEMORY_FRACTION`, `DEVICE_MEMORY_OVERHEAD`).
- At most `Num_devices` GPUs (compile time, default 4) are used: the first visible ones (device numbers 0 … k−1), so
  the GPUs can be selected with `CUDA_VISIBLE_DEVICES`.
- If the replicas do not fit, the program stops with an error. `yang_mills_gpu_plan input_file [GB_per_GPU]` (default
  64 GB) computes the same estimate without a GPU and prints the distribution and the number of GPUs needed.
- On CPU builds (and OpenACC on the host) memory is not limited and `num_gpus > 1` only simulates the splitting, to
  test the per-device loops.
- Loops over replicas run, for each GPU, one kernel over (its slots × sites) (`FOR_EACH_DEVICE`,
  `update_replica_sites` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c)). Between GPUs only scalars are exchanged
  (swap sums, charges), through the host.

## Lattice parallelization: even/odd and odd sizes

The Wilson action couples a link only to links of neighbouring sites, so all links U_μ(x) with the same direction μ
and the same parity of x can be updated at the same time. The code orders the sites so that this can be done with
simple parallel loops ([lib/geometry.c](../lib/geometry.c), comment before `shell_volume`):

- **Even sizes.** Sites are visited as all even sites, then all odd sites; since memory uses the same order,
  `d_sip_to_si` is then the identity.
- **Odd sizes.** Along an odd side L, the sites x = 0 and x = L−1 are neighbours and have the same parity, so the
  naive even/odd split is wrong. The lattice is divided into "shells": the largest sub-lattice with even sides (the
  first L−1 sites along each odd direction) and the border shells (x = L−1 along one or more odd directions). The even
  sub-lattice is updated with the even/odd scheme in parallel; the border sites are updated afterwards **sequentially**
  (on the GPU: one gang with one thread, `num_gangs(1) vector_length(1)`; for PTBC, one thread per slot, each running
  the border sites of its configuration in sequence, `update_replica_border`). This keeps
  the update correct for any size, at the cost of a sequential part proportional to the border.
- The same scheme is used by heatbath, overrelaxation and cooling.
- **PTBC update rectangles** are ordered even/odd inside the rectangle and updated in parallel; with OpenMP or OpenACC
  their size is limited to the even sides (see [ptbc.md](ptbc.md#rectangles)).
- **Trace deformation** (`yang_mills_tracedef`) sweeps the spatial lattice with its own even/odd split, which does not
  handle odd sides: with OpenMP all spatial sizes must be even when a coupling h_i ≠ 0.

Functions: `cart_to_lexeop`, `shell_*` (ordering), `init_geometry` ([lib/geometry.c](../lib/geometry.c));
`update_sites_parallel`, `update_border_sites`, `update_replica_border` ([lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c)).

## Random numbers and reproducibility

- OpenACC requires a generator that can be used independently by every thread: the Philox4x32-10 counter-based
  generator, with one stream per (configuration slot, site), is always used with OpenACC and optional on CPU
  (`--enable-rng-philox`). Its random numbers do not depend on how sites are assigned to threads or GPUs.
- The default CPU generator is dSFMT, with one state per OpenMP thread.
- Decisions taken on the host (swaps, translations, multicanonical tests, random starts) use `rng_host()`: one host
  stream with Philox, the generator of the calling thread with dSFMT (details in [algorithms.md](algorithms.md)).

Files: [include/random.h](../include/random.h), [lib/random.c](../lib/random.c).

## Portability notes

- Compilers used: gcc, Intel icc (CPU presets for Leonardo, Galileo100, Marconi), NVIDIA HPC SDK nvc (CPU and GPU).
  Presets in [`build.sh`](../build.sh); see [building.md](building.md).
- With nvc, configure appends `-Mvect=nofuse` (works around a wrong loop fusion of nvc in consecutive inlined matrix
  functions).
- Binary configuration files are written in big-endian order on every machine (portable across architectures), with an
  optional MD5 hash in the header ([output_files.md](output_files.md)).
- Timers use `omp_get_wtime` with OpenMP and `clock_gettime(CLOCK_MONOTONIC)` otherwise ([lib/timing.c](../lib/timing.c)).
