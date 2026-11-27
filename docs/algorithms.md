# Algorithms

This page describes the Monte Carlo algorithms, how they are combined in a Monte Carlo step of each program, and the
related options. Two algorithms have their own pages: parallel tempering on boundary conditions
([ptbc.md](ptbc.md)) and the multicanonical algorithm ([multicanonical.md](multicanonical.md)).

Convention: **compile time** = option of `configure`/`build.sh` ([building.md](building.md)); **runtime** = keyword in
the input file ([input_file.md](input_file.md)).

## Local link updates

All updates are local: a link U_μ(x) is updated in the field of its staple M (sum of the D−1 pairs of staples, each
multiplied by the boundary factors Z and, in PTBC, by the defect factors K; with the θ term, plus the "topological
staple" built from clovers). Implementation: staples in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c), single-link
updates in [include/u1_upd.h](../include/u1_upd.h), [include/su2_upd.h](../include/su2_upd.h),
[include/sun_upd.h](../include/sun_upd.h).

| Algorithm | Group | Method |
|---|---|---|
| Heatbath | U(1) | exact heatbath of the angle (Moriarty, Phys. Rev. D 25 (1982) 2185) |
| Heatbath | SU(2) | Kennedy–Pendleton (Phys. Lett. B 156 (1985) 393) for large coupling of the staple, Creutz otherwise (switch at k = 1.6847) |
| Heatbath | SU(N), N ≥ 3 | Cabibbo–Marinari pseudo-heatbath (Phys. Lett. B 119 (1982) 387): SU(2) heatbath on each of the N(N−1)/2 SU(2) subgroups |
| Overrelaxation | U(1), SU(2) | microcanonical reflection of the link with respect to the normalized staple |
| Overrelaxation | SU(N), N ≥ 3 | reflection on each SU(2) subgroup (Cabibbo–Marinari) |
| Metropolis | all | proposal U′ = R U or R† U (probability 1/2 each), R = projection on the group of (1 + ε X), X random; accept with min(1, e^{−ΔS}) |

- `overrelax` (runtime): number of overrelaxation sweeps after each heatbath sweep.
- `epsilon_metro` (runtime, default 0.25): ε of the Metropolis proposal. Metropolis is used only for the temporal links
  in `yang_mills_tracedef`.

### Order of a sweep and parallelization

A heatbath (or overrelaxation) sweep loops over the directions μ = 0 … D−1; for each direction it updates all the
links U_μ(x) of the even sites in parallel, then those of the odd sites in parallel (they share no staple). On lattices
with odd sides, the sites outside the largest even sub-lattice ("border" sites) are updated afterwards sequentially.
See [architecture_and_portability.md](architecture_and_portability.md) for details.

With the θ term, the clovers needed by the topological staples of direction μ are recomputed before each direction is
swept (`compute_clovers`).

## Monte Carlo step of each program

One **Monte Carlo step** (the unit of `sample`, `thermal`, `measevery`, `saveconf_*_every`) is:

| Program | One step |
|---|---|
| `yang_mills_local_agf`, `yang_mills_multilevel`, `yang_mills_multilevel_long` | 1 heatbath sweep + `overrelax` overrelaxation sweeps (`update` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c)); with the multicanonical algorithm, followed by an accept/reject test of the whole step |
| `yang_mills_local_pt_agf`, `yang_mills_tuning_pt_mc` | on every replica: full sweep as above, accept/reject (multicanonical), swaps, translation; then the hierarchical updates of the rectangles around the defect, each followed by swaps and translation (`parallel_tempering_with_hierarchical_update`; see [ptbc.md](ptbc.md)) |
| `yang_mills_tracedef` | heatbath of the spatial links; Metropolis of the temporal links, time slice by time slice; `overrelax` overrelaxation sweeps of the spatial links (`update_with_trace_def`) |

Measurements are taken when the update index is a multiple of `measevery` and at least `thermal`. With `sample 0`
the programs only measure the initial configuration.

## Parallel tempering on boundary conditions (PTBC)

Several replicas of the lattice, which differ only by the boundary conditions on a small region (the "defect"), are
evolved together and swapped with a Metropolis test. This reduces the autocorrelation of the topological charge.
Runtime: `N_replica_pt`, `defect_dir`, `defect_size`, `hierarc_upd`. Details: [ptbc.md](ptbc.md).

## Multicanonical algorithm

A bias potential V(Q) of the topological charge is added to the action and removed by reweighting. The potential can
be different for each PTBC replica and can be tuned automatically. Compile time: `--enable-use-multicanonical`;
runtime: `topo_potential_file`, `grid_step`, `grid_max`, `topo_cooling`, `topo_coolsteps`, `topo_alpha` and the
`topo_tuning_*` parameters. Details: [multicanonical.md](multicanonical.md).

## Smoothing (measurements only)

Gradient flow (fixed or adaptive step) and cooling are used to measure smoothed observables; the configuration is
restored afterwards, so they do not affect the Markov chain. Parameters and algorithms are described in
[physics_and_parameters.md](physics_and_parameters.md#smoothing-gradient-flow-and-cooling).
Implementation: `gradflow_RKstep`, `gradflow_RKstep_adaptive`, `cooling` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c).

## Multilevel algorithm

`yang_mills_multilevel` and `yang_mills_multilevel_long` measure Polyakov-loop correlators and flux tubes with the
multilevel algorithm: the time direction is divided into slabs, the links inside each slab are updated with the slab
boundaries fixed, and the slab averages of the tensor products of temporal links are multiplied together.

| Parameter | Meaning |
|---|---|
| `Num_levels` (compile time) | number of levels |
| `ml_step s_0 … s_{n−1}` | thickness of the slabs at each level, from the largest (level 0) to the smallest. `L_0` must be divisible by `s_0`, each `s_i` must divide `s_{i−1}` and be smaller, and the last must be > 1 |
| `ml_upd u_0 … u_{n−1}` | number of updates at each level |
| `multihit` | number of multihit steps for the temporal links (used when `dist_poly > 1` and `L_1 − dist_poly > 1`) |
| `ml_obs` | observable (see [physics_and_parameters.md](physics_and_parameters.md#multilevel-observables)) |
| `dist_poly`, `transv_dist`, `plaq_dir` | geometry of the observables |
| `ml_level0_repeat` | `yang_mills_multilevel_long`: number of program executions over which the level-0 average is accumulated |
| `ml_file` | file with the partial multilevel averages, written and read back by `yang_mills_multilevel_long` (the parameter is required whenever `Num_levels > 1`) |

Implementation: [lib/gauge_conf_multilevel.c](../lib/gauge_conf_multilevel.c),
`perform_multilevel_update_and_measures` in [lib/gauge_conf_meas.c](../lib/gauge_conf_meas.c).
Limits: CPU only (not built with OpenACC), not with the multicanonical algorithm, all spatial sizes equal.

> **To be verified.** This page describes the multilevel parameters only at the level of the input file; the
> detailed scheme (role of each level, which links are updated with `ml_upd`) has not been checked against
> [lib/gauge_conf_multilevel.c](../lib/gauge_conf_multilevel.c).

## Random numbers

| Generator | Selection | Streams |
|---|---|---|
| dSFMT (double-precision SIMD Mersenne Twister) | default on CPU | one generator per OpenMP thread |
| Philox4x32-10 (counter-based) | compile time `--enable-rng-philox`; always with OpenACC | one stream per lattice site and per configuration slot (`GC[i]`, stream i; the stream stays with the slot when swaps exchange the replica indices) for the link updates, plus one host stream for serial decisions |

- `randseed` (runtime): seed; `0` means "use the current time".
- The serial decisions (swaps, translations, multicanonical tests, random starts, random cooling orders) use
  `rng_host()`: with Philox a single host stream; with dSFMT the generator of the calling thread (`rng_thread()`).
- With Philox the random numbers used by a link update depend only on the seed, the slot and the site, not on the
  number of threads or GPUs. With dSFMT the stream of each thread depends on how sites are assigned to threads, so
  results change with `Num_threads`.
- With Philox and OpenMP, the few loops that draw from the host stream (initialization of the replicas, PTBC swaps)
  run serially.

Files: [include/random.h](../include/random.h), [lib/random.c](../lib/random.c), [lib/dSFMT.c](../lib/dSFMT.c).

## Starting configuration, checkpoints, restart

| `start` | Initial configuration |
|---|---|
| 0 | ordered: every link is the identity plus a small random perturbation (amplitude 0.001), projected on the group |
| 1 | random links |
| 2 | read from `conf_file` (and `twist_file`); for PTBC programs, from `conf_file_replica_<i>` and `twist_file_replica_<i>` |
| 3 | ordered, compatible with a twist on one plane (twist-eater matrices on the boundary links) |

- `saveconf_back_every n`: every n steps the configuration(s) are written to `conf_file` (and the twist file), and a
  backup copy alternately to `<file>_back0` and `<file>_back1`. The configuration is also written at the end of the run.
  `0` disables saving.
- `saveconf_analysis_every n`: every n steps the configuration (for PTBC: the physical one) is written to
  `<conf_file>_step_<index>` for offline analysis with `yang_mills_conf_analysis`.
- On restart (`start 2`) the data files are opened in append mode; otherwise they are overwritten.
- `walltime` (hours, default 24): the Monte Carlo loop stops when the elapsed time plus twice the longest step so far
  exceeds 99% of `walltime` (`wall_time_check` in [lib/timing.c](../lib/timing.c)). The final configuration is still
  saved (if `saveconf_back_every ≠ 0`).

Formats of all files: [output_files.md](output_files.md).
