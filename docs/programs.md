# Programs and how to run them

This page describes each executable: what it does, how to call it, which input sections it uses and what it writes.
The input keywords are listed in [input_file.md](input_file.md), the output formats in [output_files.md](output_files.md),
and the build options in [building.md](building.md).

## Common behaviour

- Simulation and analysis programs take an input file as first argument. Without arguments they print a short
  description (purpose, action, algorithms, measures), the usage and the compile-time configuration, and write
  `template_input.example` in the current directory. The utilities `conf_check` and `conf_txt_to_bin` print only the
  usage and the compile-time configuration; `yang_mills_gpu_plan` prints only the usage.
- At the end of the run they write `log_file` with the compile-time configuration, the parameters and the timers
  (initialization, update, measurement and step times).
- The run stops before `walltime` (see [algorithms.md](algorithms.md#starting-configuration-checkpoints-restart)).
- Restart: set `start 2`; the configuration(s) and twist file(s) are read back and the data files are appended to.

## Overview

| Program | Purpose | Build | GPU |
|---|---|---|---|
| `yang_mills_local_pt_agf` | PTBC production (optionally multicanonical, θ) | any; 4D SU(N) at run time | yes |
| `yang_mills_tuning_pt_mc` | PTBC + tuning of the multicanonical potential | `--enable-use-multicanonical`; 4D SU(N) | yes |
| `yang_mills_local_agf` | single configuration (optionally multicanonical, θ) | any | yes |
| `yang_mills_conf_analysis` | measurements on saved configurations | any | yes |
| `yang_mills_gpu_plan` | GPU memory planner for PTBC inputs | any | no GPU needed |
| `yang_mills_multilevel` | multilevel Polyakov correlators and flux tubes | not OpenACC, not multicanonical | no |
| `yang_mills_multilevel_long` | the same, split over several executions | not OpenACC, not multicanonical | no |
| `yang_mills_tracedef` | trace-deformed Yang–Mills | not OpenACC | no |
| `conf_check` | check the MD5 hash of a configuration file | any | – |
| `conf_txt_to_bin` | convert a text configuration to the binary format | N ≥ 3 | – |

## `yang_mills_local_pt_agf`

Parallel tempering on boundary conditions ([ptbc.md](ptbc.md)), with measurements on the physical replica, including
adaptive and fixed-step gradient flow and cooling.

```
yang_mills_local_pt_agf input_file
```

- Requires D = 4 and N ≥ 2 (checked at start).
- Input sections: lattice, PTBC (`defect_dir`, `defect_size`, `N_replica_pt`, `hierarc_upd`, `num_gpus`), twist / open
  boundary, multicanonical (multicanonical builds), Monte Carlo run, measurements, smoothing, files
  (`swap_acc_file`, `swap_track_file` required when `N_replica_pt > 1`).
- Each step: one PTBC step (full update, swaps, translation, hierarchical updates). Measurements on the configuration
  with replica index 0, every `measevery` steps from `thermal` on (on all replicas with `--enable-meas-replica`).
- Output: data file(s) and the other measurement files, `swap_track_file` (every step), `swap_acc_file` (end),
  `multicanonic_acc_file` (multicanonical builds), configurations `<conf_file>_replica_<a>` (and backups),
  `<conf_file>_step_<index>` for the physical replica if `saveconf_analysis_every > 0`, log file.
- On GPUs, the replicas are distributed among up to `Num_devices` GPUs (see
  [architecture_and_portability.md](architecture_and_portability.md)); the distribution is printed in the log.

Main: [src/yang_mills_local_pt_agf.c](../src/yang_mills_local_pt_agf.c).

## `yang_mills_tuning_pt_mc`

Runs PTBC like `yang_mills_local_pt_agf` and tunes the multicanonical potentials of all replicas
([multicanonical.md](multicanonical.md#tuning-the-potential-yang_mills_tuning_pt_mc)).

```
yang_mills_tuning_pt_mc input_file
```

- Requires `--enable-use-multicanonical`, D = 4, N ≥ 2.
- Input: as `yang_mills_local_pt_agf`, plus the `topo_tuning_*` keywords. `topo_potential_file` must exist (initial
  potential and weights).
- The tuning starts after `thermal` steps and ends when the convergence threshold `topo_tuning_thr` is reached, after
  `sample` steps or at the wall time. The new tuning steps are printed on standard output.
- Output: as `yang_mills_local_pt_agf`, plus `<topo_potential_file>_step_<index>` every `topo_tuning_save_every` steps
  and the final potential written to `topo_potential_file` (the input file is overwritten). The log reports the number
  of steps.

Main: [src/yang_mills_tuning_pt_mc.c](../src/yang_mills_tuning_pt_mc.c).

## `yang_mills_local_agf`

Standard simulation of one configuration: heatbath + overrelaxation, with all boundary conditions, the θ term
(θ builds) and the multicanonical algorithm (multicanonical builds), and the same measurements as the PTBC program.

```
yang_mills_local_agf input_file
```

- Any D and N (with the limits of θ and topological observables, see
  [physics_and_parameters.md](physics_and_parameters.md#summary-of-physical-limits)).
- Input sections: lattice, twist / open boundary, θ profile, multicanonical, Monte Carlo run, measurements, smoothing,
  files. PTBC keywords are not used by the simulation; leave `N_replica_pt` at its default 1, because the input checks
  (required swap files, number of columns of the potential file) depend on it. With the multicanonical algorithm the
  potential file then has one column of V.
- Output: data file and the other measurement files, `conf_file` / `twist_file` and backups, numbered configurations,
  `multicanonic_acc_file`, log file.
- Typical use: θ term with open boundary conditions and a θ profile (production inputs in
  [`prod_inputs/theta_obc`](../prod_inputs/theta_obc)).

Main: [src/yang_mills_local_agf.c](../src/yang_mills_local_agf.c).

## `yang_mills_conf_analysis`

Measures observables on configurations saved with `saveconf_analysis_every`.

```
yang_mills_conf_analysis input_file start_index stop_index
```

- Reads, for every update index i with `start_index ≤ i < stop_index`, the files `<conf_file>_step_<i>` and
  `<twist_file>_step_<i>`; indices without both files are skipped. It stops with an error if no configuration is found.
- Use the input file of the run that generated the configurations, changing the measurement and smoothing keywords as
  needed. PTBC keywords are not used for the measurements; with a PTBC input the saved configurations are those of the
  physical replica, and in multicanonical builds the reweighting factor uses the potential of replica 0.
- To avoid overwriting the files of the original run, the code appends a suffix to the output file names: `_agf`,
  `_gf`, `_cool_<cooling_type>` for each active smoothing method, or `_hot` if none is active. The output files are
  appended to if they exist (the program sets `start 2` internally).
- With the multicanonical algorithm, the charge used by the potential is recomputed for each configuration, so that
  the data file contains the reweighting factor.
- Output: measurement files (with suffix) and log file.

Main: [src/yang_mills_conf_analysis.c](../src/yang_mills_conf_analysis.c).

## `yang_mills_gpu_plan`

Estimates the memory needed by the PTBC replicas and how many GPUs are needed, without using any GPU.

```
yang_mills_gpu_plan input_file [memory_per_gpu_GB]
```

- The input file is the one of `yang_mills_local_pt_agf` or `yang_mills_tuning_pt_mc`; the memory per GPU defaults to
  64 GB. The estimate depends on the build options (θ term, multicanonical, Philox), so use the same build as the
  production program.
- Prints the memory per replica, the memory shared by every GPU, the maximum number of replicas per GPU, the
  distribution of the replicas and the number of GPUs needed (at most `Num_devices`), or that they do not fit. Honors
  `num_gpus`.

Main: [src/yang_mills_gpu_plan.c](../src/yang_mills_gpu_plan.c).

## `yang_mills_multilevel`

Simulation with heatbath + overrelaxation and multilevel measurements ([algorithms.md](algorithms.md#multilevel-algorithm)).

```
yang_mills_multilevel input_file
```

- CPU only; refuses multicanonical builds; all spatial sizes must be equal.
- Input sections: lattice, twist, Monte Carlo run, multilevel (`ml_obs`, `ml_step`, `ml_upd`, `multihit`,
  `dist_poly`, `transv_dist`, `plaq_dir`, `ml_file`), files.
- Every `measevery` steps (from `thermal` on) it runs the multilevel update and writes one line in the data file
  (columns in [output_files.md](output_files.md#data-file-data_file)). The local observables of
  [physics_and_parameters.md](physics_and_parameters.md#observables) are not measured.
- Compile-time variants `OPT_MULTIHIT` and `OPT_MULTILEVEL` print diagnostics to choose `multihit` and `ml_upd`
  instead of the measurement ([building.md](building.md)).

Main: [src/yang_mills_multilevel.c](../src/yang_mills_multilevel.c).

## `yang_mills_multilevel_long`

Multilevel measurement split over several executions, for measurements that do not fit in one job.

```
yang_mills_multilevel_long input_file
```

- Supports `ml_obs` `POLYCORR`, `TUBE_CONN` (and their `_LONG` forms) or `NONE`; same limits as `yang_mills_multilevel`.
- Each execution performs one stage and exits:
  - if no multilevel is in progress, `measevery` standard update steps;
  - otherwise one level-0 multilevel iteration; after `ml_level0_repeat` iterations the measurement is written to the
    data file and the next execution starts again with standard updates.
- The state of the multilevel (iteration and partial averages) is written to `ml_file` at the end of every execution
  and read back when `start 2`. The configuration is always saved at the end. The program is meant to be run
  repeatedly with `start 2`.

Main: [src/yang_mills_multilevel_long.c](../src/yang_mills_multilevel_long.c).

## `yang_mills_tracedef`

Yang–Mills with a trace deformation of the Polyakov loop ([physics_and_parameters.md](physics_and_parameters.md#trace-deformation)).

```
yang_mills_tracedef input_file
```

- CPU only. With OpenMP and any `htracedef` coupling ≠ 0 the spatial sizes must be even.
- Input: as `yang_mills_local_agf`, plus `htracedef` and `epsilon_metro`.
- Each step: heatbath of the spatial links, Metropolis of the temporal links, `overrelax` overrelaxation sweeps of the
  spatial links (and the multicanonical test in multicanonical builds).
- Output: as `yang_mills_local_agf`; the log also reports the mean Metropolis acceptance.

Main: [src/yang_mills_tracedef.c](../src/yang_mills_tracedef.c).

## `conf_check`

```
conf_check conf_file
```

Recomputes the MD5 hash of the links of a binary configuration file and compares it with the hash stored in its
header; prints `Check passed` or stops with an error. Works only if the code was built with OpenSSL (`HASH_MODE`);
otherwise it does nothing. Built for a given `N_c` (the link format depends on N). Source:
[src/conf_check.c](../src/conf_check.c).

## `conf_txt_to_bin`

```
conf_txt_to_bin conf_file_in conf_file_out twist_file_out twist_mu twist_nu twist_k size_0 size_1 ... size_{D-1}
```

Converts a configuration from text to the binary format of the code and writes a twist file with twist k on the plane
(μ, ν). The text file has one line per matrix element:

```
<ignored> <ignored> x_0 ... x_{D-1} mu i j (re,im)
```

with colour indices i, j from 1 to N; exactly D·N²·V lines are expected. Only for N ≥ 3. The output has update index 0.
Source: [src/conf_txt_to_bin.c](../src/conf_txt_to_bin.c).

## Test programs

`debug_rng`, `debug_su2`, `debug_sun`, `debug_u1`, `debug_twist`, `debug_agf_vs_gf`, `debug_agf_vs_delta`,
`debug_device`, `debug_device_data`, `debug_acc_routines` (OpenACC builds) are tests used during development (random
numbers, group algebra, twist, gradient-flow integrators, device data and kernels). They are not described here; see
their sources in [src/](../src).
