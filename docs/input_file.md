# Input file reference

All simulation and analysis programs read the same kind of input file, parsed by `readinput` in
[lib/gparam.c](../lib/gparam.c). This page lists every keyword, its default, its constraints and which programs use it.
Which program to run and with which arguments: [programs.md](programs.md).

## Format

- Each entry is a keyword followed by its value(s), separated by spaces or newlines: `beta 6.4881`.
- `#` starts a comment that runs to the end of the line; empty lines are ignored.
- Keywords can appear in any order. If a keyword appears twice, the last value wins.
- An unknown keyword stops the program with an error (this catches typos).
- Every program accepts all keywords; keywords not used by a program are ignored.
- Running a simulation or analysis program without arguments writes `template_input.example` in the current
  directory, with the sections that the program uses.

Notation: L_μ = lattice size along μ, D = space-time dimension, N = number of colors, n = `N_replica_pt`.
"Required" means that the code sets no default and the value must be given.

## Lattice and action

| Keyword | Values | Default | Constraints / notes |
|---|---|---|---|
| `size` | L_0 … L_{D−1} | required | each > 1; even or odd. Direction 0 is "time" |
| `beta` | β | required | > 0 |
| `theta` | θ | 0 | used only with `--enable-use-theta` (D = 4) |
| `theta_profile_dir` | μ | −1 | −1 = constant θ; otherwise direction of the θ profile (θ builds only) |
| `theta_profile_file` | file name | – | required if `theta_profile_dir ≠ −1`; L_μ numbers |
| `htracedef` | h_1 … h_{⌊N/2⌋} | all 0 | trace deformation couplings (`yang_mills_tracedef`); with OpenMP and any h ≠ 0, spatial sizes must be even |
| `k_twist` | D(D−1)/2 integers | all 0 | twist per plane, order (0,1), (0,2), …, (0,D−1), (1,2), …, (D−2,D−1) |
| `obc_dir` | μ | −1 | −1 = no open boundary; otherwise direction of the open boundary |

Details: [physics_and_parameters.md](physics_and_parameters.md).

## Parallel tempering on boundary conditions

Used by `yang_mills_local_pt_agf`, `yang_mills_tuning_pt_mc`, `yang_mills_gpu_plan`. Details: [ptbc.md](ptbc.md).

| Keyword | Values | Default | Constraints / notes |
|---|---|---|---|
| `defect_dir` | μ | 0 | 0 ≤ μ < D |
| `defect_size` | D−1 integers | all 0 | sizes along the directions ≠ μ in increasing order; each ≤ the lattice size in that direction |
| `N_replica_pt` | n c_0 … c_{n−1} | 1, c_0 = 1 | n > 0; c_0 = 1 is the physically useful case (periodic physical replica); other values are accepted, e.g. for tests, but then the translations are skipped ([ptbc.md](ptbc.md)) |
| `hierarc_upd` | n_lev R_0 … R_{n_lev−1} S_0 … S_{n_lev−1} | 0 | hierarchical levels, borders of the rectangles, sweeps per level; all ≥ 0 |
| `num_gpus` | k | 0 | 0 = smallest number of GPUs that can hold the replicas; k > 0 forces k GPUs (k ≤ `Num_devices`, k ≤ n). On CPU, k > 1 simulates the splitting |

## Monte Carlo run

| Keyword | Values | Default | Constraints / notes |
|---|---|---|---|
| `sample` | integer | required | number of Monte Carlo steps. 0 = only measure the initial configuration |
| `thermal` | integer | required | no measurements before this update index |
| `overrelax` | integer | required | overrelaxation sweeps per heatbath sweep |
| `measevery` | integer | required | measure every this many steps; use ≥ 1 |
| `start` | 0, 1, 2, 3 | required | 0 ordered, 1 random, 2 from file, 3 ordered with twist (see [algorithms.md](algorithms.md#starting-configuration-checkpoints-restart)) |
| `saveconf_back_every` | integer | required | save configuration(s) and a rotating backup every this many steps; 0 = never |
| `saveconf_analysis_every` | integer | required | save numbered configurations for `yang_mills_conf_analysis`; 0 = never |
| `randseed` | unsigned integer | required | 0 = seed from the current time |
| `walltime` | hours | 24 | the run stops safely before this time |
| `epsilon_metro` | ε ≥ 0 | 0.25 | Metropolis step size (`yang_mills_tracedef`) |

The step counters (`sample`, `thermal`, …) count Monte Carlo steps as defined in
[algorithms.md](algorithms.md#monte-carlo-step-of-each-program).

## Measurements

All flags are 0/1. Details and outputs: [physics_and_parameters.md](physics_and_parameters.md#observables).

| Keyword | Default | Needs |
|---|---|---|
| `plaquette_meas` | 1 | – |
| `clover_energy_meas` | 0 | – |
| `charge_meas` | 0 | 4D SU(N) or 2D U(1) |
| `polyakov_meas` | 0 | – |
| `polyakov_powers_meas` | 0 | – |
| `multipolyakov_order` | 0 | `multipolyakov_order n μ_1 … μ_n`, each 0 ≤ μ_i < D |
| `charge_prime_meas` | 0 | 4D SU(N) or 2D U(1) |
| `chi_prime_meas` | 0 | 4D SU(N) or 2D U(1); `chiprime_data_file` |
| `energy_density_meas` | 0 | `energy_density_file` |
| `charge_density_meas` | 0 | 4D SU(N) or 2D U(1); `charge_density_file` |
| `polyakov_density_meas` | 0 | `polyakov_density_file` |
| `energy_slices_meas` | 0 | `energy_slices_file` |
| `charge_slices_meas` | 0 | 4D SU(N) or 2D U(1); `charge_slices_file` |
| `charge_p_slices_meas` | 0 | 4D SU(N) or 2D U(1); `charge_slices_file`; direction from `test_flag` |
| `action_meas` | 0 | – |
| `meas_effective_charge` | 0 | θ builds with a θ profile |
| `test_flag` | 0 | integer; used as the second direction by `charge_p_slices_meas` |

## Smoothing

| Keyword | Default | Constraints / notes |
|---|---|---|
| `agf_length` | 0 | total flow time of the adaptive flow; the adaptive flow is active if `agf_length/agf_meas_each ≥ 1` |
| `agf_meas_each` | 1 | flow-time interval between measurements; must be > `agf_step` and > `agf_time_bin` |
| `agf_step` | 0.1 | initial step, > 0 |
| `agf_delta` | 1e−5 | tolerance on the links for the step control, > 0 |
| `agf_time_bin` | 0 | tolerance on the flow time at which measurements are taken |
| `num_gfsteps` | 0 | number of steps of the fixed-step flow; active if `num_gfsteps/gf_meas_each ≥ 1` |
| `gfstep` | 1 | step of the fixed-step flow |
| `gf_meas_each` | 1 | measure every this many flow steps, > 0 |
| `coolrepeat` | 0 | number of cooling measurements; 0 = no cooling |
| `coolsteps` | 1 | cooling sweeps between measurements |
| `cooling_type` | `LEX_DIR_LEXEO_SITE` | one of `LEX_DIR_LEXEO_SITE`, `LEXEO_SITE_LEX_DIR`, `RND_DIR_RNDEO_SITE`, `LEX_SITE_LEX_DIR`, `LEX_DIR_LEX_SITE` |

Details: [physics_and_parameters.md](physics_and_parameters.md#smoothing-gradient-flow-and-cooling).

## Multicanonical algorithm

Used only with `--enable-use-multicanonical`. Details: [multicanonical.md](multicanonical.md).

| Keyword | Default | Constraints / notes |
|---|---|---|
| `topo_potential_file` | – | required; potential file (format in [multicanonical.md](multicanonical.md#potential)) |
| `multicanonic_acc_file` | – | required; acceptance output |
| `grid_step` | 0.25 | > 0; must match the file |
| `grid_max` | 2.0 | > 0; must match the file |
| `topo_cooling` | 0 | 0 = bare charge, 1 = charge after cooling |
| `topo_coolsteps` | 0 | cooling sweeps for `topo_cooling 1` |
| `topo_alpha` | 0 | > 0: α-rounding of the charge; 0 = none |
| `topo_tuning_stp` | 0.1 | tuning program: initial step, > 0 |
| `topo_tuning_thr` | 0.1 | tuning program: stopping threshold |
| `topo_tuning_weights_every` | 1000 | tuning program: proposed swaps between weight updates |
| `topo_tuning_even` | 1 | tuning program: keep the potential even |
| `topo_tuning_save_every` | 0 | tuning program: save the potential every this many steps (0 = never) |

## Multilevel

Used by `yang_mills_multilevel` and `yang_mills_multilevel_long`. Details: [algorithms.md](algorithms.md#multilevel-algorithm).

| Keyword | Default | Constraints / notes |
|---|---|---|
| `ml_obs` | `NONE` | `NONE`, `POLYCORR`, `POLYCORR_LONG`, `TUBE_DISC`, `TUBE_CONN`, `TUBE_CONN_LONG`. It also selects the header of the data file: leave it `NONE` for the other programs |
| `ml_step` | `Num_levels` integers, `Num_levels − i` for level i | with `Num_levels > 1`: L_0 divisible by `ml_step[0]`, each divides the previous and is smaller, last > 1 |
| `ml_upd` | `Num_levels` integers, 0 | updates per level |
| `multihit` | 0 | multihit steps |
| `ml_level0_repeat` | 0 | `yang_mills_multilevel_long`: runs per level-0 measurement |
| `ml_file` | – | required with `Num_levels > 1` |
| `dist_poly` | 0 | distance between the Polyakov loops (direction 1) |
| `transv_dist` | 0 | transverse distance of the plaquette (direction 2) |
| `plaq_dir` | 0 1 | orientation of the plaquette for flux tubes |

## Files

| Keyword | Required | Content |
|---|---|---|
| `conf_file` | always | configuration (base name; suffixes are added for replicas, backups, numbered configurations) |
| `twist_file` | always | positions of twist and open boundary, translation (base name, same suffixes) |
| `data_file` | always | measurements |
| `log_file` | always | run summary written at the end |
| `swap_acc_file` | if n > 1 | swap acceptances |
| `swap_track_file` | if n > 1 | replica tracking |
| `energy_density_file`, `charge_density_file`, `polyakov_density_file` | if the corresponding flag is 1 | binary densities |
| `chiprime_data_file` | if `chi_prime_meas 1` | χ′ |
| `energy_slices_file`, `charge_slices_file` | if the corresponding flag is 1 | slices |
| `topo_potential_file`, `multicanonic_acc_file` | multicanonical builds | potential and acceptances |
| `theta_profile_file` | if `theta_profile_dir ≠ −1` | θ profile (input) |
| `ml_file` | if `Num_levels > 1` | multilevel state |

File names are limited to 149 characters (`STD_STRING_LENGTH` in [include/macro.h](../include/macro.h)), including the
suffixes added by the code. Formats: [output_files.md](output_files.md).

## Example: PTBC with twist and adaptive flow

```
size 24 16 16 24

defect_dir    0
defect_size   4 4 4
N_replica_pt  4   1.0 0.66 0.33 0.0
hierarc_upd   2   2 1   1 1

k_twist 0 0 0 1 0 0

beta 6.459
sample 100000
thermal 0
overrelax 12
measevery 10
start 3
saveconf_back_every 100
saveconf_analysis_every 0
randseed 0
walltime 24

plaquette_meas 1
clover_energy_meas 1
charge_meas 1

agf_length 20
agf_step 0.01
agf_meas_each 1
agf_delta 0.00001
agf_time_bin 0

conf_file conf.dat
twist_file twist.dat
data_file dati.dat
log_file log.dat
swap_acc_file swap_acc.dat
swap_track_file swap_track.dat
```
