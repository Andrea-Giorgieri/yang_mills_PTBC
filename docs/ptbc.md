# Parallel tempering on boundary conditions (PTBC)

PTBC is the algorithm that gives the repository its name. This page explains how it works, how to set it up and
where it is implemented.

## Background

Near the continuum limit, local updates change the topological charge Q very rarely ("topological freezing"), so
the autocorrelation time of Q grows very quickly. Open boundary conditions remove the barrier between topological
sectors but break translation invariance. Hasenbusch's parallel tempering on boundary conditions combines the two:

- N_r copies (replicas) of the lattice are simulated at the same β;
- they differ only by the boundary conditions on a small (D−1)-dimensional region of the boundary, the **defect**,
  interpolated by a coefficient c: c = 1 gives periodic boundary conditions, c = 0 open boundary conditions on the
  defect;
- after the updates, neighbouring replicas are swapped with a Metropolis test, so that a configuration can move from
  the open replica (fast topology) to the periodic one (physical) and back.

Only the replica with c = 1 is physical: measurements are taken on it.

## Action of a replica

Let μ = `defect_dir`. The defect D is the set of links in direction μ that cross the boundary of direction μ inside a
block of the boundary slice:

D = { (x, μ) : x_μ = L_μ − 1 and 0 ≤ x_ν < L^def_ν for every ν ≠ μ },

where L^def_ν are the values of `defect_size`. The links of the defect get a factor c, all the others 1. A plaquette p
gets the product K_p(c) = c^{n_p}, with n_p ∈ {0, 1, 2} the number of its links on the defect. The action of
replica a, with coefficient c_a, is

S_{c_a}[U] = β Σ_p K_p(c_a) (1 − (1/N) Re Tr[Z_p U_p])   (+ θ term, + multicanonical potential, if enabled).

The defect factors enter the staples of the heatbath and overrelaxation updates (`calcstaples_wilson_with_defect`,
`heatbath_with_defect`, `overrelaxation_with_defect` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c);
`link_defect_coeff`, `plaquette_defect_coeff` in [include/gauge_conf.h](../include/gauge_conf.h); the defect mask is
`is_on_defect` in [lib/geometry.c](../lib/geometry.c)).

PTBC can be combined with twisted and open boundary conditions, with the θ term and with the multicanonical algorithm.

## Input parameters

| Parameter | Meaning |
|---|---|
| `defect_dir μ` | direction orthogonal to the defect (default 0) |
| `defect_size L_1 … L_{D−1}` | sizes of the defect along the D−1 directions ν ≠ μ, in increasing order of ν. Each must not exceed the lattice size in that direction |
| `N_replica_pt n c_0 c_1 … c_{n−1}` | number of replicas and their coefficients. c_0 is the physical replica: c_0 = 1 (periodic) is the physically useful case; other values are accepted (useful for tests) but disable the translations (see "Translations"); usually c decreases to c_{n−1} = 0 |
| `hierarc_upd n_lev  R_0 … R_{n_lev−1}  S_0 … S_{n_lev−1}` | hierarchical update: number of levels, border R_i of the rectangle of level i, number of sweeps S_i of level i. `hierarc_upd 0` disables it |
| `swap_acc_file`, `swap_track_file` | output files, required when n > 1 |
| `num_gpus` | number of GPUs among which the replicas are distributed (0 = minimum needed); see [architecture_and_portability.md](architecture_and_portability.md) |

Example (template written by the program):

```
defect_dir    0             # space-time direction orthogonal to the topological defect
defect_size   2 2 2         # sizes of the defect along its dimensions (skipping defect_dir)
N_replica_pt  2    1.0 0.0  # number of parallel tempering replicas, followed by the boundary conditions coefficients
hierarc_upd   2    2 1    1 1
```

Programs: `yang_mills_local_pt_agf` (production) and `yang_mills_tuning_pt_mc` (multicanonical tuning). Both require
D = 4 and N ≥ 2 (checked at start).

## One parallel-tempering step

`parallel_tempering_with_hierarchical_update` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c):

1. **Full update** of all replicas (`update_with_defect`): one heatbath sweep and `overrelax` overrelaxation sweeps,
   each replica with its own coefficient.
2. Multicanonical accept/reject of each replica, if enabled ([multicanonical.md](multicanonical.md)).
3. **Swaps** between neighbouring replicas (if n > 1).
4. **Translation** of the physical replica (if n > 1 and c_0 = 1).
5. **Hierarchical update**, recursively from level 0: for level i, repeat S_i times
   - update all replicas only on the rectangle of level i (heatbath + `overrelax` overrelaxation sweeps),
   - multicanonical accept/reject on the rectangle, if enabled,
   - swaps and translation as in steps 3–4,
   - the hierarchical update of level i + 1.
6. The update index of every replica increases by 1: this whole procedure is **one** Monte Carlo step.

### Rectangles

The rectangle of level i contains the defect enlarged by R_i sites on each side in every direction: its size along
direction ν is L^def_ν + 2R_i (1 + 2R_i along μ), trimmed when it would exceed the lattice. Inside a rectangle, sites
are ordered even/odd and each parity is updated in parallel (`init_rect`, `rect_sizes`, `init_rect_utils` in
[lib/geometry.c](../lib/geometry.c)).

The trimming depends on how the rectangle is used:

- **Update rectangles** (and the cooling rectangles) are updated with the even/odd scheme. With OpenMP or OpenACC
  their size along each direction is limited to the largest even size not exceeding L_ν (L_ν − 1 along odd sides):
  along an odd side a rectangle of size L_ν would contain the neighbours x_ν = 0 and x_ν = L_ν − 1, which have the
  same parity and would be updated at the same time. In serial builds they are limited to L_ν.
- **Swap, clover and charge rectangles** are used only for sums over sites (swap action difference, clovers,
  change of the topological charge), so they are limited to L_ν in every build and always contain all the sites
  they need.

## Swaps

`swap` and `swap_group` in [lib/gauge_conf_paral_temp.c](../lib/gauge_conf_paral_temp.c).

- All n − 1 neighbouring pairs (a, a+1) are proposed, in two groups: even pairs (0,1), (2,3), … and odd pairs
  (1,2), (3,4), …; which group goes first is chosen at random, to preserve detailed balance.
- For the pair (a, b = a+1), with A the configuration that has index a and B the one that has index b, the action
  difference is
  ΔS = β [ (c_a − c_b)(P_1(A) − P_1(B)) + (c_a² − c_b²)(P_2(A) − P_2(B)) ]  (+ ΔV with the multicanonical algorithm),
  where P_n(X) is the sum of (1/N) Re Tr[Z_p U_p] over the plaquettes of X with n links on the defect. The swap is
  accepted with probability min(1, e^{−ΔS}).
- P_1 and P_2 depend only on the configuration, so they are computed once per configuration per swap call, with a
  reduction over the swap rectangle (border 1 around the defect), on the device holding the configuration
  (`swap_plaquette_sums` on the host; on the GPUs `swap_plaquette_sums_replica`, block sums per device added on the
  host in a fixed order, see [swap_translation_gpu.md](swap_translation_gpu.md)). The derivation is in
  [swap_p1p2.md](swap_p1p2.md).
- An accepted swap **exchanges the replica indices**, not the lattices: each configuration stays in its memory slot
  `GC[i]` for the whole run (with its twist factors, translation and multicanonical charge), and the map
  `replica_to_gc[a]` gives the slot of the configuration that currently has index a (`metropolis_single_swap`). On GPUs
  only scalars move between devices.
- The Metropolis tests run on the host, with the host random stream.

## Translations

After every swap call, if c_0 = 1 (within 10⁻¹³), the configuration that has index 0 is translated by one lattice
spacing in a random direction (positive verse), together with its factors Z (twist, open boundary, θ profile)
(`conf_translation` in [lib/gauge_conf_paral_temp.c](../lib/gauge_conf_paral_temp.c)). Since the defect is fixed,
this moves the defect with respect to the physical configuration. The cumulative translation is stored per
configuration, saved in the twist file and used to report slices and θ profiles in the original coordinates.

## Measurements and output

- Measurements are taken on the configuration that has index 0 (physical).
- With `--enable-meas-replica` (compile time) they are taken on every replica index a, and every output file gets the
  suffix `_replica_<a>`.
- `swap_acc_file`, written at the end of the run: for each pair, the two coefficients, the acceptance (%), its error,
  the numbers of accepted and proposed swaps. The counters start from zero at each run.
- `swap_track_file`, one line per step: the update index and, for each replica index a, the slot of the configuration
  that has index a (and its stored topological charge, with the multicanonical algorithm). It lets you follow each
  configuration through the replicas.
- Saved configurations: `<conf_file>_replica_<a>` and `<twist_file>_replica_<a>` contain the configuration that has
  index a; `start 2` reads them back in the same way. `saveconf_analysis_every` saves only the physical configuration.

Formats: [output_files.md](output_files.md).

## Parallelization

- All replicas are updated by the same loop: for each parity and direction, one parallel loop runs over
  (replica slots × sites), with OpenMP threads on CPU or one OpenACC kernel per GPU over the slots that GPU holds
  (`update_replica_sites`; `update_replica_border` for the border sites of odd sides).
- With several GPUs, each replica lives on one GPU for the whole run; the replicas are split into balanced contiguous
  blocks. Swaps need only P_1, P_2 and the multicanonical charges (scalars) on the host; after the swaps the new
  replica indices are copied to the GPUs (`copy_replica_index_to_device`).
- `yang_mills_gpu_plan input_file [GB_per_GPU]` estimates the memory per replica and the number of GPUs needed.

Details: [architecture_and_portability.md](architecture_and_portability.md).

## Files

| File | Content |
|---|---|
| [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c) | `parallel_tempering_with_hierarchical_update`, `update_with_defect`, `update_rectangle_with_defect`, `hierarchical_update_rectangle_with_defect`, staples with defect |
| [lib/gauge_conf_paral_temp.c](../lib/gauge_conf_paral_temp.c) | swaps, translations, acceptance and tracking files |
| [include/gauge_conf.h](../include/gauge_conf.h) | `Gauge_Conf` (fields `replica_index`, `device`, `translation`), `Acc_Utils` (`replica_to_gc`), defect factors |
| [lib/geometry.c](../lib/geometry.c) | defect mask, rectangles |
| [lib/gauge_conf_def.c](../lib/gauge_conf_def.c) | `init_gauge_conf_replica`, `write_replica_on_file` |
| [src/yang_mills_local_pt_agf.c](../src/yang_mills_local_pt_agf.c) | production main |
