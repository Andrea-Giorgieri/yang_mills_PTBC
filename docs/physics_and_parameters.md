# Physics: what can be simulated

This page describes the lattice theories that the code can simulate, the boundary conditions, the extra terms in the
action and the observables, with the parameters that select them. For the algorithms see [algorithms.md](algorithms.md);
for the full list of input parameters see [input_file.md](input_file.md).

Convention: **compile time** = option of `configure`/`build.sh` ([building.md](building.md)); **runtime** = keyword in
the input file.

## Gauge group and dimension

| Feature | Selection | Limits |
| --- | --- | --- |
| U(1) | compile time, `N_c=1` | – |
| SU(2) | compile time, `N_c=2` | dedicated quaternion-like representation ([include/su2.h](../include/su2.h)) |
| SU(N), N ≥ 3 | compile time, `N_c=N` | generic N×N complex matrices ([include/sun.h](../include/sun.h)) |
| Space-time dimension D | compile time, `ST_dim=D` | D ≥ 2; some features are 4D only (see below) |
| Lattice sizes L_0 … L_{D-1} | runtime, `size L0 L1 ...` | every L_μ > 1. Even and **odd** sizes are allowed (see [architecture_and_portability.md](architecture_and_portability.md)). Direction 0 is "time" |

The group is selected in [include/gauge_group.h](../include/gauge_group.h), which maps every group operation to the
U(1), SU(2) or SU(N) implementation.

## Action

### Wilson action

The basic action is the Wilson plaquette action

S_W = β Σ_p K_p (1 − (1/N) Re Tr[Z_p U_p]),

where U_p is the ordered product of links around the plaquette p, Z_p is a complex factor that implements twisted
and open boundary conditions (and is 1 otherwise), and K_p is the defect factor used by parallel tempering on boundary
conditions (K_p = 1 outside the PTBC defect, see [ptbc.md](ptbc.md)).

- β: runtime, `beta` (must be positive).
- The factors Z are stored per site and per oriented plane in `GC->Z[r][plane]` ([include/gauge_conf.h](../include/gauge_conf.h))
  and are built by `initialize_Z_with_copy` in [lib/gauge_conf_def.c](../lib/gauge_conf_def.c).
- Staples: `calcstaples_wilson` (and the variants with defect and θ term) in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c).

### Imaginary θ term

With `--enable-use-theta` (compile time) the action becomes

S = S_W − θ Q,

where `theta` (runtime) is the coefficient of the topological charge (an imaginary θ in Minkowskian conventions, so
the Boltzmann weight stays real) and Q is the clover discretization of the topological charge described below. The
term is included exactly in the link updates through "topological staples" built from clovers
(`calcstaples_with_topo` and `compute_clovers` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c)); the clovers of
the planes orthogonal to the updated direction are recomputed before each sweep of a direction.

Limits: D = 4 only (compile-time `#error` and runtime check in `readinput`), N ≥ 2.

### Position-dependent θ profile

With the θ term enabled, θ can depend on one coordinate: θ(x) = θ · f(x_μ).

| Parameter | Meaning |
| --- | --- |
| `theta_profile_dir` | direction μ of the profile; `-1` (default) = constant θ |
| `theta_profile_file` | text file with exactly L_μ numbers f(0) … f(L_μ−1) (whitespace separated; more or fewer values are an error) |
| `meas_effective_charge` | if 1, the measured topological charge density is multiplied by f(x_μ), i.e. the code measures Q_eff = Σ_x f(x_μ) q(x) |

The profile is stored as the real factor `Z[r][n_planes]` of every site and multiplies the local charge density in the
action. Reading: `read_theta_profile` in [lib/gparam.c](../lib/gparam.c); use in the staples: `calcstaples_with_topo`.
When a configuration is translated (parallel tempering, see [ptbc.md](ptbc.md)), the profile is translated with it and
the translation is saved in the twist file.

### Trace deformation

The program `yang_mills_tracedef` adds the term

S_TD = Σ_{x spatial} Σ_{i=1}^{⌊N/2⌋} h_i |Tr P(x)^i|²,

where P(x) is the Polyakov loop in direction 0 through the spatial site x and Tr is the full trace (not divided by N).

| Parameter | Meaning |
| --- | --- |
| `htracedef h_1 ... h_{⌊N/2⌋}` | couplings, exactly ⌊N/2⌋ values (runtime) |

Implementation: `update_with_trace_def`, `metropolis_with_tracedef`, `calcstaples_tracedef` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c).

Limits: CPU only (not built with OpenACC). With OpenMP and any h_i ≠ 0, all spatial sizes must be even.

## Boundary conditions

| Type | Selection | Description |
| --- | --- | --- |
| Periodic | default | – |
| Twisted ('t Hooft) | runtime, `k_twist` | Z_N twist with integer k_{μν} on each plane, applied as the phase exp(2πi k_{μν}/N) on one plaquette per plane |
| Open | runtime, `obc_dir` | open boundary along one direction |
| Defect (PTBC) | runtime, `defect_dir`, `defect_size`, `N_replica_pt` | boundary conditions interpolating between periodic and open on a region of the boundary, see [ptbc.md](ptbc.md) |

All of them act only through the plaquette factors Z_p and K_p, so they can be combined.

### Twisted boundary conditions

`k_twist` takes D(D−1)/2 integers, one per plane, in the order (0,1), (0,2), …, (0,D−1), (1,2), …, (D−2,D−1).
The phase multiplies the plaquette of the plane (μ,ν) with corner at x_μ = x_ν = 0 (or at the position read from the
twist file when restarting); the opposite orientation gets the complex conjugate.

- `start 3` initializes an ordered configuration compatible with the twist (twist-eater matrices P and Q on the
  boundary links). It handles **one** twisted plane: if several k_{μν} are non-zero, the last one is used.
- The twist position is saved in the twist file (see [output_files.md](output_files.md)).

Files: `initialize_Z_with_copy`, `init_gauge_conf_from_file_with_name` (start options) in
[lib/gauge_conf_def.c](../lib/gauge_conf_def.c).

### Open boundary conditions

`obc_dir μ` (default `-1` = periodic) cuts the lattice between the slices x_μ = x_obc and x_μ = x_obc + 1, with
x_obc = L_μ − 1 by default (or the position read from the twist file when restarting):

- plaquettes in planes containing μ that cross the cut get Z = 0;
- plaquettes in planes not containing μ that lie on the two boundary slices get Z = 1/2.

The open boundary moves with the configuration when it is translated (if using parallel tempering) and its position is
saved in the twist file.

## Observables

All measurements are switched on at runtime with 0/1 flags and are performed every `measevery` Monte Carlo steps
after `thermal` steps. They are computed on the configuration as it is ("hot" measurements) and, optionally, along the
gradient flow and during cooling (next section). Implementation: [lib/gauge_conf_meas.c](../lib/gauge_conf_meas.c),
function `perform_measures_localobs`.

Normalized trace: tr = (1/N) Re Tr below.

| Flag | Observable | Output |
| --- | --- | --- |
| `plaquette_meas` | mean spatial and temporal plaquette tr(Z_p U_p) (default **on**) | data file |
| `clover_energy_meas` | clover discretization of the energy density Σ_{μν} Tr(F_{μν}F_{μν})/2, averaged over the volume (`loc_clover_energy`) | data file |
| `charge_meas` | topological charge Q (clover discretization, see below) | data file |
| `polyakov_meas` | volume average of tr P_μ and its imaginary part, for every direction μ | data file |
| `polyakov_powers_meas` | tr P_0^k for k = 1 … ⌊N/2⌋+1 | data file (hot measurement only) |
| `multipolyakov_order n μ_1 … μ_n` | trace of the loop winding successively along μ_1, …, μ_n, averaged over the volume | data file |
| `charge_prime_meas` | Q′_μ = Σ_{x_μ} \ | Σ_{x_⊥} q(x)\ | for every direction μ | data file |
| `chi_prime_meas` | estimator (q(0)/(2D)) Σ_x d(x,0)² q(x) of χ′ (d = periodic lattice distance) | `chiprime_data_file` |
| `energy_density_meas` | local clover energy density on every site | binary `energy_density_file` |
| `charge_density_meas` | local topological charge density on every site | binary `charge_density_file` |
| `polyakov_density_meas` | local Polyakov loop (Re and Im of tr) for every direction | binary `polyakov_density_file_dir<μ>` |
| `energy_slices_meas` | clover energy averaged on each slice x_0 = const | `energy_slices_file` |
| `charge_slices_meas` | topological charge summed on each slice x_0 = const | `charge_slices_file` |
| `charge_p_slices_meas` | Fourier components in direction ν = `test_flag` of the charge on each slice x_0 (feature in development) | `charge_slices_file` |
| `action_meas` | Wilson and θ parts of the action, total action and multicanonical potential (diagnostic) | data file |

Slices take into account the translation of the configuration (they are reported in the original coordinates).

### Topological charge

In the continuum, for a general gauge group in 4 dimensions,

Q = 1/(64π²) ∫ d⁴x ε_{μνρσ} F^a_{μν} F^a_{ρσ},

with the generators normalized as Tr(T_a T_b) = K δ_{ab} such that the longest root of the representation has length 1
(C. W. Bernard, N. H. Christ, A. H. Guth, E. J. Weinberg, Phys. Rev. D 16 (1977) 2967). For SU(N) this basis is the
generalized Gell-Mann matrices divided by two, Tr(T_a T_b) = δ_{ab}/2, so that F^a_{μν} F^a_{ρσ} = 2 Tr(F_{μν} F_{ρσ}) and

Q = 1/(32π²) ∫ d⁴x ε_{μνρσ} Tr(F_{μν} F_{ρσ}).

On the lattice F_{μν} is taken from the clover C_{μν} (sum of the four plaquettes around x in the plane (μ,ν)), which
satisfies C_{μν} ≃ 4 + 4i a² F_{μν}, so that

F_{μν} = (C_{μν} − C_{μν}†)/(8i),   Tr(F_{μν} F_{ρσ}) = −(1/2⁵) Re Tr(C_{μν} [C_{ρσ} − C_{ρσ}†]).

Of the 24 terms of ε_{μνρσ} Tr(F_{μν} F_{ρσ}) only 3 are independent: each independent term appears 2³ = 8 times
(exchange μ↔ν, exchange ρ↔σ, exchange of the two pairs), and the 3 independent permutations are e.g. 0123, 0213, 0312.
Hence the local charge density used by the code is

q(x) = −1/(128π²) Σ_{independent (μνρσ)} ε_{μνρσ} Re Tr(C_{μν}(x) [C_{ρσ}(x) − C_{ρσ}(x)†]),

and Q = Σ_x q(x). Implementation: `loc_topcharge` in [lib/gauge_conf_meas.c](../lib/gauge_conf_meas.c), with the
independent permutations in `g_indep_perm_dir` ([lib/geometry.c](../lib/geometry.c)). With a θ profile and
`meas_effective_charge 1`, q(x) is multiplied by the profile f(x_μ). This derivation is taken from the README of the original
`yang-mills` package ("A note on the topological charge"), where the factor 1/π² is left implicit in the last step.

> **To be verified.** `loc_topcharge` computes (N/(128π²)) Σ_i s_i tr(C_{ρσ}(C_{μν}† − C_{μν})) with
> tr = (1/N) Re Tr, permutations (0123), (0213), (0312) and signs s = (−1, +1, −1), i.e. s_i = −ε_{μνρσ}. Written in
> the form above this is +1/(128π²) Σ ε Re Tr(C_{μν}[C_{ρσ} − C_{ρσ}†]), with the opposite overall sign. The two agree
> if the clover built by `clover` satisfies C_{μν} ≃ 4 − 4i a² F_{μν} (i.e. links U = exp(−i a A)). The convention
> has not been checked here.

Availability: topological observables (`charge_meas`, `charge_density_meas`, `charge_slices_meas`,
`charge_p_slices_meas`, `charge_prime_meas`, `chi_prime_meas`) are allowed only for D = 4 with N ≥ 2, and for D = 2
with U(1); there the charge is Σ_x arg(U_p(x))/(2π). Otherwise `readinput` stops with an error.

## Smoothing: gradient flow and cooling

Smoothing is used only for measurements: the configuration is smoothed, the selected observables are measured at
several smoothing times, then the original configuration is restored. The observables measured during smoothing are
the ones of the table above except `polyakov_powers_meas`; densities and slices are written with the index of the
smoothing step.

| Method | Parameters | Algorithm |
| --- | --- | --- |
| Wilson flow, fixed step | `gfstep` (ε), `num_gfsteps`, `gf_meas_each` | third-order Runge–Kutta of Lüscher (arXiv:1006.4518, app. C); measurements every `gf_meas_each` steps |
| Wilson flow, adaptive step | `agf_length` (total flow time), `agf_step` (initial step), `agf_meas_each` (flow-time interval between measurements), `agf_delta` (tolerance), `agf_time_bin` | adaptive step-size control of Fritzsch and Ramos (arXiv:1301.4388, app. D): the step is accepted if the maximum link distance between the 2nd- and 3rd-order solutions is below `agf_delta`; the step is shortened to land on the measurement times |
| Cooling | `cooling_type`, `coolsteps` (sweeps between measurements), `coolrepeat` (number of measurements) | local minimization of the Wilson action (with twist/open factors) |

Cooling orderings (`cooling_type`):

| Value | Order of the sweep |
| --- | --- |
| `LEX_DIR_LEXEO_SITE` (default) | directions in order; for each, even sites then odd sites (parallel) |
| `LEXEO_SITE_LEX_DIR` | even sites for all directions, then odd sites (parallel) |
| `RND_DIR_RNDEO_SITE` | random order of directions and random parity first (parallel); does not break isotropy and parity |
| `LEX_SITE_LEX_DIR` | sequential: sites in lexicographic order, all directions per site (host only) |
| `LEX_DIR_LEX_SITE` | sequential: directions in order, sites in lexicographic order (host only) |

A method is active when its number of measurements is positive: `agf_length / agf_meas_each ≥ 1`,
`num_gfsteps / gf_meas_each ≥ 1`, `coolrepeat ≥ 1`. Several methods can be active at once; they are applied in the
order adaptive flow, fixed-step flow, cooling, each starting from the unsmoothed configuration.

Checks in `readinput`: `agf_meas_each > agf_time_bin`, `agf_meas_each > agf_step`.

Files: `gradflow_RKstep`, `gradflow_RKstep_adaptive`, `cooling` in [lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c);
`perform_measures_localobs_*` in [lib/gauge_conf_meas.c](../lib/gauge_conf_meas.c).

## Multilevel observables

The programs `yang_mills_multilevel` and `yang_mills_multilevel_long` measure, with the multilevel algorithm
(see [algorithms.md](algorithms.md)):

| `ml_obs` | Observable |
| --- | --- |
| `POLYCORR` | correlator of two Polyakov loops at distance `dist_poly` along direction 1 |
| `TUBE_DISC` | Polyakov loop correlator and its disconnected correlation with a plaquette (flux tube) |
| `TUBE_CONN` | as `TUBE_DISC` plus the connected correlator |
| `POLYCORR_LONG`, `TUBE_CONN_LONG` | the same as `POLYCORR`, `TUBE_CONN`, split over several runs (`yang_mills_multilevel_long`) |
| `NONE` | no multilevel observable |

The plaquette of the flux tube is placed at distance `dist_poly/2` along direction 1 and `transv_dist` along
direction 2 from the first loop, with orientation `plaq_dir`. All spatial sizes must be equal. CPU only, not with the
multicanonical algorithm.

## Summary of physical limits

| Feature | D | Group | Programs |
| --- | --- | --- | --- |
| Wilson action, periodic/twisted/open BC | any D ≥ 2 | U(1), SU(N) | all simulation programs |
| θ term, θ profile | 4 | SU(N) | all simulation programs (θ builds) |
| Topological observables | 4 (SU(N)) or 2 (U(1)) |  | all programs with measurements |
| Multicanonical | 4 (SU(N)) or 2 (U(1)); PTBC programs 4D SU(N) only |  | `local_agf`, `local_pt_agf`, `tuning_pt_mc`, `conf_analysis`, `tracedef` |
| PTBC | 4 | SU(N), N ≥ 2 | `local_pt_agf`, `tuning_pt_mc` |
| Trace deformation | any | any (⌊N/2⌋ couplings) | `tracedef` (CPU) |
| Multilevel | any (equal spatial sizes) | any | `multilevel`, `multilevel_long` (CPU) |
