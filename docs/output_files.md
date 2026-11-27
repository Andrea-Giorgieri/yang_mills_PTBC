# Output and data files

This page describes the files read and written by the programs. File names come from the input file
([input_file.md](input_file.md)); the code adds suffixes in some cases (listed below).

Notation: D = space-time dimension, L_μ = lattice sizes, V = volume, N = number of colours, n = number of PTBC replicas.

## Suffixes added to file names

| Suffix | When |
|---|---|
| `_replica_<a>` | PTBC configuration and twist files of replica index a; measurement files of replica a with `--enable-meas-replica` |
| `_back0`, `_back1` | alternating backup copies of configuration and twist files (`saveconf_back_every`) |
| `_step_<index>` | numbered configuration and twist files (`saveconf_analysis_every`); saved multicanonical potentials (tuning) |
| `_dir<μ>` | Polyakov loop density, one file per direction |
| `_agf`, `_gf`, `_cool_<type>`, `_hot` | output files of `yang_mills_conf_analysis` |

## Configuration file (`conf_file`)

Written by `write_conf_on_file_with_name`, read by `read_gauge_conf_from_file_with_name`
([lib/gauge_conf_def.c](../lib/gauge_conf_def.c)).

1. A text header line: `D L_0 … L_{D−1} update_index md5`. The MD5 hash of the links is present when the code is
   built with OpenSSL (`HASH_MODE`) and is checked when the file is read; `conf_check` checks it from the command line.
2. The links in binary, **big-endian** doubles on every machine:
   - sites in lexicographic order, direction 0 running fastest;
   - for each site, the D links μ = 0 … D−1;
   - each link: U(1) → Re, Im (2 doubles); SU(2) → a_0, a_1, a_2, a_3 with U = a_0 + i Σ a_k σ_k (4 doubles);
     SU(N) → Re, Im of each element, row by row (2N² doubles).

When reading, D and the sizes must match the input file.

## Twist file (`twist_file`)

Text file written next to every configuration (`write_twist_on_file_with_name`), read when `start 2`:

```
D L_0 ... L_{D-1}
mu nu x_mu x_nu        # plane and position of the twisted plaquette (0 1 0 0 if there is no twist)
obc_dir x_obc          # open boundary direction (-1 if none) and position of the cut
t_0 t_1 ... t_{D-1}    # translation of the configuration modulo L (PTBC)
```

The twist position and the open boundary position follow the configuration when it is translated. The translation is
read back only in θ builds with a θ profile.

## Data file (`data_file`)

Text file, one line per measurement. Header (two comment lines):

```
# D L_0 ... L_{D-1}
# upd_index <hot columns> ( <smoothed columns> ) x <k> agfrepeat each dt = <dt> ( ... ) x <k> gfrepeat ... ( ... ) x <k> coolrepeat ... mc_topcharge mc_weight
```

Columns of a line, in order:

1. `upd_index`: update index of the configuration.
2. **Hot measurements**, for each active flag: `plaqs plaqt` (spatial and temporal plaquette), `clover_energy`,
   `charge`, `polyre_μ polyim_μ` for μ = 0 … D−1, `polyre_0^k polyim_0^k` for k = 1 … ⌊N/2⌋+1, `multipolyre
   multipolyim`, `charge_prime_μ` for μ = 0 … D−1, `action_1 action_2 action_3 potential` (Wilson part, θ part, total
   action, multicanonical potential).
3. **Smoothed measurements**, for each active method in the order adaptive flow, fixed-step flow, cooling: k blocks
   (k = number of measurements of that method), each with, for each active flag: mean plaquette (single column),
   `clover_energy`, `charge`, Polyakov loops, multi-Polyakov, `charge_prime`, action. The header gives k and the flow
   time (or cooling sweeps) between blocks.
4. Multicanonical builds: `mc_topcharge mc_weight` = charge used by the potential and e^{V(Q)} (reweighting factor).

Plaquettes and Polyakov loops are normalized traces (1/N) Re Tr and (1/N) Im Tr.

### Multilevel programs

`yang_mills_multilevel` and `yang_mills_multilevel_long` write one line per multilevel measurement, **without** the
update index. The second header line depends on `ml_obs`:

| `ml_obs` | Columns |
|---|---|
| `POLYCORR`, `POLYCORR_LONG` | `polycorr`: Polyakov loop correlator at distance `dist_poly`, averaged over space |
| `TUBE_DISC` | `polycorr_re polycorr_im polyplaq_re polyplaq_im`: correlator, and correlator with the plaquette insertion |
| `TUBE_CONN`, `TUBE_CONN_LONG` | `polycorr_re polycorr_im polyplaq_re polyplaq_im polyplaqconn_re polyplaqconn_im`: as `TUBE_DISC`, plus the correlator for the connected flux tube |

With `OPT_MULTIHIT` or `OPT_MULTILEVEL` the lines contain tuning diagnostics instead. The format of `ml_file` (state of
`yang_mills_multilevel_long`) is internal and not described here.

## Binary density files

`energy_density_file`, `charge_density_file`, `polyakov_density_file_dir<μ>`:

1. A text header line: `# D L_0 … L_{D−1}` (Polyakov density: followed by μ).
2. Records in the **native byte order** of the machine, one per measurement:
   - `long` update index, `int` measurement counter (0 = hot, k ≥ 1 = k-th measurement of a smoothing method);
   - energy and charge density: V doubles;
   - Polyakov density along μ: V/L_μ doubles (Re of the normalized trace) followed by V/L_μ doubles (Im).

For each update index the records come in the order hot, adaptive flow, fixed-step flow, cooling (as in the data file).

The sites are in the internal storage order of the code, the even/odd lexicographic index (`cart_to_lexeo` in
[lib/geometry.c](../lib/geometry.c)): with lex = x_0 + L_0 (x_1 + L_1 (x_2 + …)), even sites (Σ x_μ even) are at
lex/2 and odd sites at (lex + V)/2. For the Polyakov density along μ the index is the same even/odd lexicographic
index of the (D−1)-dimensional lattice orthogonal to μ (`lexeo_to_lexeosp_and_mu` in [lib/geometry.c](../lib/geometry.c)).

## Slices

Text files, one line per measurement:

| File | Line |
|---|---|
| `energy_slices_file` | `upd_index meas_count e_0 … e_{L_0−1}`: clover energy density averaged on each slice x_0 |
| `charge_slices_file` (`charge_slices_meas`) | `upd_index meas_count q_0 … q_{L_0−1}`: charge summed on each slice x_0 |
| `charge_slices_file` (`charge_p_slices_meas`) | for each momentum k = 0 … L_ν−1: `upd_index meas_count k Re_0 Im_0 … Re_{L_0−1} Im_{L_0−1}`, with Re/Im = Σ q(x) cos/sin(2πk x_ν/L_ν) on each slice x_0 |

Slices are given in the original coordinates of the configuration (its translation is subtracted). `meas_count` is as
in the density files.

## χ′ file (`chiprime_data_file`)

One line per measurement: `upd_index smoothing meas_count value`, with `smoothing` = 0 hot, 1 gradient flow (both
adaptive and fixed step), 2 cooling.

## Log file (`log_file`)

Written at the end of the run: program name, OpenMP/OpenACC and endianness, PTBC parameters and device distribution,
multicanonical parameters, physical and run parameters, smoothing parameters, timers. Overwritten at each run.

## PTBC files

| File | Content |
|---|---|
| `swap_acc_file` | written at the end of the run: one line per pair `from to c_1 c_2 acc(%) err_acc(%) accepted tried`. Counters refer to the current run |
| `swap_track_file` | header `# MC_step conf_labels [conf_charges]`; one line per step: update index, then for each replica index a the slot of the configuration with index a (and its stored charge in multicanonical builds). Appended on restart |

## Multicanonical files

| File | Content |
|---|---|
| `topo_potential_file` | input potential (format in [multicanonical.md](multicanonical.md#potential)); overwritten by the tuning program |
| `multicanonic_acc_file` | header `# beta L_0 … L_{D−1}`; every `measevery` steps: update index and acceptance (%) for each replica index since the previous line. Appended on restart |
| `<topo_potential_file>_step_<index>` | potentials saved during tuning, same format as the input |

## θ profile file (`theta_profile_file`)

Input only: L_μ numbers f(0) … f(L_μ−1) (μ = `theta_profile_dir`), whitespace separated.
