# Multicanonical algorithm in the topological charge

This page describes the multicanonical (bias-potential) algorithm, its combination with parallel tempering on boundary
conditions (PTBC) and the automatic tuning of the potential.

## Background

When the topological susceptibility is small (large volumes, high temperature), configurations with Q ≠ 0 are rare and
the distribution of Q is hard to sample. The multicanonical algorithm adds a bias potential V(Q) to the action, which
flattens (or reshapes) the distribution of Q, and removes it afterwards by reweighting each measurement with
e^{+V(Q)}. Here the potential can be different for each PTBC replica, V_a(Q), so that it also improves the swap
acceptances.

## Activation

| Selection | What |
|---|---|
| compile time | `--enable-use-multicanonical` (`MULTICANONICAL_MODE`) |
| runtime | `topo_potential_file` (required), `multicanonic_acc_file` (required), `grid_step`, `grid_max`, `topo_cooling`, `topo_coolsteps`, `topo_alpha` |

Programs: `yang_mills_local_agf`, `yang_mills_local_pt_agf`, `yang_mills_tracedef`, `yang_mills_conf_analysis`
(reweighting factors only), `yang_mills_tuning_pt_mc` (requires it). Not available in the multilevel programs.
Allowed only for D = 4 SU(N) (and D = 2 U(1)); the PTBC programs require D = 4 SU(N).

## Charge used by the potential

The potential is a function of a charge x computed from the configuration:

- `topo_cooling 0`: x = Q, the clover topological charge of the configuration;
- `topo_cooling 1`: x = Q after `topo_coolsteps` cooling sweeps (ordering `LEX_DIR_LEXEO_SITE`) of a copy of the
  configuration; the configuration itself is not changed;
- `topo_alpha α > 0`: α-rounding, x → round(α x); `topo_alpha 0` disables it.

The charge of the accepted configuration is stored (`stored_topcharge`) and reused in the swaps and in the output.

## Potential

V_a(x) is piecewise linear on the grid x_i = −grid_max + i·grid_step, i = 0 … n_grid − 1 with
n_grid = 2·grid_max/grid_step + 1, and constant outside the grid (equal to the value at the nearest end).
`compute_topo_potential` in [lib/gauge_conf_multicanonic.c](../lib/gauge_conf_multicanonic.c).

`topo_potential_file` format (text, whitespace separated; read by `read_topo_potential` in
[lib/gparam.c](../lib/gparam.c)):

```
grid_step_coarse  w_0 w_1 ... w_{n-1}
x_0  V_0(x_0) V_1(x_0) ... V_{n-1}(x_0)
x_1  V_0(x_1) ...
...                                     (n_grid rows)
```

- n = `N_replica_pt` (1 for the programs without parallel tempering): one column of V per replica.
- The x values must coincide with the grid defined by `grid_step` and `grid_max` (checked).
- `grid_step_coarse` and the weights w_a are used only by the tuning program (below); the other programs read them and
  ignore them.

## Where the potential enters

1. **Accept/reject after the updates.** The heatbath and overrelaxation updates do not see V. After a full update
   step (or, in PTBC, after the update of a hierarchical rectangle), the new charge x_new is computed and the whole
   step is accepted with probability min(1, e^{V_a(x_old) − V_a(x_new)}); if rejected, the configuration is restored
   to its state before the step. So the sampled distribution is ∝ e^{−S − V_a(x)}.
   - After a rectangle update, only the change of the charge inside an enlarged rectangle is computed (border
     R_i + 2, plus 4·`topo_coolsteps` with cooling), and x_new = x_old + Δx.
   - The tests are done on the host in slot order with the host random stream; charges and cooling run on the device
     holding each configuration.
2. **PTBC swaps.** The swap of the configurations A (index a) and B (index b) adds
   ΔV = V_a(x_B) + V_b(x_A) − V_a(x_A) − V_b(x_B) to the action difference ([ptbc.md](ptbc.md)).
3. **Output for reweighting.** Each line of the data file ends with `mc_topcharge mc_weight` = x and e^{V_a(x)} for the
   measured configuration (a = its replica index). Expectation values are ⟨O⟩ = ⟨O e^{V}⟩_V / ⟨e^{V}⟩_V.

Implementation: `multicanonic_metropolis_replica`, `multicanonic_metropolis_rectangle_replica`,
`multicanonic_metropolis_conf`, `delta_topo_potential_swap` in
[lib/gauge_conf_multicanonic.c](../lib/gauge_conf_multicanonic.c); calls in `update`, `accept_or_restore_replica`,
`update_with_trace_def` ([lib/gauge_conf_upd.c](../lib/gauge_conf_upd.c)) and `swap_group`
([lib/gauge_conf_paral_temp.c](../lib/gauge_conf_paral_temp.c)).

## Acceptance file

`multicanonic_acc_file`: header `# beta L_0 ... L_{D-1}`; then, every `measevery` steps, the update index followed by
the acceptance (%) of the multicanonical test for each replica index since the previous line. If the file already
exists it is appended to (whatever `start` is) and no new header is written.

## Tuning the potential: `yang_mills_tuning_pt_mc`

The program runs PTBC exactly as `yang_mills_local_pt_agf` and, after `thermal` steps, modifies the potentials at every
step (`tune_topo_potential`):

- **n ≥ 3 replicas.** Two grids are built: a coarse one V_coarse, defined at multiples of `grid_step_coarse` and
  linearly interpolated, from the charge of the physical replica (index 0), and a fine one V_fine, with step
  `grid_step`, from the charge of the last replica (index n − 1). At each step the value s is added to the grid around
  the current charge (split linearly between the two nearest grid points; also at −x if `topo_tuning_even 1`). The
  potential of replica a is V_a = w_a V_coarse + (1 − w_a) V_fine.
  Every time more than `topo_tuning_weights_every` swaps of the first pair have been proposed, the weights w_a of the
  intermediate replicas are shifted to equalize the swap acceptances (w_0 and w_{n−1} stay fixed, all w in [0, 1]),
  and the swap counters are reset.
- **n < 3 replicas.** Each V_a is updated independently with the charge of the configuration that has index a.
- All potentials are shifted so that V_a(0) = 0.
- **Step size.** Every replica starts with s = `topo_tuning_stp`. When every replica has crossed the whole grid
  (from one end to the other), each s_a is multiplied by max(0.5, 1 − ⟨ΔV⟩_a/3), with ⟨ΔV⟩_a = s_a·(steps of the
  crossing)/n_grid, and the program prints the new steps. The tuning ends when the largest ⟨ΔV⟩_a is below
  `topo_tuning_thr`, or at `sample` steps or at the wall time.
- **Output.** With `topo_tuning_save_every k > 0`, the potentials are written every k steps to
  `<topo_potential_file>_step_<index>`. At the end the final potentials (and weights) are written to
  `topo_potential_file`, **overwriting the input file**.

| Parameter | Meaning | Default |
|---|---|---|
| `topo_tuning_stp` | initial value added to the potential at each step | 0.1 |
| `topo_tuning_thr` | stop when the largest mean change per crossing is below this | 0.1 |
| `topo_tuning_weights_every` | update the weights after this many proposed swaps | 1000 |
| `topo_tuning_even` | 1: keep the potential even in x | 1 |
| `topo_tuning_save_every` | save the potentials every this many steps (0 = never) | 0 |

The tuning program needs an initial `topo_potential_file` (for example all zeros, with the chosen `grid_step_coarse`
and weights).

Implementation: `tune_topo_potential`, `tune_topo_potential_unconstrained`, `update_tuning_stp` in
[lib/gauge_conf_multicanonic.c](../lib/gauge_conf_multicanonic.c); main in
[src/yang_mills_tuning_pt_mc.c](../src/yang_mills_tuning_pt_mc.c). After tuning, the potential is copied to the GPUs
(`copy_topo_grid_to_device`).

## Limits

- The potential is a function of the global charge only.
- The bias is applied by accept/reject of whole update steps, so the acceptance drops if V changes much over the
  typical change of x in one step; monitor `multicanonic_acc_file`.
- Not available with the multilevel programs (they stop at start).
