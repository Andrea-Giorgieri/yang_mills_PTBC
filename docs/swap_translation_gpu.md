# Swaps and translation on the GPU: block sums, asynchronous queues, coalesced copies

This note describes how the swap sums $P_1$, $P_2$, the copy of `replica_index` after the swaps, and the translation of
the periodic replica (`conf_translation`) are implemented on the device, and why they are faster than the previous
version. It covers `swap_plaquette_sums_replica`, `conf_translation` (`lib/gauge_conf_paral_temp.c`),
`copy_replica_index_to_device`, `alloc_double_on_devices` (`lib/device_data.c`) (porting phase 9b). The mathematics of
the sums $P_1$, $P_2$ is in [swap_p1p2.md](swap_p1p2.md) and does not change; the physics and semantics of swaps and
translations are in [ptbc.md](ptbc.md).

## 1. Context (as in the code)

- The $N_r$ configurations live in the slots `GC[g]`, $g = 0,\dots,N_r-1$; device $d$ holds the slots
  `device_slot_begin(d)` $\le g <$ `device_slot_end(d)`. A single host thread drives all devices, one at a time
  (`select_device`).
- All compute regions are asynchronous in the `QUEUE_DEVICE` queue of their device: the host enqueues them and moves
  on; it stops only at a wait (`acc_wait`, `acc_wait_all`) or at a synchronous copy.
- A PTBC step calls `swap()` after the full update and after each sweep of the hierarchical rectangles (12 times per
  step in the B1 and B2 inputs, 33 in B3), and each time, if $c_0 = 1$, `conf_translation` on the configuration with
  replica index 0 (see [ptbc.md](ptbc.md)).
- Each `swap()` needs the sums $P_1(g)$, $P_2(g)$ of all slots on the host (the Metropolis tests use the host random
  number generator), then copies the new `replica_index` values to the devices, because the defect coefficients read
  by the kernels depend on them.
- Configuration fields (`lattice`, `Z`, their copies `lattice_copy`, `Z_copy`, and `lattice_cold`,
  `lattice_copy_cold` in the multicanonical case) are contiguous blocks of rows, one row per site; see
  [architecture_and_portability.md](architecture_and_portability.md#data-layout) for the layout. Row widths used below
  (SU(3), 4D): a `lattice` row is $W = 640$ B (4 links of 160 B, padding included), a `Z` row is
  $13 \cdot 16 = 208$ B ($n_{\rm planes} = 12$).

## 2. Swap sums (S)

### 2.1 Before

For each slot, in sequence: `swap_plaquette_sums` selects the device of the slot, launches a reduction
(`reduction(+:p1,p2)` into a host scalar) over `swap_rect` and **waits** on the queue before moving to the next slot.
The result of a reduction into a host scalar is available only after the wait, and the local variable does not outlive
the function: it is not possible to enqueue all reductions and wait only once.

Consequences, for each call of `swap()`:
- $N_r$ launch → kernel → copy → wait cycles, one after the other; the kernel of device $d+1$ starts only after the
  wait on device $d$: the devices never work concurrently;
- each cycle pays the latency of the launch, of the reduction kernel (about 20–25 µs plus the final `__red` kernel) and
  of the synchronization, even though the computation is small.

### 2.2 Now

`swap_plaquette_sums_replica` splits `swap_rect` into chunks of `SWAP_SUM_CHUNK` = 1024 sites
($n_{\rm ch} = \lceil V_{\rm rect}/1024\rceil$ chunks) and numbers the (slot, chunk) pairs as
$q = g\, n_{\rm ch} + c$. On each device:
1. **one kernel** over all pairs $q$ of its slots: `loop gang` over the chunks, inside each gang a `loop vector` with
   `reduction(+:p1,p2)` over the 1024 sites of the chunk (128 lanes, 8 sites per lane); the gang writes
   `part[2q]` $= p_1$, `part[2q+1]` $= p_2$;
2. **one copy** `update self` of its part of `part` (its slots only), in the same queue, without waiting.

Only after enqueuing kernels and copies on **all** devices does the host wait, once per device. Then it sums the chunks
in a fixed order:
$$P_n(g) = \sum_{c=0}^{n_{\rm ch}-1} \texttt{part}[2(g\,n_{\rm ch} + c) + n - 1], \qquad n = 1, 2.$$

The buffer `part` = `acc_counters->swap_partial` ($2 N_r n_{\rm ch}$ doubles) is a host array present on every device
(`alloc_double_on_devices`, `acc_create` device by device); each device uses and copies only the part of its own
slots, so the copies of the different devices write disjoint regions of the host array. It is allocated at the first
swap on the device (its size depends on $V_{\rm rect}$) and freed in `free_acc_utils`. Without data on the device
(OpenMP build, or computation forced on the host) `swap_plaquette_sums` is used for each slot, unchanged.

### 2.3 Copy of `replica_index`

Before: for each slot, `acc_wait_all` on its device and then a **synchronous** 4 B copy (`acc_update_device`).
Now: `acc_update_device_async(..., QUEUE_DEVICE)` for each slot, without waits. This is correct because:
- all compute regions use `QUEUE_DEVICE`, so in the queue the copy comes after the kernels that read the old value and
  before those that read the new one;
- the host modifies `replica_index` only in the next swap, after the waits of section 2.2 on all devices: at that point
  the previous copies have completed.

### 2.4 Benefit

Host synchronization points per call of `swap()`:

| | before | now |
|---|---|---|
| sums | $N_r$ waits, in sequence | 1 wait per device, with all devices already working |
| `replica_index` | $N_r$ `acc_wait_all` + $N_r$ synchronous copies | 0 (asynchronous copies in the queue) |
| kernel launches for the sums | $2N_r$ (reduction + `__red`) | 1 per device |

Per step, with the production replica counts on 4 GPUs: B3 (18 replicas, 33 swaps) from about 1200 to 132
synchronizations; B1 (34 replicas, 12 swaps) from about 800 to 48. (The benchmark runs of section 5 use reduced replica
counts.) The sum kernels of the four devices run concurrently and the GPU does not sit idle between one slot and the
next.

The gain matters mostly for small volumes: the phase 9a.0 nsys profile of B3 shows ~41 µs of idle GPU per launch,
mostly host waits around the swaps.

### 2.5 Rounding

The summation order changes (chunks of 1024 sites summed on the host in a fixed order, instead of the nvc reduction
tree over the whole rectangle). $P_1$, $P_2$ therefore differ by a relative $O(\varepsilon)$
($\varepsilon \approx 10^{-16}$), as discussed in [swap_p1p2.md](swap_p1p2.md), section 5, item 4. In addition:
- the result is **reproducible** from run to run (same number of chunks, same order);
- it does not depend on the number of GPUs: the chunks are per slot, and the host sum follows the slot index;
- a Metropolis decision changes only if the random number falls within $O(\varepsilon)$ of the threshold: the
  configurations stay identical in practice, but this is not guaranteed. For this reason the comparison with the
  previous version also uses the observables (section 4).

## 3. Translation (T)

### 3.1 What it does

`conf_translation` shifts the configuration by one lattice step in the direction $+\hat{d}$ (chosen at random; see
[ptbc.md](ptbc.md) for its role in PTBC):
$$U'(r) = U(r - \hat d), \qquad Z'(r) = Z(r - \hat d)$$
(and likewise for `lattice_cold` in the multicanonical case), reading from the copies (`lattice_copy`, `Z_copy`, equal
to the accepted configuration) and writing into `lattice`, `Z`. It then realigns the copies: `lattice_copy`
$\leftarrow$ `lattice`, `Z_copy` $\leftarrow$ `Z`. It runs only on the device of the configuration with replica
index 0; the other devices wait for it at the next swap, so it is on the **critical path** of the step.

### 3.2 Before

Two kernels over $(n_{\rm planes}+1)\cdot V$ threads: thread $(j, r)$ copied with `equal` the **whole link**
`lattice[r][j]` (160 B in SU(3), 20 doubles) and one $Z$ factor. Consecutive threads (sites $r$, $r+1$) access links
one row apart (640 B): each load instruction of a warp touches 32 different rows, and 20 instructions are needed per
link. The access is not coalesced and the bandwidth stays low.

nsys measurements (phase 9a.0): in B1 ($V = 1.22 \cdot 10^6$, $\sim 2.1$ GB read + written per kernel) 11.7 ms for
the gather and 12.5 ms for the copy back, i.e. $\sim 170$ GB/s, about 10 % of the A100 bandwidth. The copy back, which
has no neighbor jump, is as slow as the gather: the limit is the access pattern, not the neighbor.
Per step: B1 ~290 ms (~11 % of the update, ~60 % of the idle GPU time), B2 ~260 ms (~6 %), B3 ~98 ms (~20 %).

### 3.3 Now

Each field is a contiguous block of rows, and the translation moves **whole rows**. The two loops are flat, with one
index $t$ per double:
$$\texttt{dst}[u] = \texttt{src}[\,\texttt{nnm}(r, d)\cdot w + (u - r w)\,], \qquad r = \lfloor u / w \rfloor,$$
with $w$ = doubles per row ($w = 80$ for `lattice` in SU(3), i.e. `STDIM·sizeof(GAUGE_GROUP)/sizeof(double)`;
$2(n_{\rm planes}+1) = 26$ for `Z`). The index $t$ runs
sequentially through the segments `lattice`, [`lattice_cold`,] `Z`, so a single kernel covers all fields. The copy
back is the same loop without the neighbor: `lattice_copy[u] = lattice[u]`, `Z_copy[u] = Z[u]`.

Consecutive threads read and write consecutive doubles: the writes are always coalesced, the reads are coalesced
within each source row (640 contiguous B), and break only at the boundary between two rows. The expected cost is that
of two passes over the lattice at nearly full bandwidth: ~2.5–3 ms in B1 instead of 24 ms.

The loop is the same code on the host (OpenMP) and on the device; `_Static_assert` checks that `GAUGE_GROUP` is made
of doubles. The loop-invariant scalars are declared `firstprivate`.

### 3.4 Rounding

These are pure copies: the result is **bitwise identical** to the previous one, on CPU and on GPU.

## 4. Checks

- CPU (gcc OpenMP, commit `5fd5b2c` against the modified tree, then committed as `826d579`, same seed, 6·4²·6,
  4 replicas, 12 steps): bitwise identical in 7 cases out of 7 (SU(2), SU(3), SU(5), SU(3) θ + multicanonical;
  `local_pt_agf` and `local_agf`). The CPU path of the sums is unchanged; the translation uses the new code on the
  host too.
- `scripts/check_acc_schedule.sh` (nvc 26.1): correct schedule in all variants; sums with `loop gang` over the chunks
  and `loop vector(128)` with reduction (function excluded from the checks via `GANG_VECTOR_OK`), translation with
  two `gang, vector(128)` loops.
- OpenACC multicore (devices simulated on the host, `num_gpus` 1 and 2): bitwise identical to `5fd5b2c`.
- GPU (Leonardo, 4 A100, reference `5fd5b2c` against the new version in the same job; analysis in
  `bench8/results/b9b_summary.txt`): configurations and observables bitwise identical in B1 and B3 (SU(3),
  hierarchical PTBC), B2 (SU(5), hierarchical PTBC) and B4 (multicanonical); with these seeds no Metropolis decision
  was close enough to the threshold to change.

## 5. Results

Campaign `b9b_20261127_131556` (jobs 65356629, 65356670, 65356675); 4 GPUs unless stated otherwise. The replica
counts are the reduced benchmark ones (B1 16, B2 12, B3 8, B4 12), not the production ones of section 2.4.

| | translation per call | translation per step | update per step |
|---|---|---|---|
| B1 (SU(3), $1.22\cdot 10^6$ sites, 16 replicas) | 24.2 → 3.3 ms | 291 → 40 ms | 1.30 → 1.04 s (−20 %) |
| B2 (SU(5), 12 replicas) | | | 1.60 → 1.37 s (−14 %) |
| B3 (SU(3), $1.47\cdot 10^5$ sites, 8 replicas) | 2.95 → 0.40 ms | 98 → 13 ms | 0.282 → 0.189 s (−33 %) |
| B3, 1 GPU | | | 0.500 → 0.406 s (−19 %) |
| B4 (multicanonical, 1 swap per step) | | | 0.156 → 0.152 s (−3 %); 1 GPU −1 % |

- The translation reaches 1.09–1.43 TB/s, 67–88 % of the A100 bandwidth (previously ~170 GB/s; per site 640 B of
  links and 208 B of `Z`) and accounts for almost all of the gain: it was on the critical path.
- The block sums remove almost all context synchronizations (B3: `cuCtxSynchronize` from 1065 to 13 in 4 steps) and
  one third to one half of the event synchronizations; on the update they are worth ~1–4 %. The block kernel has a
  slightly larger GPU time than the per-slot kernels (few gangs), but it runs on all devices concurrently: negligible
  difference.
- Fraction of the update spent in update kernels: B1 74 → 94 %, B3 54 → 82 %.
- In B3 the saving is the same with 1 and with 4 GPUs (0.094 and 0.093 s per step): with one GPU it is kernel time,
  with four it is time removed from the critical path (the other GPUs were waiting for the translation).
- On the step with measurements the gain is small (B1/B2 −2 %, B3 −7/−10 %): the measurements (flow) dominate.
