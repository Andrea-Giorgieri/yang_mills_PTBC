# Swaps and translation on the GPU: block sums, asynchronous queues, coalesced copies

This note describes how the swap sums $P_1$, $P_2$, the copy of `replica_index` after the swaps, and the translation of
the periodic replica (`conf_translation`) are implemented on the device. It covers `swap_plaquette_sums_replica`,
`conf_translation` (`lib/gauge_conf_paral_temp.c`), `copy_replica_index_to_device`, `alloc_double_on_devices`
(`lib/device_data.c`). The mathematics of the sums $P_1$, $P_2$ is in [swap_p1p2.md](swap_p1p2.md); the physics and
semantics of swaps and translations are in [ptbc.md](ptbc.md).

## 1. Context (as in the code)

- The $N_r$ configurations live in the slots `GC[g]`, $g = 0,\dots,N_r-1$; device $d$ holds the slots
  `device_slot_begin(d)` $\le g <$ `device_slot_end(d)`. A single host thread drives all devices, one at a time
  (`select_device`).
- All compute regions are asynchronous in the `QUEUE_DEVICE` queue of their device: the host enqueues them and moves
  on; it stops only at a wait (`acc_wait`, `acc_wait_all`) or at a synchronous copy.
- A PTBC step calls `swap()` after the full update and after each sweep of the hierarchical rectangles (typically
  10–40 times per step), and each time, if $c_0 = 1$, `conf_translation` on the configuration with replica index 0
  (see [ptbc.md](ptbc.md)).
- Each `swap()` needs the sums $P_1(g)$, $P_2(g)$ of all slots on the host (the Metropolis tests use the host random
  number generator), then copies the new `replica_index` values to the devices, because the defect coefficients read
  by the kernels depend on them.
- Configuration fields (`lattice`, `Z`, their copies `lattice_copy`, `Z_copy`, and `lattice_cold`,
  `lattice_copy_cold` in the multicanonical case) are contiguous blocks of rows, one row per site; see
  [architecture_and_portability.md](architecture_and_portability.md#data-layout) for the layout. Row widths used below
  (SU(3), 4D): a `lattice` row is $W = 640$ B (4 links of 160 B, padding included), a `Z` row is
  $13 \cdot 16 = 208$ B ($n_{\rm planes} = 12$).

## 2. Swap sums

### 2.1 Block sums

A reduction into a host scalar is available only after a wait on the device queue, so summing slot by slot would
serialize the devices (launch → kernel → copy → wait for each slot). Instead, `swap_plaquette_sums_replica` splits
`swap_rect` into chunks of `SWAP_SUM_CHUNK` = 1024 sites ($n_{\rm ch} = \lceil V_{\rm rect}/1024\rceil$ chunks) and
numbers the (slot, chunk) pairs as $q = g\, n_{\rm ch} + c$. On each device:
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
(OpenMP build, or computation forced on the host) `swap_plaquette_sums` is used for each slot.

### 2.2 Copy of `replica_index`

The new indices are copied with `acc_update_device_async(..., QUEUE_DEVICE)` for each slot, without waits. This is
correct because:
- all compute regions use `QUEUE_DEVICE`, so in the queue the copy comes after the kernels that read the old value and
  before those that read the new one;
- the host modifies `replica_index` only in the next swap, after the waits of section 2.1 on all devices: at that point
  the previous copies have completed.

### 2.3 Benefit

Host synchronization points per call of `swap()`:

| | slot by slot | block sums |
|---|---|---|
| sums | $N_r$ waits, in sequence | 1 wait per device, with all devices already working |
| `replica_index` | $N_r$ waits + $N_r$ synchronous copies | 0 (asynchronous copies in the queue) |
| kernel launches for the sums | $2N_r$ (reduction + final reduction kernel) | 1 per device |

The sum kernels of all devices run concurrently and the GPU does not sit idle between one slot and the next. The gain
matters mostly for small volumes, where launch and synchronization latencies dominate.

### 2.4 Rounding

The chunks are summed on the host in a fixed order, instead of a compiler reduction tree over the whole rectangle.
$P_1$, $P_2$ therefore differ from a direct sum by a relative $O(\varepsilon)$ ($\varepsilon \approx 10^{-16}$), as
discussed in [swap_p1p2.md](swap_p1p2.md), section 5, item 4. In addition:
- the result is **reproducible** from run to run (same number of chunks, same order);
- it does not depend on the number of GPUs: the chunks are per slot, and the host sum follows the slot index;
- a Metropolis decision can change only if the random number falls within $O(\varepsilon)$ of the threshold.

## 3. Translation

### 3.1 What it does

`conf_translation` shifts the configuration by one lattice step in the direction $+\hat{d}$ (chosen at random; see
[ptbc.md](ptbc.md) for its role in PTBC):
$$U'(r) = U(r - \hat d), \qquad Z'(r) = Z(r - \hat d)$$
(and likewise for `lattice_cold` in the multicanonical case), reading from the copies (`lattice_copy`, `Z_copy`, equal
to the accepted configuration) and writing into `lattice`, `Z`. It then realigns the copies: `lattice_copy`
$\leftarrow$ `lattice`, `Z_copy` $\leftarrow$ `Z`. It runs only on the device of the configuration with replica
index 0; the other devices wait for it at the next swap, so it is on the **critical path** of the step.

### 3.2 Implementation: whole-row copies

Copying link by link (one thread per link, consecutive threads one site row apart) is not coalesced and reaches only a
small fraction of the bandwidth. Instead, each field is a contiguous block of rows, and the translation moves **whole
rows**. The two loops are flat, with one index $u$ per double:
$$\texttt{dst}[u] = \texttt{src}[\,\texttt{nnm}(r, d)\cdot w + (u - r w)\,], \qquad r = \lfloor u / w \rfloor,$$
with $w$ = doubles per row ($w = 80$ for `lattice` in SU(3), i.e. `STDIM·sizeof(GAUGE_GROUP)/sizeof(double)`;
$2(n_{\rm planes}+1) = 26$ for `Z`). The index runs sequentially through the segments `lattice`, [`lattice_cold`,]
`Z`, so a single kernel covers all fields. The copy back is the same loop without the neighbor:
`lattice_copy[u] = lattice[u]`, `Z_copy[u] = Z[u]`.

Consecutive threads read and write consecutive doubles: the writes are always coalesced, the reads are coalesced
within each source row (640 contiguous B), and break only at the boundary between two rows. The cost is that of two
passes over the lattice at nearly full bandwidth.

The loop is the same code on the host (OpenMP) and on the device; `_Static_assert` checks that `GAUGE_GROUP` is made
of doubles. The loop-invariant scalars are declared `firstprivate`.

### 3.3 Rounding

These are pure copies: the result is **bitwise identical** to a link-by-link copy, on CPU and on GPU.

## 4. Checks

- CPU (gcc, OpenMP) and OpenACC on the host (devices simulated, `num_gpus` 1 and 2): configurations bitwise identical
  to the slot-by-slot implementation for SU(2), SU(3), SU(5) and SU(3) with θ term and multicanonical algorithm, in
  `yang_mills_local_pt_agf` and `yang_mills_local_agf`.
- `scripts/check_acc_schedule.sh`: correct OpenACC schedule (sums with `loop gang` over the chunks and
  `loop vector(128)` with reduction, excluded from the generic check via `GANG_VECTOR_OK`; translation with two
  `gang, vector(128)` loops).
- GPU (4 × A100): configurations and observables bitwise identical to the slot-by-slot implementation for SU(3) and
  SU(5) with hierarchical PTBC and for the multicanonical algorithm (no Metropolis decision was close enough to the
  threshold to change).

## 5. Performance summary

Measured on 4 × NVIDIA A100 with SU(3) and SU(5) PTBC benchmarks, compared with the slot-by-slot sums and the
link-by-link translation:

- The translation reaches 1.1–1.4 TB/s (67–88 % of the A100 bandwidth, about 7× faster than link-by-link copies) and
  accounts for most of the gain, because it is on the critical path.
- The PTBC update step is 14–33 % faster (the largest gain on the smallest lattice); with one GPU the saving is kernel
  time, with several GPUs it is also time removed from the critical path.
- The block sums remove almost all context synchronizations; on the update step they are worth ~1–4 %.
- Steps with measurements gain less (2–10 %), because the gradient flow dominates them.
