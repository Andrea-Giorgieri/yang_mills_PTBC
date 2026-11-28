# Link memory layout: AoSoA with blocks of B sites

This page describes the storage layout of the gauge links, the reasons for choosing it, a summary of the benchmarks
behind the choice, and the meaning of the build parameter `Link_block` (B). Numbers are labeled as **measured**
(benchmark or profiler), **derived** (computed from measured numbers) or **estimated** (model or assumption).

Status: the layout is **implemented in the program** (section 6); the performance of the full program with the
new layout has not been measured yet.

## 1. Previous layout

The site order (even/odd lexicographic `si`, update order `d_sip_to_si`, shells for odd sizes) is described in
[architecture_and_portability.md, Data layout](architecture_and_portability.md#data-layout) and
[Lattice parallelization: even/odd and odd sizes](architecture_and_portability.md#lattice-parallelization-evenodd-and-odd-sizes).
Before the AoSoA layout:

- Links: `GC->lattice[si][dir]`, a contiguous `[volume][STDIM]` block of `GAUGE_GROUP` structs (AoS), each padded to a
  multiple of `DOUBLE_ALIGN` = 32 B:

  | group | data per link | `sizeof` | site row (4D) | waste |
  |---|---|---|---|---|
  | U(1) | 16 B | 32 B | 128 B | 50 % |
  | SU(2) | 32 B | 32 B | 128 B | 0 |
  | SU(3) | 144 B | 160 B | 640 B | 10 % |
  | SU(5) | 400 B | 416 B | 1664 B | 4 % |

- `Z` (twist and boundary factors): `[volume][n_planes + 1]` complex numbers, with `n_planes = STDIM (STDIM - 1)` = 12
  in 4D: 13 complex numbers = 208 B per site.
- GPU accesses: one thread per (replica, site), consecutive sites on consecutive threads. A warp loading one element
  touches 32 site rows 640 B apart (SU(3)): ~32 sectors of 32 B per request instead of the 16 needed for 32 loads of
  16 B (**measured** with ncu: 29–31).

## 2. Layouts considered and choice

Indices in units of "element" (one complex number for U(1) and SU(N ≥ 3), one double for SU(2)), `NE` elements per
link (`k` = element, for SU(N) `k = a N + b`), `d` = direction:

| layout | index of element k of link (si, d) | fastest → slowest |
|---|---|---|
| AoS (without struct padding) | `(si STDIM + d) NE + k` | element → direction → site |
| SoA | `(d NE + k) V + si` | site → element → direction |
| **AoSoA(B)** | `((p/B STDIM + d) NE + k) B + p%B` | site in block → element → direction → block |

with `p = si` plus optional padding (below). Re/Im stay adjacent (one `double complex`): a warp reading one element
reads 32 × 16 B = 512 contiguous bytes.

Why **AoSoA**:
- **B = 1 is the AoS** (without struct padding): a single formula, a single kernel code, only a constant changes; no
  directive branches between CPU and GPU.
- **Compile-time constant stride** between the elements of a link (B): address arithmetic folds into constants; the
  SoA has a runtime stride V (extra registers and 64-bit multiplications).
- **Locality**: the data of B consecutive sites lie in `B × STDIM × NE × w` contiguous bytes (SU(3), B = 32: 18 kB),
  good for TLB, prefetch and DRAM pages; the SoA scatters a site over `STDIM NE` planes V elements apart.
- **Even/odd order and checkerboard unchanged**: `si`, neighbor tables, `d_sip_to_si`, rectangle lists do not change;
  blocks do not cross the even/odd boundary if `(V/2) mod B = 0`, otherwise the even half is extended with
  `(-V/2) mod B` padding sites (extra memory ≤ (B−1)/(V/2), negligible except for tiny lattices).
- **MPI**: with replicas distributed across nodes the layout is irrelevant; with domain decomposition, face packing
  needs a kernel anyway (also with the AoS, except for faces orthogonal to the slowest direction); faces per parity are
  multiples of 4 sites in 4D with even sides, so they can start inside a block: same packing kernel, halos aligned to
  B. Hiding all addressing behind an accessor makes a local volume with halos just another V.

Choice: **AoSoA with B = 32 on CPU and GPU**, in-place link access (routines with strided operands), configuration
parameter with default 32.

## 3. The parameter B (`Link_block`)

- **Meaning**: number of consecutive sites of the even/odd order that share a memory block; inside the block the site
  index is the fastest. B = 1 is equivalent to the AoS (without the 32 B struct padding).
- **Usage** (how to set it, default 32, valid values: powers of 2 from 1 to 256, where it ends up): see
  [building.md](building.md#link-block-size).
- **Guidelines** (from the model of section 5 and the benchmarks of section 4):
  - GPU: `B × w ≥ 128 B` (a full cache line per element and per block), i.e. B ≥ 8 for U(1) and SU(N ≥ 3)
    (w = 16 B) and B ≥ 16 for SU(2) (w = 8 B). Beyond that, small gains (B 8 → 32: +0.5–4 %).
  - CPU: the gain grew up to B = 32 (measured); an estimated upper limit comes from the per-core L2 cache
    (section 5.3).
  - Exception: CPU builds of 2D SU(N ≥ 5) should use B = 1 (section 4).
  - Padding: with all sides even in 4D, V/2 is a multiple of 8, so B ≤ 8 never requires padding; B = 16/32 require it
    for lattices with V/2 not a multiple of B, at no measurable cost in the benchmarks.

## 4. Benchmark summary

A stand-alone micro-benchmark reproduces the heatbath and overrelaxation kernels of the PTBC program (same arithmetic
with defect and twist, Philox per site, one thread per (replica, site), one launch per direction and parity) and changes
only the addressing of the link and `Z` elements. All layouts give bitwise identical configurations on the GPU; on the
CPU (icc) they differ in the last bits (different contractions and vectorization), so a CPU implementation must be
validated with a tolerance. "Update" below is (heatbath + 12 overrelaxations)/13.

Hardware and compilers: 4 × NVIDIA A100 (nvc 25.11); 2 × Intel Xeon 8480+ (Sapphire Rapids), 112 OpenMP threads
(icc 2021.10).

Speed-up of the update against the previous layout of the program (AoS with struct padding), geometric mean over several
4D lattice sizes and, for 2D, over L × L lattices with L = 64, 362, 1024 (**measured**):

| group | GPU B = 8 | GPU B = 32 | CPU B = 8 | CPU B = 32 | CPU B = 64 |
|---|---|---|---|---|---|
| U(1) 4D | 3.61 | 3.76 | 1.54 | 1.62 | 1.66 |
| SU(2) 4D | 2.24 | 2.34 | 1.26 | 1.33 | 1.36 |
| SU(3) 4D | – | 1.82 | – | 1.49 | 1.57 |
| SU(5) 4D | – | 1.58 | – | 1.20 | 1.15 |
| U(1) 2D | 1.49 | 1.50 | 1.12 | 1.11 | 1.14 |
| SU(2) 2D | 1.18 | 1.21 | 1.05 | 1.03 | 1.06 |
| SU(3) 2D | 1.80 | 1.85 | 1.06 | 1.05 | 1.04 |
| SU(5) 2D | 2.27 | 2.33 | 0.939 | **0.913** | 0.880 |

- **GPU**: no loss in any group, dimension or size; B = 32 is the best variant or within 2 % of it everywhere. The gain
  saturates at B = 8 (SoA performs like AoSoA). Tiny lattices (less than one wave of threads) gain little.
- **CPU, 4D**: gains for every group, growing with B up to 32; B = 64 is slightly better for U(1), SU(2), SU(3) and
  worse for SU(5). The SoA loses on small lattices.
- **CPU, 2D SU(5)**: every block size B > 1 is slower than the AoS (B = 32 3–14 % slower depending on L), while the
  unpadded AoS (B = 1) is on par. Hence the recommendation `Link_block=1` for CPU builds of 2D SU(N ≥ 5)
  ([building.md](building.md#link-block-size)). SU(4) and 3D were not measured.
- **In-place access** (routines with strided operands) beats copying links into local matrices on both architectures.
- Removing the struct padding alone (B = 1 without padding) is worth 1.32 for U(1), 1.09 for SU(3), 1.02 for SU(5)
  and nothing for SU(2) on the GPU, in line with the padding fractions.

## 5. Why: performance model

### 5.1 The kernels are bandwidth-bound

From ncu profiles (SU(3) and SU(5), 4D), the time per link is close to `D / BW`, with D the DRAM bytes per link and BW
the bandwidth of a plain copy (1.39 TB/s on the A100). With the AoSoA both kernels run at 81–93 % of the copy bandwidth
(**derived**), against 63–78 % with the AoS. The speed-up factors into fewer bytes moved (SU(3): −35 %) and better use
of the bandwidth (SU(5), whose kernel spills registers and is latency-limited with uncoalesced requests). The
arithmetic intensity is ~1.1–1.4 flop/B (**estimated**), well below the A100 FP64 balance point (~7 flop/B): further
gains require fewer bytes per link (e.g. compressed SU(3) links), not another layout.

### 5.2 Coalescing: why the GPU saturates at B = 8

With elements of w bytes, a warp reads 32 w bytes per element, distributed in 32/B contiguous runs of B w bytes
(B ≤ 32). Sectors (32 B) and lines (128 B) touched per request (**derived** from the geometry; with w = 16 B):

| B | 1 (AoS) | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| sectors per request | 32 | 16 | 16 | 16 | 16 | 16 |
| lines per request | 32 | 16 | 8 | 4 | 4 | 4 |

Sectors are already ideal at B = 2, but the speed-up follows the **number of distinct lines** per request and stops
when each run fills a line: `B w ≥ 128 B`, i.e. B ≥ 8 for complex elements, B ≥ 16 for SU(2) (w = 8; measured knee
weaker than predicted). Neighbors not aligned to the block add at most one line per request, the same for every
B ≥ 8: they do not move the optimum. The rule depends neither on N nor on the dimension.

### 5.3 CPU

On the CPU (64 B lines, hardware prefetcher, 2 MB L2 per core) each thread sweeps a range of consecutive sites; with
the AoSoA it reads runs of B w bytes per element and direction, while with the AoS it reads whole site rows of which the
updated direction uses only a part. A reasonable upper limit on B (**estimated**): the per-thread working set, about
`(1 + 2 STDIM)` link blocks of `B × STDIM × NE × w` bytes, must stay well within L2:

| group (4D) | block with B = 32 | working set (× 9) |
|---|---|---|
| U(1) | 2 kB | 18 kB |
| SU(2) | 4 kB | 37 kB |
| SU(3) | 18 kB | 166 kB |
| SU(5) | 51 kB | 461 kB |
| SU(8) | 131 kB | 1.2 MB (at the limit) |

For N ≥ 8–10 on the CPU, B = 16 might be preferable (**estimated**, not measured). The 2D SU(5) loss is not predicted
by this model; a plausible reading (**estimated**) is that in 2D the AoS site row is used almost entirely by each
update, while the AoSoA spreads every SU(5) link over 25 element rows, more concurrent streams per thread than the
prefetchers follow.

## 6. Implementation in the program

Implemented in two steps, both validated: step A changed the storage and accessed every link through
copies into local matrices; step B replaced the copies with in-place access in the kernels (the choice of the
benchmark, section 4: in place beats local copies).

### 6.1 Storage

- **Arrays.** `GC->lattice`, `lattice_copy`, `lattice_cold`, `lattice_copy_cold` and the auxiliary lattices of the
  gradient flow (`Meas_Utils.lattice_aux[k]`) are flat arrays of `GAUGE_ELEM` (`double complex` for U(1) and SU(N),
  `double` for SU(2); `GAUGE_NE` = 1, 4, N² elements per link, defined in `include/gauge_group.h`). `Z` and `Z_copy`
  are flat arrays of `double complex`, `Z_NE` = STDIM (STDIM − 1) + 1 per site. The tables of row pointers are gone.
- **Index.** Element k of link (si, d) is at `((slot / B) STDIM + d) GAUGE_NE B + k B + slot % B`, factor k of Z at
  `((slot / B) Z_NE + k) B + slot % B` (the AoSoA formula of section 2), with `slot = si` for the even sites of the
  lexicographic even/odd order (`si < d_link_even = (V + 1) / 2`) and `slot = si + d_link_pad` for the odd ones,
  `d_link_pad = (−d_link_even) mod B`. The total `d_link_slots` is rounded up to whole blocks. The three values are
  in `GParam` (`init_link_storage_constants`, called by `init_derived_constants`). B = 1 is the AoS without the
  32-byte struct padding.
- **Accessors** (`include/link_layout.h`): `LINK(lattice, si, d, param)` is the address of element 0 of the link
  (the elements are B apart), `ZF(Z, si, k, param)` the factor (an lvalue), `load_link` / `store_link` /
  `copy_link` copy between a link and a local `GAUGE_GROUP`, `lattice_elems` / `Z_elems` the sizes,
  `link_slot_to_si` the inverse map (−1 for padding slots).
- **Allocation.** `allocate_lattice` and `allocate_Z` (`lib/memalign.c`): 256-byte alignment, padding included and
  set to zero, released with `free`. On the device every array is one block (`enter_1d`), so the views and the
  attach of the members are plain pointers. The memory estimate of the devices (`bytes_per_replica`) uses the new
  sizes: an SU(3) lattice is 10 % smaller (144 instead of 160 bytes per link), a U(1) lattice 50 % smaller.
- **Not changed.** `clover_array` (θ term only, not used in production) stays an array of structures; `rng_site`,
  the multilevel arrays and the neighbour tables are not link storage.

### 6.2 Operations on links in place

- **Strided cores.** Every group operation used with lattice operands has one implementation with strided operands,
  `<op>_<group>_s` in `include/u1.h`, `include/su2.h`, `include/sun.h`, `include/sun_upd.h` (element k of X at
  `X[k * sX]`); the function with the original name is the same core with stride 1, so the arithmetic exists once
  and local and in-place operations give the same bits. Operations: `equal_dag`, `minus_equal`,
  `minus_equal_times_real`, `times_equal`, `times_equal_dag`, `times`, `times_dag1`, `times_dag2`, `times_dag12`,
  `plus_equal_times_dag2_complex`, `plus_equal_times_complex`, `norm`, and for SU(N ≥ 3) `single_heatbath`,
  `single_overrelaxation` (with `times_element_SuN_s`, `ennetodue_times_s`). `include/gauge_group.h` dispatches them
  as `<op>_s` on `GAUGE_ELEM` pointers.
- **Named variants** (`include/link_layout.h`): the suffix lists the matrix operands in argument order, result first,
  `l` = local matrix, `f` = link of a lattice (stride B); scalars keep their position and are not listed:
  `times_dag2_lff`, `times_dag12_lff`, `plus_equal_times_dag2_complex_llf`, `plus_equal_times_complex_llf`,
  `times_equal_lf`, `times_equal_dag_lf`, `equal_dag_lf`, `times_lfl`, `times_llf`, `times_dag1_llf`, `times_fll`,
  `times_flf`, `minus_equal_lf`, `minus_equal_times_real_lf`, `minus_equal_ff`, `norm_f`, `unitarize_f`, `cool_fl`,
  `single_heatbath_fl`, `single_overrelaxation_fl`. `unitarize_f`, `cool_fl` and, in U(1) and SU(2), the heatbath
  and overrelaxation copy the link to a local matrix and back (they need the whole matrix; U(1) and SU(2) links have
  one and four elements).
- **Where.** Staples (`calcstaples_wilson`, `calcstaples_wilson_with_defect`), heatbath and overrelaxation, cooling,
  plaquettes, clovers, topological charge, Polyakov loops, gradient flow (fixed and adaptive step), multicanonical
  kernels. The staples of the θ term and of the old multicanonical topological update copy their links to local
  matrices (not production code).

### 6.3 Other parts

- **Accept and restore.** Restore of the whole lattice is a flat copy of all elements; accept unitarizes link by link
  (`accept_link`). The rectangle versions now take `param` (needed by `LINK`).
- **Translation** (`conf_translation`): one thread per storage element; the thread of element k of link (r, d)
  copies element k of link (nnm(r, dir), d) from the copy, consecutive threads on consecutive slots of a block
  (contiguous writes, reads contiguous inside runs of neighbours); padding slots are skipped. Z likewise.
- **I/O.** File format unchanged (links in lexicographic order, each link as before): writing and reading go through
  `load_link` / `store_link`; the md5 sum is computed in file order and does not depend on B.
- **Programs and tests.** `src/debug_link_layout.c` (new): for 4–6 lattice sizes (odd sides and V/2 not a multiple of
  B included) it checks that every element of every link and of Z has its own slot inside the allocation, that the
  unused slots are exactly the padding, the inverse slot map, the copies, and that every named variant gives the
  bits of the operation on local matrices. `src/debug_device_data.c`, `src/debug_twist.c`, `src/conf_txt_to_bin.c`
  adapted.

### 6.4 Device kernels and parallelism

The kernels keep one thread per link (`gang vector`, site index fastest, `ACC_VL_UPD` = 32): the mapping of the
micro-benchmark. With B = 32 a warp reads exactly one block. The optimum of B would change with several threads per
link (vector over the elements of a link, a possible future change): then the lanes read different elements of one link
and B = 1 is the coalesced layout, which `Link_block=1` restores without code changes.

### 10.5 Validation

Local (WSL, 2026-11-28), after step A and again after step B:
- **Unit test** `debug_link_layout` (gcc): all checks passed for (N, D, B) = (3, 4, 1/8/32/256), (1, 4, 32),
  (2, 4, 16/32), (5, 4, 32), (4, 4, 64), (3, 2, 32), (3, 3, 4).
- **Program against HEAD, gcc OpenMP with Philox** (4 threads): `Link_block=1` against HEAD and `Link_block=32` against
  `Link_block=1`, `yang_mills_local_pt_agf` (4D SU(2), SU(3), SU(5); inputs with V/2 a multiple of 32, with padding
  (V/2 = 1080), with an odd side, with twisted start; multicanonical and θ builds) and `yang_mills_local_agf` (also
  U(1), 2D and 3D SU(3), even and odd sides): all saved configurations, swap files and twist files **bitwise
  identical**. `dati.dat` was bitwise identical or equal within 1e-13 relative: the same binary run twice also differs
  in the last digit of some measurements (OpenMP reductions without a fixed order), so the measurements are
  compared with a tolerance. `yang_mills_conf_analysis` on configurations written by HEAD, `debug_twist` and
  `debug_device_data` (73 checks): identical output, all passed.
- **OpenACC schedule** (`scripts/check_acc_schedule.sh`, nvc 26.1): 0 variants with problems (2, 3, 4, 5, 8, 5mc,
  5th). `update_replica_sites` stack/spill stores/spill loads in bytes: SU(3) 176/0/0 (225 registers, 224 before),
  SU(5) 1520/2736/2240 (before 1520/2648/2176), SU(8) 3568/2616/3244 (before 3568/2600/3100).

To do on Leonardo: GPU bitwise check (B = 1 against HEAD, B = 32 against B = 1; the micro-benchmark was bitwise
identical on the GPU for every layout), benchmarks B1, B2, B3 (estimate of section 6.4: update −33 to −36 %),
translation time, an `ACC_VL_UPD` 32/64/128 sweep (the occupancy optimum may move with fewer bytes per link), icc
CPU build with a tolerance check.