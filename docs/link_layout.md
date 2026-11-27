# Link memory layout: AoSoA with blocks of B sites

Note from phase 9c of the OpenACC port (November 2026). It collects the layout choices, the tests, the results and their
interpretation, and documents the build parameter `Link_block` (B). Throughout the text every number is labeled as
**measured** (from a job or a profile), **derived** (computed from measured numbers with an explicit formula) or
**estimated** (assumption or unverified model).

Status: layout chosen and measured with a micro-benchmark; **not yet implemented in the program** (design review
first). The `Link_block` parameter already exists in `configure.ac`/`build.sh`, but the code does not use it yet.

## 1. Current layout

The site order (even/odd lexicographic `si`, update order `d_sip_to_si`, shells for odd sizes) and the padded struct
sizes are described in [architecture_and_portability.md, Data layout](architecture_and_portability.md#data-layout) and
[Lattice parallelization: even/odd and odd sizes](architecture_and_portability.md#lattice-parallelization-evenodd-and-odd-sizes).
What matters here:

- Links: `GC->lattice[si][dir]`, a contiguous `[volume][STDIM]` block of `GAUGE_GROUP` structs (AoS), each padded to a
  multiple of `DOUBLE_ALIGN` = 32 B (**measured** with gcc and nvc):

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

## 2. Proposed layouts and choice

Indices in units of "element" (one complex number for U(1) and SU(N ≥ 3), one double for SU(2)), `NE` elements per
link (`k` = element, for SU(N) `k = a N + b`), `d` = direction:

| layout | index of element k of link (si, d) | fastest → slowest |
|---|---|---|
| AoS (today; without struct padding) | `(si STDIM + d) NE + k` | element → direction → site |
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
- **MPI**: with replicas distributed across nodes (phase 10 idea) the layout is irrelevant; with domain decomposition,
  face packing needs a kernel anyway (also with the AoS, except for faces orthogonal to the slowest direction); faces
  per parity are multiples of 4 sites in 4D with even sides, so they can start inside a block: same packing kernel,
  halos aligned to B. Hiding all addressing behind an accessor makes a local volume with halos just another V.

Choice (user, 2026-11-27): **AoSoA with B = 32 on CPU and GPU**, in-place link access (routines with strided operands),
configuration parameter with default 32.

## 3. The parameter B (`Link_block`)

- **Meaning**: number of consecutive sites of the even/odd order that share a memory block; inside the block the site
  index is the fastest. B = 1 is equivalent to the AoS (without the 32 B struct padding).
- **Usage** (how to set it, default 32, valid values: powers of 2 from 1 to 256, where it ends up): see
  [building.md](building.md#link-block-size).
- **Constraints and guidelines** (from the model in section 7, with the available checks):
  - GPU: `B × w ≥ 128 B` (a full cache line per element and per block), i.e. B ≥ 8 for U(1) and SU(N ≥ 3)
    (w = 16 B) and B ≥ 16 for SU(2) (w = 8 B). Beyond that, no measured gain (B 8 → 32: +0.5–2 %).
  - CPU: the gain grew up to B = 32 (measured); an estimated upper limit comes from the per-core L2 cache
    (section 7).
  - Padding: with all sides even in 4D, V/2 is a multiple of 8, so B ≤ 8 never requires padding; B = 16/32 require it
    in 18 of the 66 inputs in `prod_inputs/` (few sites per parity, no measurable cost in case M5).
- **Status**: defined but not yet used by the code (layout to be implemented after the design review).

## 4. Tests performed

Micro-benchmark `bench9/layout/layout_bench.c` (independent of the program):
- heatbath and overrelaxation arithmetic with defect and twist **transcribed from the program** (same operations in the
  same order: `include/{random,u1,u1_upd,su2,su2_upd,sun,sun_upd}.h`, `calcstaples_wilson_with_defect`,
  `heatbath_with_defect`, `overrelaxation_with_defect`), Philox per site, one thread per (replica, site) and one launch
  per direction and parity like `update_replica_sites`; only the addressing of link and `Z` elements changes;
- variants `<group>[d<D>]_<layout>_<access>`: layout `b1p` (AoS with struct padding: the program's layout), `b1` (AoS
  without padding), `b2`…`b64` (AoSoA), `soa`; access `inplace` (links read and updated in place, as now) or `local`
  (copied into local matrices, as before 9a.2a);
- output: one CSV row per (run, kernel) with all parameters (group, dimension, volume, replicas, layout, block, access,
  CPU/GPU, device, compiler, ns per link with median/min/max over 3 repetitions, copy bandwidth, configuration hash,
  plaquette);
- validation: bit hash of all links after thermalization (equal for every layout if the configurations are bitwise
  identical) and plaquette;
- "upd" = (hb + 12 or)/13, the production update (12 overrelaxations per heatbath in all PTBC scale-setting inputs,
  57 of the 66 in `prod_inputs/`; θ with OBC and multicanonical use 4–8).

Main campaign (`B9L_SET=main`): SU(3) and SU(5) 4D, cases from `prod_inputs/` with replicas = production / 4 GPUs:
S3 24·16²·24 ×5, M3 36·20²·36 ×5, L3 46·24²·46 ×9 (= B1), S5 10·8²·10 ×3, M5 26·14²·26 ×6 (V/2 not a multiple of 16),
L5 40·18²·40 ×8; layouts `b1`…`b32`, `soa`, in place and local. GPU job 65365326 (1 Booster node, 4 A100, nvc 25.11,
84 runs + 28 ncu profiles) and CPU job 65365395 (1 DCGP node, 2 × Xeon 8480+, 112 OpenMP threads, icc 2021.10,
84 runs), debug queues, both `JOB RESULT: OK`, nothing skipped. Note: in this campaign the reference is `b1` (AoS
**without** struct padding), slightly more favorable than the program's true AoS; the comparison with `b1p` is in the
edge-case campaign (section 7).

Reproduction (from the repository root, on Leonardo):

    bash scripts/build_bench9_layout_leonardo.sh                  # GPU (nvc) and CPU (icc): build_b9_layout/{gpu,cpu}
    sbatch scripts/job_leonardo_bench9_layout_gpu.sh              # timings (4 GPUs in parallel) + ncu
    sbatch scripts/job_leonardo_bench9_layout_cpu.sh              # DCGP, 112 threads

Analysis and plots (locally, Python 3 with matplotlib for the plots):

    python bench9/layout/analyze_layout.py bench9/results/b9c_layout_gpu_b9c_20261127_155314_65365326.csv \
        bench9/results/b9c_layout_cpu_b9c_20261127_155523_65365395.csv --out bench9/results/b9c_layout_summary_20261127.csv
    python bench9/layout/plot_layout.py bench9/results/b9c_layout_gpu_b9c_20261127_155314_65365326.csv \
        bench9/results/b9c_layout_cpu_b9c_20261127_155523_65365395.csv \
        --ncu bench9/results/b9c_layout_gpu_b9c_20261127_155314_65365326_ncu.csv --outdir bench9/results/plots_b9c_20261127

## 5. Results (main campaign)

Speed-up of the update (hb + 12 or)/13 versus `b1_inplace`, geometric mean over the 3 cases per N (**measured**;
spread across repetitions: median 0.2 % GPU, 0.4–1.2 % CPU):

| | B 2 | B 4 | B 8 | B 16 | B 32 | SoA | B 32 local |
|---|---|---|---|---|---|---|---|
| GPU SU(3) | 1.31 | 1.55 | 1.69 | 1.71 | 1.71 | 1.71 | 1.59 |
| GPU SU(5) | 1.26 | 1.44 | 1.55 | 1.56 | 1.57 | 1.56 | 1.51 |
| CPU, 6 cases | 1.09 | 1.23 | 1.26 | 1.30 | 1.33 | 1.05 | 1.29 |

Absolute times of `b1_inplace` (ns per link, upd): GPU SU(3) 2.94–3.03, SU(5) 16.4–17.2; CPU SU(3) 9.0–12.4,
SU(5) 22.6–33.7. The micro-benchmark kernel is ~20–25 % lighter than the program's (SU(3) hb 3.32 versus
4.37 ns per link of B6 in 9a.2a), with the same register, occupancy and spill regime.

![Speed-up vs B, GPU](../bench9/results/plots_b9c_20261127/speedup_vs_block_gpu.png)

![Speed-up vs B, CPU](../bench9/results/plots_b9c_20261127/speedup_vs_block_cpu.png)

![ns per link vs volume, GPU](../bench9/results/plots_b9c_20261127/ns_vs_volume_gpu.png)

![GPU vs CPU](../bench9/results/plots_b9c_20261127/gpu_vs_cpu.png)

![ncu metrics](../bench9/results/plots_b9c_20261127/ncu.png)

- **Validation**: on the GPU all variants are bitwise identical after thermalization; on the CPU (icc) the variants
  other than `b1_inplace` differ in the last bits (plaquette equal to 15 digits): different contractions and
  vectorization with different addresses. A CPU implementation therefore requires validation with a tolerance.
- **GPU**: the gain saturates at B = 8; SoA = AoSoA; SU(5) gains almost as much as SU(3).
- **CPU**: no regression for B ≥ 2 (minimum +3–5 % in the tiny case S5); the gain grows up to B = 32
  (1.44–1.46 on large lattices); the SoA loses on small lattices (0.86 in S5).
- **Access**: in place beats local copies with AoSoA on both architectures.

ncu (one launch, direction 1, even sites, cases M3 and M5; **measured**):

| N | layout | kernel | ns/link | sectors per load request | DRAM bytes per link | DRAM % | L2 % | registers |
|---|---|---|---|---|---|---|---|---|
| 3 | b1 | hb | 3.33 | 28.9 | 3357 | 61.7 | 91.4 | 218 |
| 3 | b8 | hb | 1.94 | 15.8 | 2195 | 69.3 | 86.1 | 227 |
| 3 | b32 | hb | 1.91 | 15.8 | 2195 | 70.2 | 86.3 | 227 |
| 3 | b1 | or | 2.81 | 30.8 | 3070 | 67.1 | 95.9 | 218 |
| 3 | b8 | or | 1.66 | 16.4 | 2056 | 76.2 | 91.0 | 227 |
| 5 | b1 | hb | 20.6 | 29.2 | 18104 | 53.8 | 72.5 | 255 |
| 5 | b8 | hb | 12.1 | 15.9 | 14896 | 75.7 | 78.0 | 255 |
| 5 | b1 | or | 17.1 | 31.3 | 16477 | 59.4 | 80.4 | 255 |
| 5 | b8 | or | 10.9 | 16.4 | 14097 | 79.5 | 78.4 | 255 |

## 6. Interpretation

### 6.1 Bandwidth limit (roofline)

The time per link cannot drop below `t_min = D / BW`, with D the DRAM bytes per link (ncu) and BW the bandwidth of a
plain copy of the same array measured by the benchmark (1.39 TB/s on the A100, ~85 % of the theoretical peak of
1.63 TB/s; measured with host timing, launch included: the copy kernels of the 9b translation, timed by nsys, reach
1.43 TB/s, see [swap_translation_gpu.md](swap_translation_gpu.md)).
Efficiency `η = t_min / t` (**derived** from measured numbers; ncu times are at locked clocks, so the efficiency is
approximate):

| N | kernel | AoS: t_min / t → η | AoSoA B 8: t_min / t → η | byte ratio D_AoS / D_AoSoA | η ratio | speed-up |
|---|---|---|---|---|---|---|
| 3 | hb | 2.41 / 3.33 → 72 % | 1.58 / 1.94 → 81 % | 1.53 | 1.12 | 1.72 |
| 3 | or | 2.21 / 2.81 → 78 % | 1.48 / 1.66 → 89 % | 1.49 | 1.14 | 1.69 |
| 5 | hb | 13.0 / 20.6 → 63 % | 10.7 / 12.1 → 89 % | 1.22 | 1.40 | 1.70 |
| 5 | or | 11.9 / 17.1 → 69 % | 10.1 / 10.9 → 93 % | 1.17 | 1.35 | 1.57 |

The speed-up factors exactly into (fewer bytes) × (better bandwidth use):
- **SU(3)**: mostly **fewer bytes** (−35 %). With the AoS a warp brings 32 partial links into cache per load; the parts
  of the links not yet used are probably evicted from L1/L2 before they are needed and reloaded (**estimated**: not
  verified with reuse metrics).
- **SU(5)**: mostly **better bandwidth use** (63–69 % → 89–93 %): the kernel has 255 registers and ~7.4 KB of spill per
  link, with few requests in flight per thread; coalesced requests (16 sectors instead of 32) reduce the latency of
  each one, so that the limit becomes the DRAM bandwidth rather than the latency.
- After the change both kernels run at 81–94 % of the copy bandwidth: they are **bandwidth-bound**. Further gains
  require fewer bytes per link (compressing SU(3) links to 12 reals, more staple reuse across threads), not another
  layout.
- **Arithmetic intensity** (**estimated**: hand count of operations): SU(3) ~3 kflop per link (staples ~2.4 kflop,
  heatbath ~0.8) over 2.2 KB → ~1.4 flop/B; SU(5) ~16 kflop over 15 KB → ~1.1 flop/B. The A100 balance point in
  FP64 (9.7 TFLOP/s / 1.39 TB/s) is ~7 flop/B: the kernels are memory-bound by a factor ~5, consistent with
  SM at 20–35 %.

### 6.2 Coalescing: why the GPU saturates at B = 8

With elements of w bytes, a warp reads 32 w bytes per element, distributed in 32/B contiguous runs of B w bytes
(B ≤ 32). Sectors (32 B) and lines (128 B) touched per request (**derived** from the geometry; with w = 16 B):

| B | 1 (AoS) | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| sectors per request | 32 | 16 | 16 | 16 | 16 | 16 |
| lines per request | 32 | 16 | 8 | 4 | 4 | 4 |
| measured GPU SU(3) speed-up | 1 | 1.31 | 1.55 | 1.69 | 1.71 | 1.71 |

Sectors are already ideal at B = 2, but the speed-up follows the **number of distinct lines** per request (L1 work,
L2 requests, DRAM access granularity) and stops when each run fills a line: `B w ≥ 128 B`.
Rule (**derived**, verified only for w = 16 B): B ≥ 128 / w, i.e. B ≥ 8 for complex elements, B ≥ 16 for SU(2) (w = 8).
Neighbors not aligned to the block (shifts of 0/1 in direction 0, of L0/2 in direction 1, …) add at most one line per
request, the same for every B ≥ 8: they do not move the optimum.

### 6.3 CPU

On the CPU (64 B lines, hardware prefetcher, 2 MB L2 cache per core) each thread sweeps a range of consecutive sites;
with the AoSoA it reads runs of B w bytes per element and per direction, while with the AoS it reads for each site a
whole row of which the updated direction uses only a part. The gain grows with B up to 32 (**measured**); the
explanation (prefetch on longer runs, lines fully used) is **estimated**, without hardware counters. A reasonable
upper limit (**estimated**): the per-thread working set, about `(1 + 2 STDIM)` link blocks
(`B × STDIM × NE × w` bytes each), must stay well within L2:

| group (4D) | block with B = 32 | working set (× 9) |
|---|---|---|
| U(1) | 2 kB | 18 kB |
| SU(2) | 4 kB | 37 kB |
| SU(3) | 18 kB | 166 kB |
| SU(5) | 51 kB | 461 kB |
| SU(8) | 131 kB | 1.2 MB (at the limit) |

For N ≥ 8–10 on the CPU, B = 16 might be preferable (**estimated**, not measured).

### 6.4 Expected effect on the program

**Estimated** (not measured): the update kernels are 82–94 % of the update after 9b; with the micro-benchmark speed-up
for B = 32 (L3 1.63 for B1, S3 1.81 for B3, M5 1.54 for B2, which is not among the benchmark cases; for B2 the ~92 %
share of update kernels measured in 9a.0 is assumed) the update would drop by ~36 % in B1 and B3 and by ~33 % in B2.
The program has the AoS **with** padding (`b1p`), so the true gain could be somewhat larger. Flow and measurements
read the same arrays (9d), not measured here.

## 7. Choice of B in edge cases: model and targeted tests

The model of sections 6.2–6.3 gives rules without simulations:
1. **GPU**: B ≥ 128 / w (w = 16 B for U(1) and SU(N ≥ 3), 8 B for SU(2)); beyond that, flat. It depends neither on N
   nor on the dimension. B = 32 satisfies all groups.
2. **CPU**: B as large as possible as long as `(1 + 2 STDIM) × B × STDIM × NE × w` stays ≪ L2 per core; in 2D the
   working set is half (fewer directions and neighbors).
3. **Lattice geometry**: matters only through padding ((−V/2) mod B sites per parity, relevant only if V/2 is a few
   thousand sites) and unaligned neighbors (at most one extra line per request, the same for every B ≥ 8). The extent
   along the fastest direction (half row L0/2) does not enter the choice of B, because the sites of one parity are
   contiguous in `si` also across rows.
4. **Tiny volumes** (less than one wave of threads on the GPU, ~27 000 threads with 255 registers): limited by launches
   and latencies; the layout matters little.

Model limits: the real bandwidth efficiency (η) is not predicted by the model, only bytes and lines; the CPU behavior is
semi-empirical (no counters). The edge-case tests (section 7.1) confirm rules 1, 3 and 4 on the GPU (with a weaker
SU(2) knee than predicted) and rule 2 in 4D, and contradict the 2D part of rule 2 for SU(5) on the CPU.

Targeted tests (`B9L_SET=edge`, in-place access only, 4 replicas, defect 4, debug queues), designed to check exactly
the model predictions not already covered by the main campaign:
- **U(1) and SU(2) in 4D** on the shapes of S3/M3/L3 (`U1_S/M/L`, `SU2_S/M/L`): GPU prediction: U(1) saturates at
  B = 8 and gains a lot from removing the padding (today 32 B per link, of which 16 B are data); SU(2) saturates only
  at B = 16.
- **2D with L × L, L = 64, 362, 1024** for U(1), SU(2), SU(3), SU(5) (same L for all groups): L = 64 tiny
  (V/2 = 2048), L = 362 with odd half rows (181) and V/2 = 65 522 not a multiple of 4 (padding for every B ≥ 4),
  L = 1024 large and aligned. Prediction: same thresholds as in 4D, reduced gains for L = 64.
- **Layouts**: `b1p` (reference: the program today), `b1`, `b2`…`b64`, `soa`; plus `b1p`, `b1`, `b32`, `b64` on cases
  M3, L3, M5, L5 (effect of struct padding in SU(3)/SU(5) and B = 64 for the CPU).
- 178 runs per architecture; GPU estimated at ~5 minutes plus 24 ncu profiles (U1_M, SU2_M, SU3D2_1024 with b1p, b8,
  b16, b32); CPU in layout priority order with a time check.

Commands (on Leonardo):

    B9L_SET=edge bash scripts/build_bench9_layout_leonardo.sh      # build_b9_layout_edge/{gpu,cpu}
    export B9_TAG=b9c_edge_<date>                                   # same tag for the two jobs
    B9L_SET=edge sbatch scripts/job_leonardo_bench9_layout_gpu.sh
    B9L_SET=edge sbatch scripts/job_leonardo_bench9_layout_cpu.sh

Analysis: `python bench9/layout/analyze_layout.py --ref b1p_inplace bench9/results/b9c_edge_*_<tag>_*.csv` (section
"B = 32 against b1p_inplace": flags LOSS if B = 32 is more than 3 % slower than the program's layout) and
`plot_layout.py` with `--ref b1p_inplace`.

### 7.1 Results of the edge cases

GPU job 65371450 (178 runs + 24 ncu profiles) and CPU job 65371466 (178 runs), tag `b9c_edge_20261127`, both
`JOB RESULT: OK`, nothing skipped; analysis in `bench9/results/b9c_edge_summary.txt`. Reproduce the tables and figures
with:

    python bench9/layout/analyze_layout.py --ref b1p_inplace bench9/results/b9c_edge_gpu_b9c_edge_20261127_65371450.csv \
        bench9/results/b9c_edge_cpu_b9c_edge_20261127_65371466.csv --out bench9/results/b9c_edge_summary_20261127.csv
    python bench9/layout/plot_layout.py --ref b1p_inplace bench9/results/b9c_edge_gpu_b9c_edge_20261127_65371450.csv \
        bench9/results/b9c_edge_cpu_b9c_edge_20261127_65371466.csv \
        --ncu bench9/results/b9c_edge_gpu_b9c_edge_20261127_65371450_ncu.csv \
        --variants b1p_inplace,b8_inplace,b32_inplace,soa_inplace --outdir bench9/results/plots_b9c_edge_20261127

Upd speed-up against the program's layout (`b1p`), geometric mean over the cases (**measured**):

| group | GPU b8 | GPU b32 | GPU best | CPU b8 | CPU b32 | CPU b64 |
|---|---|---|---|---|---|---|
| U(1) 4D | 3.61 | 3.76 | 3.77 (b64) | 1.54 | 1.62 | 1.66 |
| SU(2) 4D | 2.24 | 2.34 | 2.34 (b32) | 1.26 | 1.33 | 1.36 |
| SU(3) 4D (M3, L3) | – | 1.82 | 1.83 (b64) | – | 1.49 | 1.57 |
| SU(5) 4D (M5, L5) | – | 1.58 | 1.59 (b64) | – | 1.20 | 1.15 |
| U(1) 2D | 1.49 | 1.50 | 1.53 (soa) | 1.12 | 1.11 | 1.14 |
| SU(2) 2D | 1.18 | 1.21 | 1.21 (b32) | 1.05 | 1.03 | 1.06 |
| SU(3) 2D | 1.80 | 1.85 | 1.88 (soa) | 1.06 | 1.05 | 1.04 |
| SU(5) 2D | 2.27 | 2.33 | 2.35 (b64) | 0.939 | **0.913** | 0.880 |

![Speed-up against B, GPU, edge cases](../bench9/results/plots_b9c_edge_20261127/speedup_vs_block_gpu.png)

![Speed-up against B, CPU, edge cases](../bench9/results/plots_b9c_edge_20261127/speedup_vs_block_cpu.png)

- **GPU**: no loss in any group, dimension or size; B = 32 is the best variant or within 2 % of it everywhere. Tiny
  lattices behave as predicted: U(1) 2D with L = 64 (8192 threads per launch, less than one wave) gains only 1.04.
- **Model check on the GPU** (rule 1): the knee is at B = 8 for complex elements, as predicted, but B 8 → 32 still adds
  1–4 % (more than "flat"). For SU(2) (w = 8 B) the predicted knee at B = 16 is weak: +3 % from B = 8 to 16, because
  ncu shows that B = 8 already reaches the same load sectors per request (11.1 against 11.0). 2D follows the same
  thresholds as 4D.
- **Padding of the program's structs** (`b1` against `b1p`): U(1) 1.32, SU(3) 1.09, SU(5) 1.02, SU(2) 1.00 on the GPU,
  in line with the padding fractions (50 %, 10 %, 4 %, 0). Against the real layout of the program the main-set gains
  become SU(3) 1.82 and SU(5) 1.58 on the GPU, SU(3) 1.49 and SU(5) 1.20 on the CPU.
- **CPU, 4D**: gains for every group. B = 64 is slightly better than 32 for U(1), SU(2), SU(3) and worse for SU(5)
  (1.15 against 1.20), consistent with the working-set estimate of section 6.3 (rule 2).
- **CPU, 2D**: small gains for U(1), SU(2), SU(3); a **loss for SU(5)**: B = 32 at 0.973 (L = 64, cache-resident),
  0.906 (L = 362) and 0.863 (L = 1024), every AoSoA block size and SoA losing; the unpadded AoS (`b1`) is at 1.01.
  The model did not predict it (it expected 2D to be favourable). A plausible reading (**estimated**, no hardware
  counters): in 2D the AoS row of a site (2 links) is used almost entirely by the update of either direction, so the
  AoS wastes little, while the AoSoA spreads every SU(5) link over 25 element rows 512 B apart, more concurrent streams
  per thread than the prefetchers follow. The U(1) 2D case L = 362 (0.95) is cache-resident (copy bandwidth
  5.6 TB/s) and its differences are small.

## 8. U(1) and SU(2) groups

They use different structs (`U1`: one complex number, `Su2`: four reals) and their own routines (`u1.h`, `su2.h`,
`u1_upd.h`, `su2_upd.h`), transcribed in the micro-benchmark with the same generic addressing. In the AoSoA layout the
element is the complex number (U(1)) or the double (SU(2)): same scheme, same B.

Result (**measured**, section 7.1): no loss on either architecture, and the largest gains of all groups on the GPU
(U(1) 3.6–3.9× in 4D: the struct is 32 B for 16 B of data, and DRAM bytes per link drop from 1365 to 385; SU(2)
2.2–2.4×). No fallback layout and no different B are needed for these groups.

The only configuration that loses is a CPU build of 2D SU(5) (section 7.1). Decision (user, 2026-11-27): the default
stays B = 32 everywhere, and the build documentation recommends `Link_block=1` for CPU builds of 2D SU(N ≥ 5)
([building.md](building.md#link-block-size)). SU(4) and 3D were not measured.

## 9. Files

- Micro-benchmark: `bench9/layout/layout_bench.c`; analysis `bench9/layout/analyze_layout.py`; plots
  `bench9/layout/plot_layout.py`.
- Scripts: `scripts/bench9_layout_common.sh`, `scripts/build_bench9_layout_leonardo.sh`,
  `scripts/job_leonardo_bench9_layout_{gpu,cpu}.sh`.
- Results: `bench9/results/b9c_layout_*` (main campaign), `b9c_layout_summary_20261127.csv`,
  `plots_b9c_20261127/`, `b9c_summary.txt`; edge cases: `b9c_edge_*`, `b9c_edge_summary_20261127.csv`,
  `plots_b9c_edge_20261127/`, `b9c_edge_summary.txt` (index: `bench9/results/README.md`).
- Parameter: `configure.ac` (`Link_block`), `include/macro.h` (`LINK_BLOCK`), `build.sh` (`LINK_BLOCK`).
