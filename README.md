# yang_mills_PTBCMC

Monte Carlo simulations of lattice SU(N) Yang–Mills theories with **parallel tempering on boundary conditions (PTBC)** and **multicanonical sampling of the topological charge (MC)**.
The U(1) group and any number of dimensions are also supported with standard algorithms. A single C99 source runs on CPUs with OpenMP and on NVIDIA GPUs with OpenACC.

Detailed documentation: [docs/index.md](docs/index.md).

## Key features

- **Parallel tempering on boundary conditions** (4D SU(N)). Replicas that differ only by the boundary conditions on a small
  defect, interpolating between periodic and open, are swapped with a Metropolis test; hierarchical updates around
  the defect and random translations of the physical replica. Also works with the θ term and the multicanonical algorithm. → [docs/ptbc.md](docs/ptbc.md)
- **Multicanonical algorithm in the topological charge**, with a different bias potential for each PTBC replica and a
  program that tunes the potentials automatically. → [docs/multicanonical.md](docs/multicanonical.md)
- **Twisted ('t Hooft), open and periodic boundary conditions**, an **imaginary θ term** included in the update, also
  with a position-dependent θ profile.
- **Gradient flow** with fixed or adaptive step size, and **cooling**, for smoothed observables.
- **One source, two architectures:** OpenMP on CPU or OpenACC on one or more GPUs per process, with even/odd parallel
  updates that also support odd lattice sizes.

## Supported physical parameters

| Gauge groups | U(1), SU(2), SU(N) for any N ≥ 3 (compile time) |
| Dimension | any D ≥ 2 (compile time; PTBC only in 4D); lattice sizes even or odd (runtime) |
| Action | Wilson plaquette action; optional imaginary θ term (4D) with optional θ profile; trace deformation of the Polyakov loop (dedicated program) |
| Boundary conditions | periodic, twisted on any plane, open along any direction, PTBC defect |
| Observables | plaquette, clover energy density, topological charge (4D SU(N), 2D U(1)) with its density, time slices, Q′ and χ′, Polyakov loops (also powers and multi-direction loops), Polyakov-loop correlators and flux tubes with the multilevel algorithm |
| Smoothing | Wilson flow (fixed or adaptive step), cooling (five sweep orderings) |

Details and limits: [docs/physics_and_parameters.md](docs/physics_and_parameters.md).

## Algorithms

- Heatbath: Kennedy–Pendleton/Creutz for SU(2), Cabibbo–Marinari over SU(2) subgroups for SU(N), exact heatbath for
  U(1); overrelaxation; Metropolis for the trace-deformed temporal links.
- PTBC with hierarchical updates, swaps and translations.
- Multicanonical accept/reject on the topological charge (optionally measured after cooling), with reweighting factors
  written to the data file.
- Multilevel algorithm for Polyakov-loop correlators and flux tubes.
- Random numbers: dSFMT (one generator per thread) or the counter-based Philox (one stream per site, always used on GPU).

Details: [docs/algorithms.md](docs/algorithms.md).

## Architecture and portability

- **OpenMP** (`--enable-use-openmp`) or **OpenACC** (`--enable-use-openacc`, NVIDIA GPUs) from the same source: every
  parallel loop carries both directives, and the same kernels run on the host or on the GPU.
- **Several GPUs per process:** each PTBC replica lives on one GPU for the whole run, replicas are spread evenly over
  the GPUs, and only scalars move between GPUs. `yang_mills_gpu_plan` estimates how many GPUs an input needs.
- **Even/odd parallelization with odd sizes:** the largest even sub-lattice is updated in parallel by parity, the
  remaining border sites sequentially.
- Compilers: gcc, Intel icc, NVIDIA HPC SDK (nvc); presets for the CINECA clusters Leonardo (CPU and GPU),
  Galileo100 and Marconi in `build.sh`.

Details: [docs/architecture_and_portability.md](docs/architecture_and_portability.md).

## Building

The gauge group, the dimension and the optional features are fixed at compile time:

```bash
autoreconf -fi
./configure N_c=3 ST_dim=4 --enable-use-openmp Num_threads=16 CC=gcc CFLAGS="-O3"
make -j
```

Main options: `N_c`, `ST_dim`, `--enable-use-openmp` + `Num_threads`, `--enable-use-openacc` + `Num_devices`,
`--enable-use-theta`, `--enable-use-multicanonical`, `--enable-meas-replica`, `--enable-rng-philox`, `--enable-debug`.
Alternatively, set the options at the top of `build.sh` and run it.
All options and limits: [docs/building.md](docs/building.md).

`Link_block` (default 32, a power of 2 from 1 to 256) sets the block size B of the planned AoSoA storage of the gauge
links (not used by the code yet): usage in [docs/building.md](docs/building.md#link-block-size), rationale and
benchmarks in [docs/link_layout.md](docs/link_layout.md).

## Quick start

1. Build the code (see above).
2. Run a program without arguments: it prints its usage and the compile-time configuration, and writes a template
   input file `template_input.example` in the current directory.
```bash
   ./yang_mills_local_pt_agf
```
3. Edit the template (lattice size, β, replicas, measurements, file names) and run:
```bash
   ./yang_mills_local_pt_agf input_file
```

Input keywords: [docs/input_file.md](docs/input_file.md). Output formats: [docs/output_files.md](docs/output_files.md).
Production inputs used so far are in [`prod_inputs/`](prod_inputs).

## Main programs

PTBC simulations (`yang_mills_local_pt_agf`, `yang_mills_tuning_pt_mc`), single-configuration simulations
(`yang_mills_local_agf`), analysis of saved configurations (`yang_mills_conf_analysis`), GPU planning
(`yang_mills_gpu_plan`), multilevel and trace-deformation programs (CPU only) and configuration-file utilities. The
list of programs, what each does and its command line: [docs/programs.md](docs/programs.md).

## Documentation

All pages are listed, with their content, in [docs/index.md](docs/index.md).

## Credits and license

Based on the `yang-mills` package by Claudio Bonati (claudio-bonati/yang-mills). PTBC and the χ′ routines by Claudio Bonanno;
twisted boundary conditions, new observables/smoothing routines and GPU porting by Andrea Giorgieri.

Released under the GNU General Public License v3 ([LICENSE](LICENSE)).
