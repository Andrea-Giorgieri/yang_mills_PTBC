# GitHub repository description: proposals

Draft text for the GitHub "About" field (tagline, max 100 characters), the short repository
description (max 350 characters) and the repository topics.
Audience: lattice field theorists and HPC practitioners.
Each claim was checked against the code, `AUTHORS`, `README.md` and the git history.

## Variant A: physics first

**Tagline (97 chars):**
Lattice Yang-Mills Monte Carlo with parallel tempering on boundary conditions, OpenMP and OpenACC

**Description (344 chars):**
> Monte Carlo simulations of SU(N) and U(1) Yang-Mills theories on the lattice with the Wilson action. Its main feature is parallel tempering on boundary conditions (PTBC) against topological freezing, optionally with multicanonical sampling. One C source runs on CPUs with OpenMP and on GPUs with OpenACC. By C. Bonati, C. Bonanno, A. Giorgieri.

## Variant B: HPC first

**Tagline (89 chars):**
SU(N) and U(1) lattice gauge theory in C, for CPUs (OpenMP) and GPUs (OpenACC), with PTBC

**Description (348 chars):**
> C code for lattice Yang-Mills simulations (heatbath, overrelaxation) that builds for multicore CPUs (OpenMP) or GPUs (OpenACC) from one source. It implements parallel tempering on boundary conditions (PTBC); the replicas can be split over several GPUs. Topological observables with gradient flow and cooling. By C. Bonati, C. Bonanno, A. Giorgieri.

## Variant C: feature list

**Tagline (79 chars):**
Yang-Mills on the lattice: PTBC, multicanonical, twisted BC; OpenMP and OpenACC

**Description (339 chars):**
> Lattice Monte Carlo for SU(N) and U(1) gauge theories: Wilson action; periodic, open and twisted boundary conditions; imaginary theta; trace deformation; gradient flow and cooling. Topological freezing is addressed by PTBC and multicanonical sampling. CPU (OpenMP) and GPU (OpenACC) from one source. By C. Bonati, C. Bonanno, A. Giorgieri.

## GitHub topics

`lattice-gauge-theory`, `yang-mills`, `lattice-qcd`, `monte-carlo-simulation`,
`parallel-tempering`, `openacc`, `openmp`, `gpu-computing`

The code is pure gauge, not QCD. `lattice-qcd` is included because that topic is where most of
this audience looks; drop it to be strict.

## Check before publishing

- **The GPU support isn't on `master` yet.** OpenACC exists only on the `openacc_gpu_port` branch.
  If the description goes up before the merge, the OpenACC claim describes a branch, not the
  default code.
- **GPU testing.** The checks behind the feature list ran on the CPU only (OpenACC compiled for
  the host). Back the GPU claim with production runs on real GPUs (e.g. on Leonardo) before stating it.
- **Credits.** They follow `README.md`: base code (`yang-mills`) by C. Bonati; PTBC and
  multicanonical by C. Bonanno; twisted boundary conditions by A. Giorgieri. The GPU-port commits
  are by A. Giorgieri, so Variant B could add "GPU port: A. Giorgieri" if space allows.
- **Dimensions.** "D dimensions" is left out on purpose. The PTBC program describes itself as 4D,
  so a general-D claim next to PTBC would overstate it.
