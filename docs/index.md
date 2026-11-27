# Documentation index

Documentation of **yang&#95;mills&#95;PTBCMC**, a C code for Monte Carlo simulations of lattice Yang–Mills theories with
parallel tempering on boundary conditions (PTBC) and multicanonical sampling of the topological charge, running on
CPUs (OpenMP) and NVIDIA GPUs (OpenACC) from a single source.

Topics are organized in pages, linked to each other. For an overview of the project and its features see the top-level [README.md](../README.md).

| Page | Content |
| --- | --- |
| [building.md](building.md) | Compile-time options (`configure`, `build.sh`, `Link_block`), compiler presets, which programs are built, compile-time limits |
| [physics_and_parameters.md](physics_and_parameters.md) | Gauge groups and dimensions, action (Wilson, θ term and θ profile, trace deformation), boundary conditions (periodic, twisted, open), observables, gradient flow and cooling, physical limits |
| [algorithms.md](algorithms.md) | Heatbath, overrelaxation, Metropolis; structure of a Monte Carlo step in each program; multilevel; random numbers; start, checkpoints and restart |
| [ptbc.md](ptbc.md) | Parallel tempering on boundary conditions: defect, replicas, hierarchical update, swaps, translations, outputs, parallelization |
| [multicanonical.md](multicanonical.md) | Multicanonical bias potential in the topological charge, its combination with PTBC, automatic tuning |
| [architecture_and_portability.md](architecture_and_portability.md) | Single source for OpenMP and OpenACC, data layout, device memory management, distribution of replicas among GPUs, even/odd parallelization and odd lattice sizes |
| [programs.md](programs.md) | Usage of each executable: command line, input sections, outputs |
| [input_file.md](input_file.md) | Reference of all input keywords with defaults and constraints |
| [output_files.md](output_files.md) | Formats of configuration, twist, data, density, slice, PTBC and multicanonical files |
| [swap_p1p2.md](swap_p1p2.md) | Derivation of the PTBC swap action difference in terms of two sums per configuration |
| [swap_translation_gpu.md](swap_translation_gpu.md) | GPU implementation of the swap sums, replica index copies and translation: block sums, asynchronous queues, coalesced copies |
| [link_layout.md](link_layout.md) | Storage layout of the gauge links (planned AoSoA with blocks of B sites): choice, benchmarks, performance model, meaning of `Link_block` |

Benchmark results of the GPU port: [bench8/results/README.md](../bench8/results/README.md) and
[bench9/results/README.md](../bench9/results/README.md).


## Where to start

- **Run a simulation:** [building.md](building.md) → [programs.md](programs.md) → [input_file.md](input_file.md).
- **Understand what is simulated:** [physics_and_parameters.md](physics_and_parameters.md), [ptbc.md](ptbc.md), [multicanonical.md](multicanonical.md).
- **Work on the code or port it:** [architecture_and_portability.md](architecture_and_portability.md), [algorithms.md](algorithms.md).


## Topic owners

| Topic | Authoritative page |
| --- | --- |
| Feature summary, project overview, credits | [README.md](../README.md) |
| Configure options, `build.sh`, `Link_block` usage, compile-time limits | [building.md](building.md) |
| List of programs, command lines, what each program does | [programs.md](programs.md) |
| Input keywords and defaults | [input_file.md](input_file.md) |
| Output file formats | [output_files.md](output_files.md) |
| Physical limits (groups, dimensions, observables) | [physics_and_parameters.md](physics_and_parameters.md) |
| Update algorithms, Monte Carlo step, random numbers, start and restart | [algorithms.md](algorithms.md) |
| PTBC semantics (defect, swaps, translations, rectangles) | [ptbc.md](ptbc.md) |
| Maths of the swap sums | [swap_p1p2.md](swap_p1p2.md) |
| GPU implementation of swaps and translation | [swap_translation_gpu.md](swap_translation_gpu.md) |
| Data layout, struct sizes, even/odd order, device memory, multi-GPU | [architecture_and_portability.md](architecture_and_portability.md) |
| Planned link layout and the choice of B | [link_layout.md](link_layout.md) |
