# APCEMM documentation index

This directory contains source-tree, physics, numerics, and issue-audit notes for the APCEMM C++ codebase under `Code.v05-00/`.

## Start here

- [Architecture](ARCHITECTURE.md) - high-level map of the executable, source tree, build structure, major subsystems, and OpenMP/case-flow architecture.
- [Physical calculation flow](physics.md) - end-to-end description of the active LAGRID plume-model path, from inputs and EPM through transport, microphysics, meteorology, remapping, diagnostics, and chemistry.
- [Numerical methods](numerics.md) - summary of the discretizations, solvers, integration methods, remapping, interpolation, and conservation/stability points used by the physical workflow.

## Subsystem notes

- [Core](CORE.md) - top-level application layer: `Input`, `OptInput`, `LAGRIDPlumeModel`, meteorology, diagnostics, aircraft/fuel/emission objects, timestep wrappers, and statuses.
- [EPM](EPM.md) - Early Plume Model interfaces and implementations, including the original built-in model and external NetCDF EPM input path.
- [AIM](AIM.md) - aerosol and ice microphysics: size distributions, coagulation, growth/sublimation, nucleation, settling, and diagnostics.
- [FVM_ANDS](FVM_ANDS.md) - finite-volume advection/diffusion transport solver, boundary handling, operator splitting, and sparse solve structure.
- [LAGRID](LAGRID.md) - adaptive/free-coordinate remapping, mass-box representation, overlap mapping, grid initialization, and buffers.
- [KPP](KPP.md) - KPP-generated chemistry, stiff forward/adjoint integrators, generated globals, and thread-safety considerations.
- [YamlInputReader](YAML_INPUT_READER.md) - YAML/default input loading, parser helpers, parameter sweeps, Monte Carlo expansion, and configuration flow.
- [Util](UTIL.md) - shared physical functions, meteorological helpers, vector/mask utilities, YAML helpers, random sampling, and error utilities.

## Issue and risk audits

- [Numerical issues](numerics-issues.md) - potential convergence and accuracy risks in EPM integration, aerosol moments, coagulation, settling, transport, SOR diffusion solves, ice growth, met interpolation, pressure-coordinate conversion, LAGRID remapping, diagnostics, and chemistry.
- [Code issues](code-issues.md) - broader robustness risks: bad parameter values, indexing/off-by-one hazards, OpenMP/shared-state issues, file/output handling, error handling, memory/resource concerns, and unit/data-model risks.

## Suggested reading order

1. [Architecture](ARCHITECTURE.md)
2. [Physical calculation flow](physics.md)
3. [Numerical methods](numerics.md)
4. Relevant subsystem note(s) from the list above
5. [Numerical issues](numerics-issues.md) and [Code issues](code-issues.md) before modifying solver, remapping, parallelism, or input-handling code
