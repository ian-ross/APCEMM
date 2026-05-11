# APCEMM Architecture

This document is an initial high-level architecture map for the C++ source tree in `Code.v05-00/`. APCEMM is the Aircraft Plume Chemistry, Emissions, and Microphysics Model. The current executable simulates an aircraft exhaust plume/contrail using an early plume model followed by a 2-D LAGRID plume simulation with transport, meteorology, and ice microphysics.

## Source tree overview

- `Code.v05-00/src/Main.cpp` - CLI entry point and case driver.
- `Code.v05-00/include/` - public/internal headers grouped by subsystem.
- `Code.v05-00/src/Core/` - top-level domain objects, simulation orchestration, meteorology, diagnostics, aircraft/fuel/emission models, timestep and simulation option wrappers.
- `Code.v05-00/src/EPM/` - Early Plume Model implementations and serialization of EPM output.
- `Code.v05-00/src/AIM/` - aerosol/ice microphysics: aerosol size distributions, growth, coagulation, nucleation, settling.
- `Code.v05-00/src/FVM_ANDS/` - finite-volume advection/diffusion numerical solver built on Eigen sparse matrices.
- `Code.v05-00/src/LAGRID/` - adaptive/free-coordinate grid remapping utilities used by the plume model.
- `Code.v05-00/src/KPP/` - KPP-generated chemistry and adjoint integration code.
- `Code.v05-00/src/YamlInputReader/` - YAML configuration parsing and parameter-sweep expansion.
- `Code.v05-00/src/Util/` - shared physical/math/vector/YAML utilities and random-number seeding.
- `Code.v05-00/defaults/` - default input data embedded into generated headers at configure time.
- `Code.v05-00/tests/` - Catch2-based unit/integration tests for major numerical and input components.

## Build and packaging

The project is built with CMake (`Code.v05-00/CMakeLists.txt`) as a C++20 application named `APCEMM`.

The build creates static libraries for each subsystem:

```text
APCEMM executable
  links: FVM_ANDS, Util, AIM, KPP, EPM, Core, YamlInputReader, LAGRID
```

External dependencies are resolved through vcpkg and CMake `find_package` calls:

- Boost headers
- Catch2
- Eigen3
- fmt
- FFTW3
- netCDF C/C++
- yaml-cpp
- OpenSSL
- OpenMP

At configure time, CMake also generates includable headers from default text/YAML files:

- `defaults/input.yaml` -> `include/Defaults/Input.hpp`
- `defaults/Ambient.txt` -> `include/Defaults/Ambient.hpp`
- `defaults/ENG_EI.txt` -> `include/Defaults/Engine_EI.hpp`

## Runtime control flow

The main simulation path is:

```text
main(argc, argv)
  -> YamlInputReader::readYamlInputFiles(...)
  -> setSeed(...)
  -> YamlInputReader::generateCases(...)
  -> for each case, optionally in OpenMP parallel loop:
       -> construct Core::Input case object
       -> construct LAGRIDPlumeModel
       -> LAGRIDPlumeModel::runFullModel()
       -> write status_caseN file
```

`Main.cpp` currently hardcodes `model = 1`, so the active path is the plume model. Box model and adjoint-only branches are present as placeholders but not implemented in the entry-point switch.

Outputs are written under `OptInput::SIMULATION_OUTPUT_FOLDER`. Forward/adjoin/box/microphysical output filenames are assembled per case, usually including a zero-padded case number. Diagnostic fields are saved as NetCDF files through `Core/Diag_Mod`.

## Main simulation orchestrator

`Core/LAGRIDPlumeModel` is the central runtime coordinator. It owns or references most state needed for a case:

- immutable run options: `OptInput`
- case-specific input: `Input`
- aircraft/fuel/emissions: `Aircraft`, `Fuel`, `Emission`
- simulation/timestep wrappers: `MPMSimVarsWrapper`, `TimestepVarsWrapper`
- meteorology: `Meteorology`
- aerosol state: `AIM::Grid_Aerosol`
- grid coordinates and edges
- water vapor, contrail tracer, diffusion fields, settling velocities

High-level `runFullModel()` behavior:

1. Configure OpenMP thread count.
2. Run the EPM (`runEPM()`), returning either an `EPM::Output` or an early `SimStatus`.
3. Initialize the 2-D simulation grid and H2O field from EPM output.
4. Initialize settling velocities if gravitational settling is enabled.
5. Enter the time loop:
   - run transport with `FVM_ANDS` when scheduled;
   - update temperature perturbations and meteorology;
   - run ice growth with `AIM::Grid_Aerosol::Grow` when scheduled;
   - update contrail tracer/masks;
   - remap/regrid retained fields with `LAGRID` utilities;
   - save diagnostics/time-series outputs as configured;
   - stop early when no contrail-containing cells remain or other status criteria apply.

## Subsystems

### Core

Detailed subsystem documentation: [Core](CORE.md).

`Core` contains the domain model and high-level glue code:

- `Input` represents one generated case, including location/time, aircraft, emissions, background species, and output filenames.
- `Input_Mod`/`OptInput` holds global options parsed from YAML.
- `Aircraft`, `Engine`, `Fuel`, and `Emission` compute aircraft and exhaust properties.
- `Meteorology` constructs or reads meteorological profiles/time series, interpolates them to the simulation grid, and updates temperature, humidity, pressure, shear, vertical velocity, and air density.
- `Diag_Mod` writes NetCDF diagnostics.
- `MPMSimVarsWrapper` and `TimestepVarsWrapper` translate input options into simulation flags and timestep schedules.
- `Status` defines case completion/failure status values.

### YamlInputReader

Detailed subsystem documentation: [YamlInputReader](YAML_INPUT_READER.md).

`YamlInputReader` parses one or more YAML input files into `OptInput`, performs validity checks, and expands parameter sweeps/Monte Carlo inputs into per-case parameter maps consumed by `Input`.

### EPM

Detailed subsystem documentation: [EPM](EPM.md).

The Early Plume Model computes near-field plume evolution and initial conditions for the later 2-D plume simulation. The public interface is `EPM::make_epm(...)` plus `EPM::Output`, which stores initial temperature, ice/soot/aerosol quantities, water vapor, sulfur species, plume area, and engine geometry fields. Model implementations live under `src/EPM/Models/`.

### AIM

Detailed subsystem documentation: [AIM](AIM.md).

AIM is the aerosol/ice microphysics package. It provides:

- `AIM::Aerosol` for one-dimensional size distributions.
- `AIM::Grid_Aerosol` for binned aerosol fields over the 2-D plume grid.
- coagulation kernels, nucleation, settling velocity, and ice growth routines.

The plume model primarily uses `Grid_Aerosol` for ice particle number/volume/PDF fields and derived diagnostics such as extinction, optical depth, radius, effective radius, and ice water content.

### FVM_ANDS

Detailed subsystem documentation: [FVM_ANDS](FVM_ANDS.md).

`FVM_ANDS` implements advection/diffusion transport. `FVM_Solver` wraps `AdvDiffSystem`, boundary conditions, Eigen sparse matrices, and iterative BiCGSTAB solution. It supports implicit/explicit solves and operator splitting, and is used by `LAGRIDPlumeModel::runTransport()`.

### LAGRID

Detailed subsystem documentation: [LAGRID](LAGRID.md).

`LAGRID` supplies adaptive/free-coordinate grid remapping. It represents active-mass boxes (`FreeCoordBoxGrid`, `MassBox`) and remaps fields after transport/meteorological updates so that the domain can track the contrail region rather than keeping a fixed full grid.

### KPP

Detailed subsystem documentation: [KPP](KPP.md).

`KPP` contains generated chemical mechanism code and C-compatible integration entry points (`INTEGRATE`, `INTEGRATE_ADJ`, rate updates, photolysis, heterogeneous chemistry). In the current LAGRID plume path, KPP is still required for background ambient spin-up during EPM initialization. Some KPP state is global (`RCONST`, `PHOTOL`, `HET`, etc.), which affects thread-safety considerations.

### Util

Detailed subsystem documentation: [Util](UTIL.md).

`Util` contains shared helper code for physical constants/functions, meteorological interpolation, plume-model utility functions, vector utilities/masking, YAML helpers, errors, and random-number generation.

## Physical calculation flow

A detailed end-to-end trace of the model physics is available in [Physical calculation flow](physics.md).
Numerical methods used by those physical processes are summarized in [Numerical methods](numerics.md).

## Data flow

```text
YAML input files + embedded defaults
  -> OptInput
  -> generated case parameter maps
  -> Input per case
  -> Aircraft/Fuel/Emission + Meteorology
  -> EPM::Output initial plume state
  -> LAGRIDPlumeModel grid fields:
       ice aerosol PDF, H2O, contrail tracer, met fields
  -> repeated transport / growth / met update / remap
  -> NetCDF diagnostics and status files
```

Meteorological data can come from user/default input values or NetCDF met files/time series. Simulation diagnostics and outputs are written with netCDF C++ APIs.

## Parallelism and numerical libraries

OpenMP is used in the case loop and within numerical kernels. The case loop in `Main.cpp` uses dynamic scheduling when `PARALLEL_CASES` is enabled. Per-case `LAGRIDPlumeModel` also sets the OpenMP thread count from input options.

Eigen is used for sparse matrices and iterative solves in transport and regridding. FFTW3 is linked as a dependency, although its exact use should be verified in lower-level numerical code before making changes that affect linkage.

## Testing

Tests are configured from `Code.v05-00/tests/CMakeLists.txt` and use Catch2. Existing tests cover areas such as:

- aerosol and microphysics routines
- advective/diffusive solver behavior
- LAGRID/remapping/grid functionality
- meteorology and physical functions
- nucleation and kernel construction
- YAML input parsing
- aircraft logic

## Architectural notes and constraints

- `LAGRIDPlumeModel` is the key integration point and currently carries substantial state and orchestration logic.
- The active executable path is plume-model-only (`model = 1` in `Main.cpp`); other model modes are stubs in the CLI switch.
- KPP-generated code exposes global variables required by chemistry routines. This should be considered when changing parallel execution or adding new concurrent chemistry use.
- Many shared data structures are typedef-style vectors (`Vector_1D`, `Vector_2D`, `Vector_3D`) from `Util/ForwardDecl.hpp`; modules exchange these directly.
- Output and diagnostics are tightly coupled to netCDF and filesystem paths held in `Input`/`OptInput`.
