# Core subsystem

`Core` is the top-level application layer. It turns parsed inputs into case objects, constructs aircraft/emission/meteorology state, runs the active plume simulation, and writes diagnostics/status.

## Main code

- `include/Core/Input.hpp`, `src/Core/Input.cpp` - per-case scalar inputs and output filenames.
- `include/Core/Input_Mod.hpp`, `src/Core/Input_Mod.cpp` - global `OptInput` options populated from YAML/defaults.
- `include/Core/LAGRIDPlumeModel.hpp`, `src/Core/LAGRIDPlumeModel.cpp` - primary case orchestrator for the current executable path.
- `Aircraft`, `Engine`, `Fuel`, `Emission` - aircraft/fuel/emission-index calculations.
- `Meteorology` - ambient profiles, met-file/time-series loading, interpolation to the plume grid, perturbation updates.
- `Diag_Mod` - NetCDF diagnostic output.
- `MPMSimVarsWrapper`, `TimestepVarsWrapper` - option-to-runtime flag/timestep translation.
- `Status` - normal and early-stop/failure result codes.

## LAGRIDPlumeModel responsibilities

`LAGRIDPlumeModel` owns the coupled plume state for one case:

- immutable `OptInput` and mutable `Input` references;
- aircraft/fuel/emission objects;
- timestep and model-switch wrappers;
- `Meteorology` profiles/fields;
- `AIM::Grid_Aerosol` for ice size distribution on the 2-D grid;
- coordinate vectors/edges for the current adaptive grid;
- water vapor (`H2O_`), contrail tracer, diffusion coefficients, fall speeds, and diagnostic scalars.

`runFullModel()` follows this sequence:

1. Set the per-case OpenMP thread count.
2. Run EPM through `runEPM()` to get near-field plume initial conditions or an early `SimStatus`.
3. Initialize a 2-D grid and map EPM aerosol/water/tracer quantities into it.
4. Initialize gravitational settling velocities when enabled.
5. Iterate until the requested simulation time or an early stop:
   - schedule transport, ice growth, meteorology, remapping, and output using `TimestepVarsWrapper`;
   - run finite-volume transport for H2O, contrail tracer, and each aerosol bin;
   - update temperature perturbations and met fields;
   - call `Grid_Aerosol::Grow` for ice deposition/sublimation;
   - use masks of ice/H2O/tracer content to decide retained cells;
   - remap all retained fields to a new LAGRID grid;
   - save NetCDF diagnostics/time-series as requested.

## Meteorology

`Meteorology` can either construct ambient fields from scalar inputs/default lapse-rate assumptions or read NetCDF met fields/time series. It maintains vertical-coordinate-dependent pressure, altitude, shear, vertical velocity, and 2-D fields of temperature, water vapor, and air molecular density. It can:

- interpolate met profiles to the simulation `y` grid;
- apply prescribed updraft displacement;
- update diurnal/turbulent temperature perturbations;
- regenerate fields after LAGRID changes the grid;
- expose reference-altitude values used by EPM.

Humidity conversions rely on `Util/PhysFunction` saturation vapor-pressure relationships. Pressure is treated as vertically varying only.

## Physics and coupling

Core itself is mostly orchestration, but it defines important coupling assumptions:

- EPM provides the near-field plume area, temperature perturbation, water vapor, and initial aerosol PDFs.
- The later plume is modeled as a 2-D cross-section in horizontal/vertical coordinates, moving with the aircraft track.
- Ambient shear and vertical velocity enter transport/advection and grid-relative plume displacement.
- Ice mass and water vapor are coupled through AIM growth; temperature and pressure fields come from `Meteorology`.
- Diagnostics depend on derived aerosol moments such as ice water content, optical depth, extinction, effective radius, and total number.

## Integration notes

- `LAGRIDPlumeModel` is the highest-coupling class in the codebase; changes to grids, aerosol array shapes, met update frequency, or transport vectorization usually touch it.
- Most exchanged arrays are `Vector_1D`, `Vector_2D`, or `Vector_3D` typedefs from `Util/ForwardDecl.hpp`.
- Diagnostics and outputs are tied to paths and flags in `Input`/`OptInput`.
