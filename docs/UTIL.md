# Util subsystem

`Util` contains shared math, physics, meteorological, vector, YAML, error, and random-number helpers used throughout APCEMM.

## Main code

- `include/Util/ForwardDecl.hpp` - common aliases such as `Vector_1D`, `Vector_2D`, `Vector_3D`, `UInt`.
- `PhysConstant.hpp` - physical constants.
- `PhysFunction.hpp`, `src/Util/PhysFunction.cpp` - saturation vapor pressures, air/particle properties, diffusion coefficients, humidity conversions.
- `MetFunction.hpp`, `src/Util/MetFunction.cpp` - meteorological interpolation/helpers.
- `PlumeModelUtils.hpp`, `src/Util/PlumeModelUtils.cpp` - timestep, advection, and diffusion-parameter helpers for plume dynamics.
- `VectorUtils.hpp`, `src/Util/VectorUtils.cpp` - vector/mask operations and field reductions.
- `YamlUtils.hpp`, `src/Util/YamlUtils.cpp` - YAML utility functions.
- `MC_Rand.hpp`, `src/Util/MC_Rand.cpp` - random seeding/sampling support.
- `Error.hpp`, `src/Util/Error.cpp` - shared error handling.

## Physical utilities

`physFunc` centralizes formulas used by EPM, AIM, Core meteorology, and transport setup:

- saturation vapor pressures over liquid water, ice, sulfuric acid, and nitric acid;
- derivatives of saturation pressure;
- air density, viscosity, kinematic viscosity, mean free path;
- particle mass, terminal fall speed, Knudsen/Reynolds/Schmidt/Stokes numbers;
- slip-flow correction, particle diffusion, gas diffusion coefficients;
- latent heat of sublimation, thermal conductivity, Kelvin factor;
- conversions among H2O concentration, RH over water, and RH over ice.

Using these shared formulas avoids inconsistent thermodynamic assumptions across subsystems.

## Plume/vector utilities

`PlumeModelUtils` provides scheduling and simple plume kinematics:

- `UpdateTime` and `BuildTime` for dynamic timestep schedules;
- `AdvGlobal` for updraft-driven displacement/velocity;
- `DiffParam` for time-dependent diffusion parameters.

`VectorUtils` is heavily used by `LAGRIDPlumeModel` to compute masks, extrema, and reductions over `Vector_2D` fields. Masks carry both selected-cell indices and bounding information used during remapping.

## Implementation notes

- Unit conventions come from callers; utility functions document expected units in headers where possible.
- Because vector typedefs are ubiquitous, changes to `ForwardDecl.hpp` ripple through nearly every subsystem.
- Keep utility code dependency-light; it is linked by most static libraries.
