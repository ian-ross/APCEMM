# EPM subsystem

The Early Plume Model (EPM) computes near-field aircraft-exhaust evolution and produces the initial conditions for the 2-D LAGRID plume model.

## Main code

- `include/EPM/EPM.hpp`, `src/EPM/EPM.cpp` - public factory and serializable output.
- `include/EPM/Models/Base.hpp` - common model interface.
- `include/EPM/Models/Original.hpp` and `src/EPM/Models/Original/` - built-in original EPM implementation.
- `include/EPM/Models/External.hpp`, `src/EPM/Models/External.cpp` - read initial conditions from an external NetCDF EPM file.
- `include/EPM/Utils/Physics.hpp`, `src/EPM/Utils/Physics.cpp` - EPM-specific entrainment, condensation/deposition/freezing helpers.

`EPM::make_epm(...)` chooses the implementation from `OptInput::SIMULATION_EPM_TYPE`. `EPM_NEW_PHYSICS` is declared but not implemented in the factory.

## Output contract

`EPM::Output` is the boundary between EPM and the later plume model. It contains:

- final plume temperature and water/sulfur gas concentrations;
- sulfuric-acid and ice `AIM::Aerosol` PDFs;
- representative ice radius/density and soot density;
- plume area, bypass area, and core exit temperature.

`Output::write` and `Output::read` serialize these fields to NetCDF. Scalar fields are stored individually; sulfur and ice PDFs use their own radius-center/radius-edge dimensions.

## Physics represented

EPM covers the early jet/vortex/dispersion phase that cannot be resolved by the later 2-D grid. It estimates dilution and thermodynamic evolution of fresh exhaust through:

- entrainment/mixing of exhaust with ambient air;
- plume temperature evolution, including vortex-regime parameterizations;
- water vapor condensation/deposition and possible freezing;
- sulfur partitioning into gas/liquid aerosol;
- soot/ice particle initialization;
- ambient chemistry spin-up through KPP where required.

The original implementation has ODE/RHS components for gas-aerosol evolution and microphysics routines that decide whether a contrail survives into the resolved plume stage. If EPM determines non-survival or invalid initial conditions, it can return `SimStatus` instead of `Output`.

## Integration notes

- EPM depends on `Core` aircraft/emission/meteorology objects but should not own long-lived plume-grid state.
- The plume model assumes EPM output units exactly as documented in `EPM::Output` (`molec/cm^3`, `m`, `kg/m^3`, `m^2`).
- External EPM files must provide the same NetCDF variables used by `Output::read`.
