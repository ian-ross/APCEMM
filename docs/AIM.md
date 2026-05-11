# AIM subsystem

AIM is APCEMM's aerosol and ice microphysics package. It represents particle size distributions, evolves them by coagulation and growth, and computes optical/mass diagnostics.

## Main code

- `include/AIM/Aerosol.hpp`, `src/AIM/Aerosol.cpp` - 1-D and grid aerosol distributions.
- `include/AIM/Coagulation.hpp`, `src/AIM/Coagulation.cpp` - coagulation kernel object and application.
- `include/AIM/buildKernel.hpp`, `src/AIM/buildKernel.cpp` - kernel construction.
- `include/AIM/Nucleation.hpp`, `src/AIM/Nucleation.cpp` - sulfuric-acid/water nucleation parameterizations.
- `include/AIM/Settling.hpp`, `src/AIM/Settling.cpp` - ice-particle terminal fall velocities.

## Data structures

`AIM::Aerosol` represents one size distribution with bin centers, bin edges, bin sizes, volume centers, and a PDF. The PDF convention is `dN / d(ln r)` in `#/cm^3`-style concentration units.

`AIM::Grid_Aerosol` extends this to the plume grid:

```text
pdf[y][x][bin]
bin_VCenters[y][x][bin]
```

It stores the number PDF and per-cell bin volume centers so transported/grown particles can keep cell-local representative sizes.

## Microphysical processes

AIM provides:

- **Moments and diagnostics**: number, volume, area, radius, effective radius, standard deviation, ice water content, extinction, optical depth, and integrated totals.
- **Coagulation**: bin-to-bin particle collisions using precomputed kernels and remapping of collision products into size bins.
- **Ice growth/sublimation**: `Grid_Aerosol::Grow` couples ice volume changes to vapor-field changes. It uses effective vapor diffusivity, thermal/latent-heat corrections, and bin-particle flux between bins.
- **Nucleation**: helper functions for cluster composition, nucleation rate, total cluster number, cluster radius, and threshold criteria.
- **Settling**: `SettlingVelocity` computes terminal fall speeds by particle radius using physical functions for air viscosity, slip correction, and spherical-particle fall speed.

## Coupling to plume model

`LAGRIDPlumeModel` primarily uses `Grid_Aerosol` for ice. During each active timestep it may:

1. transport every bin number and volume field with FVM_ANDS;
2. call `UpdateCenters` after transport so volume/number stay consistent;
3. call `Grow(dt, H2O, T, P)` to exchange mass with vapor;
4. derive masks and diagnostics from `TotalNumber`, `IWC`, `Extinction`, `EffRadius`, etc.;
5. remap PDFs when LAGRID changes the structured grid.

## Implementation notes

- Several routines support symmetry flags and reduced-domain calculations; check `CheckCoagAndGrowInputs` before changing grid dimensions.
- Units are mixed by historical convention: radii are meters, PDFs are effectively per `cm^3`, vapor is often `molec/cm^3`, and diagnostics convert as needed. Preserve unit conversions carefully.
