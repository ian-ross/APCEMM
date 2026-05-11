# LAGRID subsystem

`LAGRID` provides adaptive/free-coordinate remapping so APCEMM can keep resolution around the contrail instead of simulating a large fixed domain.

## Main code

- `include/LAGRID/FreeCoordBoxGrid.hpp`, `src/LAGRID/FreeCoordBoxGrid.cpp` - mass-box representation of active cells.
- `include/LAGRID/RemappingFunctions.hpp`, `src/LAGRID/RemappingFunctions.cpp` - structured-grid remapping and initialization helpers.

## Core concepts

`MassBox` represents the mass of a scalar field in an arbitrary rectangular box:

```text
mass = phi * dx * dy
```

It stores top-left and bottom-right coordinates plus mass. `FreeCoordBoxGrid` is a collection of these boxes, usually created from a structured grid and a mask selecting active cells. It also records boundary boxes and min/max coordinates to accelerate remapping.

`twoDGridVariable` bundles a `Vector_2D` field with x/y coordinates and inferred uniform `dx`, `dy`. Buffers can be added around a remapped variable with a fill value.

## Remapping workflow

Typical plume-model use:

1. Build a mask from ice number, water-vapor perturbation, or contrail tracer.
2. Convert retained structured-grid cells into a `FreeCoordBoxGrid`.
3. Choose a target `Remapping` with origin, spacing, and dimensions.
4. Compute overlap area between old mass boxes and new grid cells with `coveredArea`.
5. Map mass back to concentration-like `phi` on the new structured grid.
6. Add buffers and regenerate meteorology/transport fields on the new coordinates.

The overlap-area approach is designed to conserve mass for concentration fields when changing grids, subject to mask and buffer choices.

## Initialization helpers

LAGRID also initializes fields from idealized distributions:

- rectangular plume distributions;
- Gaussian distributions;
- bimodal vertical distributions;
- generic weighted distribution via a user-supplied weight function.

These are used when placing EPM output mass/number into the initial 2-D plume grid.

## Coupling and constraints

`LAGRIDPlumeModel::remapAllVars()` remaps H2O, contrail tracer, aerosol PDFs/volumes, meteorology-dependent fields, and diffusion arrays. Any field added to the plume state must either be remapped there or explicitly regenerated from met/input data after grid changes.
