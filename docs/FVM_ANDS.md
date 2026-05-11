# FVM_ANDS subsystem

`FVM_ANDS` is the finite-volume advection/diffusion solver used to transport water vapor, contrail tracer, and aerosol-bin fields on the current structured plume grid.

## Main code

- `include/FVM_ANDS/AdvDiffSystem.hpp`, `src/FVM_ANDS/AdvDiffSystem.cpp` - discrete equation assembly and RHS construction.
- `include/FVM_ANDS/FVM_Solver.hpp`, `src/FVM_ANDS/FVM_Solver.cpp` - solver wrapper and update API.
- `include/FVM_ANDS/BoundaryCondition.hpp`, `src/FVM_ANDS/BoundaryCondition.cpp` - boundary/ghost-point representation.
- `include/FVM_ANDS/FVM_ANDS_HelperFunctions.hpp` - vector-index and conversion helpers.
- `include/FVM_ANDS/FVM_ANDS_enums.hpp` - directions, boundary flags, vector format.

## Numerical model

The subsystem discretizes a 2-D advection-diffusion equation for scalar concentration-like fields:

```text
∂phi/∂t + u ∂phi/∂x + v ∂phi/∂y = ∂/∂x(Dh ∂phi/∂x) + ∂/∂y(Dv ∂phi/∂y) + source
```

`AdvDiffParams` supplies base horizontal velocity, vertical velocity, shear, horizontal/vertical diffusion, and timestep. Diffusion can be scalar, `Vector_2D`, or Eigen-vector fields.

## Discretization and solver

`AdvDiffSystem`:

- builds interior and ghost `Point` objects;
- applies top/left/right/bottom boundary conditions;
- assembles an Eigen row-major sparse coefficient matrix;
- computes RHS/source vectors;
- constructs velocity vectors including shear effects;
- provides forward-Euler advection and implicit/operator-split paths.

Advection uses a minmod flux limiter for total-variation-diminishing behavior. `FVM_Solver` solves implicit systems with Eigen `BiCGSTAB` and an optional diagonal preconditioner. It also exposes explicit and operator-split solves, half-timestep advection, Courant-limited splitting, and matrix reuse through shared sparse-matrix pointers.

## Coupling to plume model

`LAGRIDPlumeModel::runTransport()` creates/updates solvers for plume fields, sets boundary conditions, updates advection/diffusion/timestep values, and converts between APCEMM `Vector_2D`/`Vector_3D` arrays and Eigen vectors. The same matrix can be reused for fields with identical grid/transport parameters.

## Implementation notes

- Grid changes require updating spacing, `nx`, `ny`, field vector sizes, and boundary-condition vectors together.
- The solver contains both physical boundary points and ghost points; vector lengths are larger than `nx * ny`.
- Courant number and flux-limiter behavior are central to stable advection changes.
