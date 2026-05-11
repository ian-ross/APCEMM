# Potential numerical convergence and accuracy issues

This document reviews the methods summarized in [Numerical methods](numerics.md), in the context of the active plume calculation described in [Physical calculation flow](physics.md). Items below are potential issues to verify with convergence tests, conservation checks, and case-specific sensitivity studies; not all are necessarily active failures for the default configuration.

## 1. Case expansion, initialization, and grids

- **Grid resolution is fixed by input or LAGRID heuristics, not by error estimates.** Initial Cartesian spacing and LAGRID remap spacing are not selected from truncation-error indicators, so plume gradients, shear layers, or narrow humidity structures can be under-resolved.
- **Log-radius bin moments use midpoint quadrature.** For broad or sharply evolving size distributions, moments such as volume, area, effective radius, and extinction can have bin-resolution bias.
- **Parameter sweeps and Monte Carlo cases are independent but not a convergence study.** They vary inputs, not numerical parameters, so they cannot by themselves establish time-step, grid, or bin convergence.

## 2. EPM near-field integration

- **Explicit RKF78 may struggle with stiffness.** EPM combines dilution, condensation/deposition, nucleation, and coagulation-like terms. If nucleation or phase-change rates become stiff, an explicit adaptive method can take very small steps, fail, or produce tolerance-dependent results.
- **Physics is split at coarse logarithmic output intervals.** Sulfate self-coagulation is applied once before each EPM log interval, and newly nucleated sulfate is inserted after each interval. This operator splitting can make results sensitive to the 301-point logarithmic output grid even if the ODE tolerance is tight.
- **Discontinuous or piecewise RHS functions can reduce order.** Piecewise entrainment and temperature-ramp formulas introduce derivative discontinuities. Adaptive RK methods remain usable, but high formal order is not realized near breakpoints.
- **Positivity is partly enforced after the fact.** Several rate evaluations use clipping such as `max(..., 0)`. This avoids nonphysical calls but can hide tolerance/step-size problems and introduce nonsmooth dynamics.

## 3. Aerosol moments and diagnostics

- **Moment accuracy depends strongly on bin placement.** Midpoint log-bin quadrature can under-represent distribution tails, especially for high-order moments used in area/volume/extinction.
- **Diagnostic reductions inherit transport/remap errors.** Optical depth, line totals, and effective radii combine several approximate operations: size-bin quadrature, cell-area integration, transport diffusion, remapping, and pressure-coordinate unit conversions.
- **Ratios of moments can be ill-conditioned.** Effective radius and standard-deviation-like quantities can become noisy in cells with very small particle number or near complete sublimation.

## 4. Coagulation

- **Semi-implicit loss treatment is stable but only first-order in time.** The update `(old + dt*P)/(1 + dt*L)` can be too diffusive or biased for large coagulation timesteps.
- **Production terms use old-state information.** Strong coagulation over one step may require subcycling or timestep sensitivity checks, because only the loss is implicit.
- **Destination-bin mapping is discrete.** Volume-conserving bin fractions reduce mass error, but number and shape of the size distribution can still be grid-dependent.

## 5. Settling velocities

- **Fall speeds are computed from bin centers and reference met state.** If particle centers move within bins or density/met conditions vary strongly, fixed bin-center velocities can bias vertical transport.
- **The active Stokes/slip formula may be outside its range for large crystals.** The alternative Reynolds-number-dependent branch is disabled, so large-particle settling should be checked against expected Reynolds numbers.
- **Settling is operator-split from growth/remapping.** Growth changes particle size before/after transport, but transport uses discrete bin velocities, which can lag rapidly changing crystals.

## 6. Finite-volume transport

- **Explicit advection has only a CFL stability control.** The code subcycles by Courant number, but accuracy still depends on the selected CFL and global model timestep.
- **Zero or tiny velocity may be a corner case.** The advection substep calculation divides by the current Courant number in `FVM_Solver::operatorSplitSolve`; cases with exactly zero velocity should be verified.
- **Strang splitting is only second-order when operators are smooth and accurately solved.** Boundary conditions, nonlinear coefficient updates, remapping, and source/growth splitting can reduce the practical order.
- **The TVD limiter reduces to first order near extrema and boundaries.** This prevents oscillations but adds numerical diffusion and may broaden contrail edges.
- **Boundary conditions can dominate budgets.** Dirichlet ghost-cell boundaries are supported in the active path; outflowing plume material or humidity anomalies can be artificially lost or constrained when the domain is too small.
- **A possible limiter implementation hazard should be reviewed.** In `AdvDiffSystem.hpp`, one south-face limiter branch appears to test `neighbor_point(...)` directly in a boolean expression, which may force the limiter to zero for many valid points and increase numerical diffusion.

## 7. Implicit diffusion and SOR convergence

- **The SOR solve has no maximum iteration cap.** `sor_solve` iterates until a residual threshold is met; if the system converges slowly or stalls, the model can hang.
- **Relative residual can be problematic for small right-hand sides.** The residual is divided by `||rhs||_2`; nearly zero RHS values can produce noisy, infinite, or NaN residuals.
- **Default SOR parameters are conservative and may be slow.** The default relaxation factor is `omega = 1.0` and only a few sweeps are done between residual checks. Large grids or strong diffusion can require many iterations.
- **Matrix reuse assumes unchanged coefficients and grid.** Sharing a diffusion matrix is efficient, but any mismatch in grid, timestep, boundary conditions, or diffusion coefficients would make the solve inconsistent.
- **Solver tolerance is not clearly tied to physical error.** A relative linear residual threshold does not directly bound concentration, mass, or optical-depth error.

## 8. Ice growth, sublimation, and size remapping

- **APC growth is stable but not fully error-controlled.** The vapor-coupled analytical update avoids explicit instability, but the model still needs timestep checks for rapidly changing supersaturation.
- **Clipping can violate smooth convergence.** Ice volume is clipped to nonnegative values and to maximum per-particle bin volume. This preserves physical bounds but can introduce timestep-dependent mass/number changes.
- **Particles below the minimum bin are removed.** Sublimating crystals that fall below the smallest volume edge disappear; results may depend on the minimum radius and bin resolution.
- **Particles above the maximum bin are retained in the largest bin.** The code comments flag this as a TODO. Large crystals can pile up in the last bin, biasing fall speed, effective radius, and extinction.
- **Moving-center remapping is conservative for selected quantities, not shape-preserving.** Mapping all particles from an old bin into one destination bin can introduce numerical diffusion or artificial narrowing in size space.

## 9. Meteorology interpolation and hydrostatic grid motion

- **Linear/nearest interpolation can miss sharp vertical structure.** Humidity, temperature, and shear inversions are especially sensitive to vertical met resolution.
- **Hydrostatic mapping assumes piecewise-constant lapse rates.** Accuracy depends on met-level spacing and on the validity of the lapse-rate approximation between levels.
- **Pressure-coordinate displacement is operator-split.** Applying `omega * dt` to pressure edges separately from transport/growth can cause timestep-dependent vertical placement and humidity exposure.
- **Random temperature perturbations are grid-cell noise.** Cellwise uncorrelated perturbations can add grid-scale variability rather than physically correlated turbulence; convergence with grid refinement may be nonstandard.

## 10. Pressure-coordinate unit conversion

- **The local number-density proxy can amplify errors.** Dividing fields by `(p_edge[j] - p_edge[j+1])/(y_edge[j+1] - y_edge[j])` is sensitive to small layer thicknesses or pressure-edge inconsistencies.
- **Unit-conversion/remap conservation should be tested directly.** Because fields are converted to mixing-ratio-like quantities before grid motion, integrated mass conservation depends on the consistency of the pressure-edge update, altitude inversion, and subsequent remapping.

## 11. LAGRID adaptive remapping

- **The area-overlap remap is first order.** Piecewise-constant old cells are conservative but diffusive; repeated remapping can broaden gradients and reduce peaks.
- **Mass can be lost if the active mask excludes relevant tails.** Conservation holds for retained boxes. Low-concentration ice, tracer, or H2O anomaly outside the mask/buffer may be discarded or filled from defaults.
- **Grid-spacing heuristics can cause abrupt resolution changes.** Hard limits and mask-size formulas may create timestep-to-timestep changes in truncation error.
- **No monotone high-order reconstruction is used.** This avoids overshoot but limits remap accuracy for smooth fields.

## 12. Chemistry numerics

- **KPP Rosenbrock tolerances affect coupled initialization.** Stiff chemistry is handled by appropriate integrators, but tolerances and Jacobian updates should be checked when chemistry materially affects EPM/background state.
- **Threadprivate/global state requires care in parallel cases.** Parallel case execution can expose hidden coupling if generated chemistry state is not fully isolated.
- **Coupling cadence matters if resolved-plume chemistry is re-enabled.** Transport, heterogeneous rates, and chemistry would need splitting-error and tolerance checks.

## 13. Recommended verification tests

- Run **time-step convergence** for EPM interval count, plume transport timestep, ice-growth timestep, and pressure-grid update cadence.
- Run **grid convergence** for initial grid spacing and LAGRID remap limits, including cases with strong shear and sharp humidity gradients.
- Run **size-bin convergence** for sulfate and ice moments, settling, growth/sublimation, and extinction.
- Track **global budgets** for H2O, ice mass, particle number, sulfate mass, and passive tracer before/after transport, growth, pressure conversion, and remapping.
- Add **linear-solver diagnostics**: SOR iteration counts, final residuals, zero-RHS handling, and a maximum iteration failure mode.
- Compare selected cases against **alternative solver settings**: tighter EPM/KPP tolerances, smaller CFL, disabled/remap-less fixed-grid runs, and finer met interpolation grids.
