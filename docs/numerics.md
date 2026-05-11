# Numerical methods

This document summarizes the numerical methods used for the physical processes described in [Physical calculation flow](physics.md). It is based on the current active LAGRID plume path and the code inspected in `Code.v05-00/src` and `Code.v05-00/include`.

## 1. Case expansion and initialization numerics

Input parsing is not a physical solver, but it determines the numerical experiment. Parameter sweeps and Monte Carlo inputs are expanded into case maps by `YamlInputReader::generateCases`; each case is then run independently. The main executable can parallelize cases with OpenMP dynamic scheduling.

Initial grids are structured Cartesian cross-sections. `LAGRIDPlumeModel::initializeGrid()` builds uniformly spaced `x` and `y` edge/center arrays from `ADV_GRID_*` options. Aerosol bins are logarithmic in radius/volume in EPM/AIM code, with moments evaluated by summing over log-radius bin widths.

## 2. EPM near-field integration

### Time grid and ODE solver

The original EPM model (`src/EPM/Models/Original/Integrate.cpp`) advances a 13-variable gas/aerosol/particle state from `1e-4 s` to `2e3 s` on a logarithmically spaced output grid of 301 times. Within each log interval it uses Boost Odeint:

- stepper type: `runge_kutta_fehlberg78<Vector_1D>`;
- controller: `controlled_runge_kutta` with `EPM_ATOLS` and `EPM_RTOLS`;
- integration call: `integrate_adaptive(..., currTimeStep / 100.0, observer)`.

So the EPM state uses an adaptive explicit embedded Runge-Kutta-Fehlberg 7/8 method, initialized with a nominal substep of one hundredth of the current log-time interval.

### EPM equations and couplings

The RHS (`src/EPM/Models/Original/RHS.cpp`) is evaluated directly from parameterized rates:

- entrainment/dilution rate is a piecewise analytic function of time;
- vortex temperature perturbation is a piecewise linear ramp with analytic derivative;
- water condensation/deposition rates use direct formula evaluation;
- homogeneous sulfate nucleation uses AIM nucleation parameterizations;
- sulfate-soot coagulation rates are accumulated by summing over the sulfate PDF bins.

Sulfate aerosol self-coagulation is applied once per EPM log interval before the ODE integration over that interval. Newly nucleated sulfate aerosol is added as a lognormal distribution after each interval when the diagnosed new liquid sulfate exceeds a threshold.

### EPM size discretization

EPM constructs sulfate and ice bins from minimum/maximum radii and a constant volume ratio. Edges are geometric in radius via `radius_ratio = volume_ratio^(1/3)`. Initial and produced distributions use AIM's analytic lognormal PDF sampled at bin centers, where the stored quantity is `dN/dln(r)`.

## 3. Aerosol distribution moments and diagnostics

`AIM::Aerosol` and `AIM::Grid_Aerosol` compute moments with midpoint quadrature in log-radius space:

```text
M_n ≈ Σ_i ln(r_edge[i+1] / r_edge[i]) * r_center[i]^n * pdf[i]
```

Grid diagnostics apply the same bin summations in each cell, then integrate over cell areas where line totals are needed. Radius, effective radius, standard deviation, volume, area, ice water content, extinction, and optical-depth diagnostics are algebraic functions of these moments.

## 4. Coagulation numerics

AIM coagulation (`Aerosol::Coagulate`, `Grid_Aerosol::Coagulate`) uses the volume-conserving discrete scheme cited in comments as Jacobson, *Fundamentals of Atmospheric Modeling* (2005).

For each bin, it builds:

- a production term from smaller/source bins coagulating into the target bin;
- a loss coefficient from target-bin particles coagulating with all other bins.

The update is semi-implicit in the loss term:

```text
v_new[i] = (v_old[i] + dt * P[i]) / (1 + dt * L[i])
```

where `v` is particle-volume concentration. This is more stable and mass-conserving than a fully explicit update. Coagulation destination fractions `f[k][j][i]` are precomputed from bin volumes; physical kernels include Brownian, differential-emission/turbulence-like terms, and for ice additional gravitational/turbulent components depending on phase.

In the active LAGRID plume loop, ice coagulation is present as AIM capability but the observed current loop primarily schedules transport and ice growth.

## 5. Settling numerics

`AIM::SettlingVelocity` computes one terminal fall speed per ice radius bin. The active branch uses Stokes-law fall speed with Cunningham slip correction:

```text
v_fall = 2 g (rho_ice - rho_air) r^2 / (9 mu_air) * Cc(Kn)
```

The vector of fall speeds is computed once from bin centers and reference met state, then passed to transport as bin-dependent downward velocity. An alternative Sölch/Kärcher-style Reynolds-number formula is present but disabled by a local constant.

## 6. Finite-volume transport numerics

Transport is implemented by `FVM_ANDS` and called from `LAGRIDPlumeModel::runTransport()`.

### Equation and grid

Fields are transported on the current structured `x-y` grid. The discretized scalar equation is an advection-diffusion equation with source support:

```text
∂phi/∂t + u ∂phi/∂x + v ∂phi/∂y
  = ∂/∂x(Dh ∂phi/∂x) + ∂/∂y(Dv ∂phi/∂y) + source
```

Velocity is represented by per-cell vectors. Horizontal velocity includes shear:

```text
u(j) = u_base - y(j) * shear
v(j) = v_base
```

For ice bins, `v_base` is set to negative fall speed. H2O and contrail tracer use zero settling velocity.

### Operator splitting

The active plume code uses `FVM_Solver::operatorSplitSolve2DVec`, which performs Strang splitting:

1. explicit advection for a half timestep;
2. implicit diffusion for the full timestep;
3. explicit advection for the second half timestep.

The advection substep is subdivided by a Courant limit. The code computes the current Courant number, chooses a maximum explicit advection timestep from `courant_max`, then takes `ceil(0.5 * dt / dt_adv)` equal substeps for each half step.

### Explicit advection scheme

Advection uses forward Euler in time for each advection substep. Face values are reconstructed with upwind-biased, minmod-limited slopes. This gives a TVD-style second-order correction in smooth regions while falling back toward first-order upwind near extrema/discontinuities.

Boundary faces use boundary values directly. Current implemented boundary condition support is Dirichlet through ghost cells; periodic branches exist but throw if used.

### Implicit diffusion scheme

Diffusion uses a sparse five-point finite-volume stencil assembled in `AdvDiffSystem::buildCoeffMatrix(true)`:

```text
C = 1 + 2 dt (Dh/dx^2 + Dv/dy^2)
E,W = -Dh dt/dx^2
N,S = -Dv dt/dy^2
```

For operator splitting, advection terms are omitted from the implicit matrix. The matrix includes ghost-point equations for boundary conditions. In the plume loop, the diffusion matrix is built once for a representative solver and shared across aerosol-bin solvers when grid and diffusion coefficients match.

Although `FVM_Solver` contains an Eigen BiCGSTAB path, the operator-split diffusion path currently calls the custom sparse SOR routine. `sor_solve` performs Gauss-Seidel/SOR sweeps until relative residual `||A phi - rhs||_2 / ||rhs||_2` falls below the threshold, with a fixed number of inner row sweeps between residual checks.

### Field-specific details

- Ice aerosol PDF bins use enhanced plume-region diffusion. `updateDiffVecs()` sets enhanced coefficients where ice number exceeds `1e-4` of the maximum ice number, otherwise background diffusivity.
- H2O transport solves for the anomaly `H2O - H2O_background` with zero boundary anomaly, then adds the background field back. Numerically this prevents diffusion from eroding the ambient humidity profile.
- Contrail tracer uses the same transport method as H2O but directly on the tracer field.
- Very small fields are skipped if their vector norm is below `1e-100`, avoiding precision problems in Eigen/norm calculations.

## 7. Ice growth and sublimation numerics

`AIM::Grid_Aerosol::Grow` uses an Analytical Predictor of Condensation (APC) scheme, documented in code as Jacobson (1997), for deposition/sublimation coupled to gas-phase water.

For each active grid cell:

1. Convert current ice volume in all bins to equivalent water molecules and add to gas H2O to form total water.
2. Compute per-bin growth coefficients:

```text
k_i = 1e6 * N_i * 4π r_i * D_eff(r_i, T, P, H2O)
```

3. Sum `k_i` and Kelvin-corrected `k_i * Kelvin(r_i)` over bins.
4. Analytically update gas molar water concentration with a backward/predictor expression:

```text
C_qt = (H2O/Na + dt * Σ(k_i Kelvin_i) * C_sat_ice)
       / (1 + dt * Σ k_i)
```

5. Update each bin's ice molar/volume concentration from the new gas concentration.
6. Clip ice volume to nonnegative values and to the maximum allowed bin volume times particle number.
7. Reset gas H2O to total water minus updated ice water.

This scheme couples all bins in a cell through the shared vapor reservoir and is more stable than a purely explicit growth update.

### Moving-center remapping in size space

After growth/sublimation, `ComputeBinParticleFlux` calculates each old bin's new particle volume per particle and finds the destination volume bin with `lower_bound` on bin volume edges. Particles below the minimum volume are removed; particles above the maximum are retained in the largest bin.

`ApplyBinParticleFlux` then conservatively sums particle number and ice volume into each destination bin. The destination bin's local moving volume center is set to `iceVol / icePart`, clipped to that bin's edges; the PDF is recomputed as:

```text
pdf = particle_number / ln(r_edge[i+1] / r_edge[i])
```

## 8. Meteorology interpolation and hydrostatic grid numerics

Meteorology is effectively one-dimensional in the vertical for pressure, altitude, shear, humidity, and vertical velocity, then broadcast in `x` where appropriate.

### Time interpolation

Met time series are stored in memory and interpolated by `interpMetTimeseriesData`. The code chooses the enclosing time index from `simTime_h / MET_DT`, clamps to the last record, and linearly interpolates between records. Non-time-series fields use the first record.

### Vertical interpolation

Spatial interpolation to the simulation grid uses either nearest neighbor or linear interpolation (`met::linInterpMetData`) depending on input flags. Relative humidity over ice is converted to H2O concentration using `physFunc::RHiToH2O`; temperature perturbations are added cellwise after base interpolation.

### Hydrostatic pressure/altitude mapping

`estimateMetDataAltitudes`, `estimateSimulationGridPressures`, and `recalculateSimulationGrid` use hydrostatic formulas with piecewise-constant lapse rate between met levels. For near-zero lapse rate, they switch to the logarithmic isothermal limit.

Pressure velocity/updraft is applied explicitly by adding `omega * dt` to all pressure centers and edges. After this pressure-coordinate displacement and met update, `recalculateSimulationGrid()` inverts pressure edges back to altitudes/y-edges using the same piecewise-hydrostatic formulas.

### Temperature perturbation numerics

`updateTempPerturb()` samples two uniform random values in `[-1, 1]` per cell and multiplies them with the configured amplitude. This is a cellwise random perturbation update, not a spatially correlated turbulence solver.

## 9. Pressure-coordinate unit conversion

Before the vertical met-grid update, the plume loop divides H2O, contrail tracer, and aerosol PDFs by a local number-density proxy:

```text
localND = (p_edge[j] - p_edge[j+1]) / (y_edge[j+1] - y_edge[j])
field /= localND
```

This is a discrete hydrostatic conversion tied to the pressure-edge grid. It makes subsequent grid movement/remapping operate on mixing-ratio-like quantities rather than raw concentrations. After remapping/regeneration, fields continue through the next transport/growth cycle with the model's expected units.

## 10. LAGRID adaptive remapping numerics

LAGRID remapping is a first-order conservative overlap remap.

### Mask and target grid

The plume model builds a binary mask from contrail tracer/ice/H2O criteria and `VectorUtils::MaskInfo` stores the bounding box. `remapVariable()` chooses new spacings from mask extent with hard limits:

```text
dx_new = max(20 m, min(width / 50, 50 m))
dy_new = max(5 m,  min(depth / 50, 7 m))
```

It then creates enough cells to cover the active region plus buffer cells.

### Conservative area-overlap remap

`rectToBoxGrid` converts old active structured cells into `MassBox` records with:

```text
mass = phi * old_cell_area
```

`mapToStructuredGrid` loops over each old mass box and every new grid cell it overlaps. For overlap area `A_overlap`, the contribution to the new concentration is:

```text
phi_new += (mass_old * A_overlap / area_old) / area_new
```

Thus the method is conservative for cell-integrated mass when all overlapping boxes are retained. It is first-order: old-cell values are piecewise constant, and no slope reconstruction is used during remapping.

`getUnusedFraction` computes how much of each new cell was not covered by old boxes. Buffers are added by padding coordinates/arrays with configured fill values.

## 11. Diagnostics/output numerics

Diagnostics are post-processing reductions rather than prognostic solvers. They combine:

- log-bin aerosol moments;
- cell-area integration;
- 1-D/2-D optical-depth accumulations;
- NetCDF writes of scalar and gridded arrays.

Because diagnostics are derived from model state after transport/growth/remap, their accuracy is limited by the transport grid, size-bin quadrature, remapping conservativeness, and unit conversions described above.

## 12. Chemistry numerics

KPP-generated chemistry uses stiff Rosenbrock integrators (`src/KPP/KPP_Integrator.cpp`) for `INTEGRATE`, with generated Jacobian/sparse linear-algebra support. The file includes several L-stable or stiffly stable Rosenbrock method definitions and adaptive error control using absolute/relative tolerances.

In the active LAGRID path, KPP is mainly involved in EPM/background initialization. If full resolved-plume chemistry is re-enabled or expanded, its Rosenbrock tolerances, global/threadprivate state, and heterogeneous-rate update cadence become part of the numerical coupling.

## 13. Main numerical conservation and stability points

- EPM uses adaptive high-order explicit ODE integration but also applies coagulation and nucleated-aerosol insertion at the coarser log-time interval boundaries.
- Coagulation and LAGRID remapping are designed to conserve integrated volume/mass-like quantities.
- Transport is finite-volume and conservative in form, but boundary conditions, masks, remapping retention, and field-specific background subtraction affect global budgets.
- Advection stability is controlled by explicit CFL subcycling; diffusion stability is handled implicitly.
- Ice growth uses an analytical/semi-implicit vapor-coupled update to avoid explicit condensation instability.
- Meteorological vertical motion is operator-split from transport/growth and applied through pressure-edge displacement plus hydrostatic remapping.
