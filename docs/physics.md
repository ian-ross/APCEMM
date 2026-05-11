# Physical calculation flow

This document traces the physical calculations in APCEMM's active plume-model path, as represented by `Code.v05-00/src/Main.cpp` and `Core/LAGRIDPlumeModel`. It focuses on the currently executed LAGRID plume simulation rather than dormant box-model or adjoint-only branches.

For discretization, time integration, solver, remapping, and interpolation details for these processes, see [Numerical methods](numerics.md).

## 1. Inputs become physical case state

YAML defaults and user files are parsed into `OptInput`; parameter sweeps are expanded into per-case maps consumed by `Input`. Each `Input` case carries the physical scenario:

- ambient pressure, temperature, relative humidity, horizontal/vertical diffusivity;
- latitude, longitude, day/time for solar and meteorological context;
- aircraft parameters: speed, mass, wingspan, engine count, fuel flow;
- emission indices for NOx, CO, hydrocarbons, SO2, sulfate conversion, soot mass, and soot radius;
- background species concentrations.

`LAGRIDPlumeModel` then constructs aircraft, fuel, emission, timestep, and meteorology objects. `Meteorology` either builds default ambient profiles or reads NetCDF met data/time series, then interpolates temperature, humidity, shear, vertical velocity, pressure, altitude, and air number density to the simulation grid.

Shared thermodynamic relationships come from `Util/PhysFunction`: saturation vapor pressures over water/ice, humidity conversions, air density/viscosity, gas diffusion, particle fall speeds, and related nondimensional numbers.

## 2. Aircraft and emissions preprocessing

Before EPM is run, the model computes emitted particle scales from aircraft and emission data. In `LAGRIDPlumeModel::runEPM()`:

1. Soot particle volume is computed from the soot radius.
2. Soot particle mass uses the soot material density.
3. The soot mass emission index is converted to an emitted soot number emission index, `#/kg_fuel`.
4. Aircraft fuel burn per distance converts this to emitted soot number per flight distance, `#/m`.

The aircraft/vortex objects also supply geometric plume scales used later, including mature vortex/plume rectangle width, depth, and vertical center displacement.

## 3. Early Plume Model near-field evolution

The Early Plume Model represents the unresolved jet and early vortex/dispersion regime. The factory `EPM::make_epm(...)` selects either the original built-in model or an external NetCDF EPM result.

The original EPM uses aircraft, emissions, ambient meteorology, and model switches to estimate near-field plume evolution. Its physics includes:

- exhaust dilution and entrainment of ambient air;
- jet/vortex temperature evolution;
- sulfur gas/liquid aerosol partitioning;
- soot, sulfate aerosol, and ice particle initialization;
- condensation, deposition, and freezing checks;
- chemistry/background spin-up where required through KPP.

EPM returns `EPM::Output`, which is the physical interface to the resolved plume model:

- plume temperature and water/sulfur gas concentrations;
- sulfate and ice aerosol size distributions;
- ice radius/density and soot density;
- plume area and engine geometry quantities.

Because EPM is effectively run for one engine, APCEMM multiplies `EPM::Output::area` by aircraft engine count before converting EPM concentrations to line totals.

## 4. Post-jet and vortex survival parameterization

After EPM, APCEMM computes the post-jet ice crystal count:

```text
N_postjet = IceAer.Moment(0) * epm_area * 1e6    [#/m]
```

The factor converts from the aerosol concentration convention to per-meter line count. An advanced option can override this post-jet ice number by scaling the EPM ice PDF.

The aircraft vortex parameterization then estimates a survival fraction for vortex sinking losses:

```text
survival_fraction = Aircraft::VortexLosses(N_postjet, exhaust_water, reference_N, reference_wingspan)
```

If the fraction is zero or negative, the case ends with `NoSurvivalVortex`. Otherwise, the ice PDF is scaled by the survival fraction. This gives the mature-plume ice population used to initialize the 2-D simulation.

## 5. Initial 2-D plume grid and fields

`initializeGrid()` creates an initial structured `x-y` cross-section from advanced grid options. The grid is centered around the aircraft plume coordinate system.

The mature contrail is initialized as a rectangular plume using aircraft vortex dimensions:

- width: `aircraft_.vortex().width_rect_mature()`;
- depth: `aircraft_.vortex().depth_mature()`;
- vertical center: `-aircraft_.vortex().z_center()`.

For each ice size bin, EPM bin number is converted to a mass/number line quantity and placed into the rectangle with `LAGRID::initVarToGridRectangular`. The result is an `AIM::Grid_Aerosol` PDF on the 2-D grid. Initial total particle number and ice mass are derived by integrating over cell areas.

`initH2O()` initializes water vapor:

1. Start from the meteorological background `H2O` field.
2. Identify cells containing initialized ice.
3. Compute remaining emitted water vapor as exhaust water minus initialized ice mass.
4. Convert water mass to molecules.
5. Spread emitted water evenly over ice-containing cells, normalized by cell area.
6. Initialize a contrail tracer to 1 in ice-containing cells and 0 elsewhere.

If gravitational settling is enabled, AIM computes bin fall speeds from ice bin radius, reference temperature, and pressure.

## 6. Main resolved plume timestep

The resolved plume loop advances until final time or early extinction. The recurring physical sequence is:

```text
transport
  -> temperature perturbation update
  -> ice growth/sublimation
  -> contrail tracer update
  -> unit conversion for pressure-coordinate/met update
  -> vertical advection/met update
  -> conservative remapping/adaptive grid update
  -> diagnostics and early-stop checks
```

### 6.1 Transport of ice, vapor, and tracer

`runTransport()` solves advection-diffusion equations with `FVM_ANDS`. Transport is applied to:

- each ice aerosol number/PDF bin;
- water vapor perturbation relative to background humidity;
- contrail tracer.

The generic scalar equation is approximately:

```text
∂phi/∂t + u ∂phi/∂x + v ∂phi/∂y
    = ∂/∂x(Dh ∂phi/∂x) + ∂/∂y(Dv ∂phi/∂y) + source
```

Important details:

- Representative shear comes from the meteorology at the plume's strongest optical-depth column.
- Ice-bin vertical advection includes gravitational fall speed when settling is enabled.
- Ice uses enhanced diffusion inside high-ice-number regions and background diffusivity elsewhere.
- H2O and contrail tracer use background diffusivity and no settling.
- H2O diffusion is applied to `H2O - H2O_background`, then the background is added back, so ambient vertical humidity gradients are not artificially smoothed away.
- FVM_ANDS uses sparse finite-volume matrices, operator splitting, and limited advection to preserve stability/monotonicity.

### 6.2 Temperature perturbations and meteorology

If turbulent/diurnal temperature perturbations are enabled, `Meteorology::updateTempPerturb()` updates the temperature perturbation field. The code comments note that, with LAGRID remapping every transport timestep, temperature perturbations physically make sense at the transport cadence.

The model computes mid-timestep solar and simulation times, then later updates meteorology with:

```text
met_.Update(transport_dt, solar_time, simulation_time)
```

This can interpolate met time series, update diurnal temperature effects, and update the ambient fields used by subsequent physics.

### 6.3 Ice growth and sublimation

When scheduled, `AIM::Grid_Aerosol::Grow(dt, H2O, T, P)` evolves ice bins and water vapor together. For each grid cell and particle bin, AIM evaluates depositional growth or sublimation using local:

- water vapor concentration;
- temperature;
- pressure;
- particle radius/volume center;
- effective vapor diffusivity and heat-transfer corrections.

Growth removes vapor and increases ice volume; sublimation returns vapor and shrinks ice. Particles can move between size bins through bin flux calculations. Derived quantities such as total ice mass, ice water content, extinction, and effective radius come from moments of the updated distribution.

### 6.4 Contrail influence tracer

After ice growth, the contrail tracer is reinforced wherever ice remains. The model adds the local ice number converted from `#/cm^3` to `#/m^3` and clamps the tracer into `[0, 1]`:

```text
Contrail = clamp(Contrail + TotalNumber * 1e6, 0, 1)
```

A high tracer threshold defines the mask of cells retained for adaptive remapping. If no cells remain, the simulation completes early.

### 6.5 Pressure-coordinate/met coupling

Before vertical advection and met-grid recalculation, fields are converted from concentration-like units to mixing-ratio-like units using the local air number density implied by hydrostatic pressure-edge differences:

```text
localND = (pressure_edge[j] - pressure_edge[j+1]) / (y_edge[j+1] - y_edge[j])
field /= localND
```

The code comment identifies the hydrostatic assumption:

```text
dp/dz = -rho g = -(n/V) M g
```

This conversion is applied to contrail tracer, H2O, and aerosol PDFs before the meteorological vertical-coordinate update. If pressure vertical velocity/updraft is active, `met_.applyUpdraft()` modifies pressure/altitude alignment. `met_.recalculateSimulationGrid()` then determines updated altitude/y edges after vertical advection and temperature changes.

### 6.6 Adaptive conservative remapping

LAGRID remapping lets the computational domain follow the contrail. The mask generated from the tracer/ice/H2O state identifies active cells and bounds. `remapAllVars()` then conservatively maps fields from the old structured grid to a new structured grid:

1. Convert retained cells to `MassBox` objects where `mass = phi * dx * dy`.
2. Select new grid spacing and dimensions from active-mask extent, with min/max resolution limits and buffers.
3. Compute overlap area between old boxes and new cells.
4. Deposit mass into the new grid and recover concentration-like fields.
5. Add buffers around the active region.
6. Regenerate or remap meteorological and diffusion fields on the new coordinates.

This preserves integrated quantities as much as possible while reducing domain size and following plume deformation/growth.

## 7. Diagnostics and early termination

At output times, diagnostics derive physical plume properties from the current aerosol, vapor, tracer, and meteorology fields. AIM computes integrated particle number, ice mass, ice water content, optical depth/extinction, radius, effective radius, and size distributions. `Diag_Mod` writes configured NetCDF outputs.

The model stops early when either:

- no retained contrail-containing cells remain; or
- total remaining particles fall below a small fraction of the initial mature-plume particle number.

Otherwise, the timestep advances and the cycle repeats.

## 8. Chemistry role in this flow

KPP provides generated gas-phase and heterogeneous chemistry integrators, photolysis, and rate updates. In the currently active LAGRID plume path, KPP is primarily part of EPM/background initialization rather than a dominant per-timestep plume process. Its global/threadprivate state remains important if chemistry is expanded in the resolved plume loop.

## 9. Unit conventions and conservation points

Several historical unit conventions are important when following the physics:

- particle radii are in meters;
- aerosol PDFs are treated as concentration distributions, commonly `#/cm^3` per `ln r`;
- water vapor is commonly `molec/cm^3`;
- plume initialization and diagnostics often integrate over cell area to obtain line quantities, e.g. `#/m` or `kg/m`;
- remapping conserves `phi * area` for transported scalar fields;
- humidity and saturation calculations use SI thermodynamic quantities internally and convert at API boundaries.

The main conservation-sensitive transitions are EPM concentration to line total, grid initialization over cell areas, vapor/ice exchange in AIM growth, transport with finite-volume discretization, and LAGRID remapping by overlap area.
