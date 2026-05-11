# Potential general code issues

This document records potential non-numerical code risks found while reviewing the active APCEMM/LAGRID path described in [Architecture](ARCHITECTURE.md), [Physical calculation flow](physics.md), and [Numerical methods](numerics.md). Numerical convergence/accuracy risks are covered separately in [Potential numerical convergence and accuracy issues](numerics-issues.md).

These items should be treated as audit findings: some may be harmless for current default inputs, but they are worth testing or hardening.

## 1. Input validation and bad parameter values

- **Several physical/grid parameters are parsed but not fully range-checked.** Examples include zero grid cell counts/spacings, reversed or zero-width domains, zero/negative diffusivities, invalid humidity/temperature/pressure, and size-bin limits/ratios. Some checks exist in `YamlInputReader`, but many parsed doubles are accepted as-is.
- **Grid limit validation only checks non-negativity.** `ADV_GRID_XLIM_LEFT`, `ADV_GRID_XLIM_RIGHT`, `ADV_GRID_YLIM_UP`, and `ADV_GRID_YLIM_DOWN` are required to be nonnegative, but there is no obvious check that the resulting domain has positive width/depth or enough cells for all stencil operations.
- **Parameter sweep step values are not guarded.** `parseParamSweepInput()` accepts `min:step:max` without checking that `step` is positive and directionally consistent with `min`/`max`. A zero step can produce an infinite loop; a negative step with `min < max` can also fail to terminate.
- **Parameter sweep empty-vector edge case.** For `min:step:max`, if the loop inserts no values, the code reads `paramVector[paramVector.size()-1]`, which is out of bounds.
- **Monte Carlo parsing appears inconsistent with its error message.** The code says Monte Carlo accepts `min:max` or a singular constant, but `parseParamSweepInput()` rejects a one-token value because it checks `colon_split_tokens.size() != 0 && != 2`. This should be verified against tests and documentation in [YAML input reader](YAML_INPUT_READER.md).
- **Monte Carlo min/max ordering is unchecked.** `fRand(min, max)` will return values on the interval implied by the arithmetic even if `min > max`; this may silently invert user intent.
- **Initialization distributions can divide by zero.** `LAGRID::initVarToGrid()` scales by `mass / newMass`. If a rectangle/Gaussian has zero overlap with the grid, zero width/depth, or a zero weight integral, `newMass` can be zero.
- **Environment assumptions in output metadata.** `CreateREADME()` streams `SLURM_JOB_USER` and `SLURM_NODELIST` directly. On non-SLURM systems these environment variables may be null, which should be handled explicitly.

## 2. Off-by-one, bounds, and indexing risks

- **Stencil code requires adequate ghost/interior dimensions.** Transport limiters and neighbor lookups in `FVM_ANDS` assume neighboring points exist. Very small `Nx`/`Ny` values may exercise boundary paths not covered by normal cases.
- **Potential limiter logic typo.** In `AdvDiffSystem.hpp`, one south-face limiter branch tests `neighbor_point(FaceDirection::NORTH, pointID)` directly in a boolean expression. A valid nonzero point index will evaluate true and can force the limiter to zero. See the related numerical note in [numerics-issues.md](numerics-issues.md#6-finite-volume-transport).
- **Mask bounding boxes can be wrong for first selected cells.** `VectorUtils::mask(..., xEdges, yEdges, ...)` uses `else if` when updating `maxX`/`maxY` after `minX`/`minY`. For the first masked cell, updating the min prevents updating the max in the same coordinate. This can understate the active region when only a few cells are masked.
- **`numeric_limits<double>::min()` is used where `lowest()` is intended.** `VectorUtils::max()` and `FreeCoordBoxGrid::setMinMaxCoords()` initialize maxima with the smallest positive normalized double, not the most negative double. If coordinates or field values are negative, maxima may be wrong. `VectorUtils::mask()` correctly uses `lowest()` and can be used as a model.
- **Axis convention is easy to misuse.** `Vector_2D` is stored as `[y][x]`, while some helper names and arguments refer to `n_x`, `n_y`. This is documented in comments but remains a common source of transposed indexing bugs.
- **Remapping index conversion truncates through double-to-int.** `std::floor(...)` results are clamped with `std::max/std::min` against doubles and assigned to `int`. This is probably safe for normal grids, but boundary cells exactly on remap edges should be regression-tested.
- **Status reporting can use an uninitialized status.** In `Main.cpp`, `SimStatus case_status;` is declared before a switch. For unimplemented model values or unexpected switch paths, `case_status` may be read before assignment.

## 3. Multithreading and shared-state issues

- **`Input_Opt` is shared and mutated inside the parallel case loop.** `Main.cpp` uses `#pragma omp parallel for ... shared(Input_Opt, ...)`, then assigns `Input_Opt.TS_AERO_FILENAME` per case. This is a data race and can make per-case aerosol time-series filenames nondeterministic.
- **Global C RNG is not thread-safe or stream-safe.** `setSeed()` calls `srand()` once and `fRand()` uses `rand()`. Case generation currently happens before the parallel case loop, but any future random draws inside parallel physics would race and/or produce schedule-dependent sequences.
- **OpenMP nested parallelism is intentionally conditional but fragile.** Several kernels use `if (!PARALLEL_CASES)` to avoid nesting when cases are parallelized. This relies on the global `PARALLEL_CASES` state being correct and consistently visible.
- **KPP global/threadprivate state remains a risk.** [KPP](KPP.md) notes mutable generated globals, some declared `threadprivate`. Every thread that uses KPP must initialize its own state; adding resolved-plume chemistry or new parallel calls can easily introduce cross-thread contamination.
- **Global diagnostic flags may affect all cases.** `Diag::set_storePSD(true)` appears process-global. If future inputs allow per-case changes, shared diagnostic state will need isolation.
- **Console output is only partly protected.** Many messages are wrapped in `omp critical`, but not all prints are. Interleaved logs are mostly cosmetic, but they can obscure failures in parallel runs.
- **Exceptions in OpenMP regions need careful handling.** Several utility functions throw on invalid values. Uncaught exceptions escaping OpenMP parallel regions can terminate the process rather than marking only one case as failed.

## 4. File and output handling

- **Existence check plus write is not atomic.** The parallel case loop checks whether an output file exists and later writes it. Unique case filenames reduce risk, but restarts or duplicate case numbers can still race.
- **Output directory creation is old-style and non-recursive.** `mkdir()` does not create parent directories and error handling does not distinguish existing files, permission failures, or invalid paths. The code itself has a TODO to move to `std::filesystem`.
- **Status files have no extension and are overwritten.** This may be intentional, but it makes post-processing conventions less explicit and can hide repeated attempts.
- **NetCDF/thread safety should be confirmed.** Parallel cases may write separate files. That is usually safe, but it depends on the NetCDF library build and any shared APCEMM wrapper state.

## 5. Error handling and failure modes

- **Many parser catch blocks collapse different failures into generic messages.** This helps users but can hide which YAML key was missing or malformed unless parameter names are consistently propagated.
- **Some fatal errors call `exit(1)`.** Direct process exit in `Main.cpp` prevents cleanup and makes the code harder to embed in tests or larger workflows.
- **Unsupported code paths print “Not implemented yet” and continue.** Unimplemented model branches should either throw or set a clear failure status to avoid uninitialized values and misleading completion messages.
- **NaN/Inf checks are uneven.** `VectorUtils::{min,max}` throw on NaN/Inf, but many physics and remapping paths do not check intermediate values before using them in array indexing or file output.
- **Linear solver hangs are possible.** As noted in [numerics-issues.md](numerics-issues.md#7-implicit-diffusion-and-sor-convergence), SOR has no maximum iteration cap; this is both a numerical and operational code robustness issue.

## 6. Memory and resource management

- **Manual allocation appears in utility/output conversion code.** Functions such as `vect2float` return raw pointers according to `Core/Util.hpp`; ownership and deallocation should be audited.
- **Large 3-D aerosol arrays are copied in microphysics.** `Grid_Aerosol::Coagulate()` builds `v = Volume()` and `v_new = v`. This is clear but memory-heavy for large grids/bins and can amplify parallel memory pressure.
- **Stack arrays with runtime size are non-standard C++.** `double P[nBin]; double L[nBin];` in `Grid_Aerosol::Coagulate()` are variable-length arrays accepted by some compilers as an extension. Prefer `std::vector<double>` for portability.
- **Recursive case generation can explode memory.** `generateCasesHelper()` builds the full Cartesian product of sweep values in memory. Large sweeps can exhaust memory before any case starts.

## 7. Data model and unit-consistency risks

- **Mixed units are pervasive.** The code intentionally mixes SI, cgs aerosol conventions, per-meter line totals, and mixing-ratio-like remap variables. This is described in [Numerical methods](numerics.md#9-pressure-coordinate-unit-conversion), but it increases the chance of future unit bugs.
- **Pressure/altitude/grid arrays must remain synchronized.** LAGRID remapping and pressure-coordinate updates change multiple arrays that are assumed to have matching lengths and ordering. Add invariant checks after grid regeneration.
- **Boundary condition support is partial.** Periodic and other branches throw. Input options should prevent unsupported combinations before the solver is constructed.
- **Case-level mutable input objects blur configuration and state.** The same `OptInput` object carries parsed input, generated settings, and mutable per-run filenames. Separating immutable configuration from per-case runtime state would reduce accidental sharing.

## 8. Suggested hardening checks

- Add parser validation for positive grid sizes, positive/finite physical constants, monotone domains, positive time steps, positive bin ratios, and valid parameter-sweep step direction.
- Add unit tests for `parseParamSweepInput()`, including single-value Monte Carlo, `min == max`, zero/negative step, and empty tokens.
- Replace `numeric_limits<double>::min()` maxima initializers with `lowest()` where negative values are possible.
- Fix or test `VectorUtils::mask()` bounding-box updates so min and max are updated independently.
- Make `Input_Opt` immutable inside the parallel case loop; move per-case filenames into `Input` or a local copy.
- Replace `rand()/srand()` with per-thread/per-case C++ RNG objects and recorded seeds.
- Add assertions/invariant checks for vector dimensions after initialization, transport, remap, and met-grid updates.
- Add maximum-iteration/error-return behavior to SOR and wrap per-case exceptions so one failed case does not necessarily abort a whole sweep.
- Run ThreadSanitizer/AddressSanitizer/UndefinedBehaviorSanitizer builds on small representative cases.
