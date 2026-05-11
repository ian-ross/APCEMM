# YamlInputReader subsystem

`YamlInputReader` converts user YAML plus embedded defaults into the `OptInput` configuration and per-case parameter sets consumed by `Core::Input`.

## Main code

- `include/YamlInputReader/YamlInputReader.hpp`
- `src/YamlInputReader/YamlInputReader.cpp`
- defaults generated from `Code.v05-00/defaults/input.yaml` into `include/Defaults/Input.hpp`

## Read path

The executable calls:

```text
YamlInputReader::readYamlInputFiles(optInput, filenames)
YamlInputReader::generateCases(optInput)
```

`readYamlInputFiles` loads defaults and user-supplied YAML files, then dispatches to menu readers:

- simulation menu -> output, parallelism, seeding, model/EPM type, save flags;
- parameter menu -> scalar values, sweeps, Monte Carlo settings;
- transport menu -> transport/updraft timestep settings;
- chemistry menu -> chemistry and J-rate options;
- aerosol menu -> settling, coagulation, ice-growth options;
- meteorology menu -> met-file/time-series loading options;
- diagnostic/timeseries menus -> output frequencies/variables;
- advanced menu -> grid, EPM, and numerical tuning options.

After parsing, `performOtherInputValidnessChecks` enforces cross-field constraints that are not local YAML type checks.

## Case generation

`OptInput::PARAMETER_PARAM_MAP` stores each parameter name as a vector of numeric values. `generateCases` expands these into `std::vector<std::unordered_map<std::string,double>>`, one map per generated case. `Input` then pulls the case-indexed values from this structure.

Parameter strings support normal scalar parsing plus sweep/Monte-Carlo modes through helpers such as `parseParamSweepInput`, `parseDoubleString`, `parseBoolString`, and vector parsing functions.

## Implementation notes

- Path handling goes through `parseFileSystemPath`; the reader tracks an `INPUT_FILE_PATH` so relative paths can be resolved relative to the input file.
- Numeric string validation rejects malformed scientific notation before calling `std::stod`.
- This subsystem should remain independent of physics implementation details: it should validate and normalize options, not perform model calculations.
