# Snowball Earth Coupled Climate-Ice Model

This repository contains a zonally averaged Snowball Earth model that couples:

- a one-dimensional atmospheric energy-balance model;
- a saline slab ocean;
- a dynamically evolving floating ice shelf/sea glacier;
- a subgrid fractional ice-cover parameterization.

Each experiment uses a fixed additive CO2 longwave forcing. Positive `FCO2` warms the model and negative `FCO2` cools it; forcing units are W m^-2.

## Requirements

- MATLAB R2025b or later is recommended. The release is tested with R2025b.
- The model requires MATLAB only.
- Expect coupled climate-ice experiments to require substantial wall time and storage.

## Quick Start

Clone the repository, open MATLAB in its root folder, and run:

```matlab
setup_snowball
run_single_experiment;
```

Override selected run values with a scalar structure:

```matlab
run_single_experiment(struct( ...
    'FCO2',35, ...
    'solarConstant',1285, ...
    'atmosphericEmissivity',0.70, ...
    'experimentNumber',1, ...
    'flowMode','flow'));
```

See [Configuration](docs/CONFIGURATION.md) for all run options and initial-condition selection.

## Output

Runs write native restart and diagnostic products beneath `EBMRestart/`, `FISRestart/`, `FISEBMRestart/`, `EBMFigures/`, `FISFigures/`, and `Figures/`.

Generated state, figures, tables, and logs are excluded from Git by `.gitignore`.

## Repository Layout

```text
src/                 Model equations, solvers, and coupling code
src/EBMInput/        Available initial-condition scenario scripts
examples/            User-facing experiment entry points
tests/               Fast MATLAB tests
docs/                Model, configuration, and reproducibility notes
setup_snowball.m     MATLAB path setup
```

The publication release intentionally excludes restart products, figures, autosave files, notebooks used for exploration, and superseded `_old` or numbered solver copies.

## Tests

Run the fast test suite with:

```matlab
setup_snowball
results = run_all_tests;
```

The checks performed while preparing this release are recorded in [Validation](docs/VALIDATION.md).

## Initial-Condition Files

The original working directory contains scenario files `exp_01`, `exp_03`, `exp_04`, and `exp_10` through `exp_22`; those are the scenarios included here. Requests for absent scenario files fail with a clear error instead of silently selecting another state.

## Scientific Lineage

The fractional-cell and coupled sea-glacier formulation follows the model family described by Pollard and Kasting (2005), [doi:10.1029/2004JC002525](https://doi.org/10.1029/2004JC002525), with asynchronous climate/sea-glacier coupling related to Pollard et al. (2017), [doi:10.1002/2017JD026621](https://doi.org/10.1002/2017JD026621).

## Citation and License

Use `CITATION.cff` as the repository citation record and update it with the associated paper DOI when available. No reuse license has been selected in this release folder; choose an appropriate license before making the GitHub repository public.

See [Code availability](docs/CODE_AVAILABILITY.md) for a manuscript-ready statement.
