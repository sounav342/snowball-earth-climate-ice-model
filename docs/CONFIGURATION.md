# Configuration

Use `run_single_experiment` to run one asynchronously coupled atmosphere-ocean/ice experiment at fixed CO2 forcing.

## Run Options

| Field | Default | Meaning |
|---|---:|---|
| `FCO2` | 35 | Additive longwave forcing in W m^-2 |
| `solarConstant` | 1285 | Solar constant in W m^-2 |
| `atmosphericEmissivity` | 0.70 | Atmospheric longwave emissivity |
| `experimentNumber` | 1 | Initial-condition script number in `src/EBMInput` |
| `flowMode` | `flow` | Use `flow` or `no-flow` for ice dynamics |
| `startCycle` | 1 | First asynchronous coupling-cycle index |

Positive `FCO2` produces warming and negative `FCO2` produces cooling in the model equations. The value remains fixed throughout a run.

## Example

```matlab
setup_snowball
run_single_experiment(struct( ...
    'FCO2',20, ...
    'solarConstant',1292, ...
    'atmosphericEmissivity',0.70, ...
    'experimentNumber',1, ...
    'flowMode','flow', ...
    'startCycle',1));
```

Unknown option names raise an error to prevent misspelled settings from being ignored.

## Initial Conditions

Initial conditions are defined by the numbered scripts in `src/EBMInput`. Select one with `experimentNumber`; for example, `experimentNumber=1` loads `exp_01.m`. Each script defines the ice edges, temperature-profile choice, critical ice thickness, and other scenario-specific values used by `mkrestart`.

To add a scenario, copy an existing input script to a new `exp_XX.m` file, edit its values, and use the matching experiment number. Do not overwrite a published scenario when reproducibility matters.

## Physical Parameters

Atmosphere-ocean parameters are defined in `set_EBM_parameters.m`, and ice-model parameters are defined in `set_FIS_parameters.m`. Record any edits to these files with the experiment metadata.

