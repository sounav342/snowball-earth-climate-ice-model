# Model Architecture

## Coupling Sequence

`FISEBM` performs repeated asynchronous coupled cycles at the configured forcing:

1. Run `EBM` and `integrate_EBM_1d_sphere`.
2. Write the atmospheric-ocean state to the native EBM restart.
3. Run `fis_1D_sphere` and `integrate_h_1d_sphere` from that state.
4. Write the ice-model state and coupled diagnostics.
5. Use the saved state to initialize the next native coupling cycle.

## State Variables

The restart state carries atmospheric temperature `Ta`, ocean temperature `To`, ice-surface temperature `Ts`, atmospheric humidity `qa`, ice thickness `h`, fractional ice cover `R`, salinity `Sal`, pressure-dependent freezing temperature `T_f`, ocean depth `h_o`, and ice-flow fields. Three time columns retain previous, current, and next levels where required by the native integrators.

## Radiation and Surface Fluxes

Fractional ice cells use the existing `R` field. The effective surface albedo is

```text
alpha_surface = R alpha_i + (1 - R) alpha_o.
```

The implementation retains separate per-area ice and ocean shortwave fluxes and applies the area weighting once in the grid-cell energy budget. Longwave and sensible heat fluxes use the same ice/ocean area weighting. Atmospheric shortwave absorption remains separate.

CO2 forcing enters additively in the longwave energy tendency. Its model sign convention is positive warming and negative cooling.

## Principal Files

- `integrate_EBM_1d_sphere.m`: atmosphere, slab-ocean, humidity, salinity, and surface-temperature integration;
- `integrate_h_1d_sphere.m`: ice thickness, ice flow, subgrid ice fraction, and ice/ocean coupling;
- `FISEBM.m`: asynchronous atmosphere-ocean/ice coupling driver;
- `EBM.m`: atmospheric-ocean component entry point;
- `fis_1D_sphere.m`: ice component entry point;
- `mkrestart.m`: initial restart-state construction.

## Release Cleanup

The equations and parameter values are retained from the working model. Cleanup is limited to removing inactive/commented alternatives, editor backups, duplicate old solvers, dead diagnostics, and interactive debugger stops. Invalid states now raise MATLAB errors, which is safer for unattended batch runs.
