# Validation Record

Release 1.0.0 was checked with MATLAB R2025b on 2026-09-05.

## Automated Tests

- Fast MATLAB tests: 3 passed, 0 failed.
- MATLAB parse errors across source, examples, and tests: 0.
- Interactive `keyboard` calls in publication source: 0.
- Required product: MATLAB only.
- YAML parsing for `CITATION.cff` and GitHub Actions: passed.

## Scientific Equivalence

The compact `mkrestart` implementation was compared against the working model for every available input scenario: `exp_01`, `exp_03`, `exp_04`, and `exp_10` through `exp_22`. All eight saved state fields (`h`, `qa`, `R`, `Ta`, `To`, `Ts`, `T_f`, and `Sal`) were identical for all 16 scenarios, covering 128 array comparisons.

The physical equations, grid, native asynchronous coupling sequence, and default parameter values were retained during release cleanup. A full-duration climate-ice integration was not repeated as a packaging test because of its substantial runtime.

