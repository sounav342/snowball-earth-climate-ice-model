# Reproducibility Checklist

Before archiving a release with a paper:

1. Record the Git commit hash and MATLAB release.
2. Preserve the exact structure passed to `run_single_experiment`.
3. Preserve the selected `src/EBMInput/exp_XX.m` initial-condition file.
4. Record any edits to `set_EBM_parameters.m` and `set_FIS_parameters.m`.
5. Archive the generated restart files, figures, and run log.
6. Run `run_all_tests` from a clean clone.
7. Add the paper DOI and repository archive DOI to `CITATION.cff`.
8. Select and add a software license before public release.

Full physical integrations can take prohibitively long for continuous integration. Automated checks therefore cover core numerical helpers and initial-state construction; publication runs should retain their complete native outputs.

