# Contributing

Keep changes narrowly scoped and preserve the physical equations unless a scientific change is explicitly documented. New behavior should include a focused MATLAB test and a note describing whether existing equilibria or restart compatibility may change.

Before opening a pull request, run:

```matlab
setup_snowball
run_all_tests
```

Do not commit generated restart files, figures, tables, or run logs.
