# pyprimal 2.0 Documentation

`pyprimal` is a native Python package for sparse learning via the parametric
simplex method, mirroring the R `PRIMAL` package.

If you are new, follow this order:

1. [Installation](installation.md)
2. [Quick Start](quickstart.md)
3. [Function Manual](man/index.md)

## Package capabilities

- Solvers: `dantzig_solver`, `sparse_svm_solver`, `compressed_sensing_solver`,
  `quantile_regression_solver`
- Result methods: `PrimalResult.summary`, `PrimalResult.coef`, `PrimalResult.plot`
- Environment probe: `pyprimal.test()`

## Runtime model

`pyprimal` uses a compiled C++ core (`libpsm`) wrapped via `ctypes`:

- Core algorithm in C++ with Eigen for linear algebra
- Python layer handles input validation, result packaging, and visualization
- No R runtime dependency

## Documentation map

- Practical pages: [Installation](installation.md), [Quick Start](quickstart.md)
- Reference pages: [API Reference](api.md), [Function Manual](man/index.md)
- [Citation](citation.md)
