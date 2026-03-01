# `PrimalResult`

## Description

Result object returned by all PRIMAL solvers. This is a Python `dataclass`
that corresponds to the R S3 class `"primal"`.

## Fields

- `type`: `str` -- problem type (`"Dantzig"`, `"SparseSVM"`, `"Compressed sensing"`, `"Quantile Regression"`)
- `data`: `ndarray` of shape `(n, d)` -- input data matrix
- `response`: `ndarray` of shape `(n,)` -- input response vector
- `beta`: `ndarray` of shape `(d, iterN)` -- coefficient estimates; each column corresponds to a regularization parameter along the solution path
- `beta0`: `ndarray` of shape `(iterN,)` or `None` -- intercept estimates (only for SparseSVM)
- `df`: `ndarray` of shape `(iterN,)` -- degrees of freedom (number of nonzero coefficients) along the solution path
- `value`: `ndarray` of shape `(iterN,)` -- optimal objective function values along the path
- `iterN`: `int` -- number of iterations performed
- `lambda_`: `ndarray` of shape `(iterN,)` -- regularization parameters along the solution path
- `runtime`: `float` -- elapsed time in seconds

## Methods

- [`summary()`](summary.md) -- return a formatted summary string
- [`coef(n=None)`](coef.md) -- extract coefficients at a specific iteration
- [`plot(n=None)`](plot.md) -- plot the regularization path

## Notes

- `print(result)` calls `summary()` internally.
- `lambda_` uses a trailing underscore because `lambda` is a Python keyword.
  In R, the corresponding field is `$lambda`.

## See also

[`dantzig_solver`](dantzig_solver.md),
[`sparse_svm_solver`](sparse_svm_solver.md),
[`compressed_sensing_solver`](compressed_sensing_solver.md),
[`quantile_regression_solver`](quantile_regression_solver.md)
