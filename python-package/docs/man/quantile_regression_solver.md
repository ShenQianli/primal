# `quantile_regression_solver`

## Usage

```python
quantile_regression_solver(
    X,
    y,
    max_it=50,
    lambda_threshold=0.01,
    tau=0.5,
) -> PrimalResult
```

## Description

Solve the sparse quantile regression problem using the parametric simplex
method. The method estimates the conditional quantile of the response
variable while encouraging sparsity in the coefficient vector.

## Arguments

- `X`: `ndarray` of shape `(n, d)` -- data matrix.
- `y`: `ndarray` of shape `(n,)` -- response vector.
- `max_it`: `int` -- maximum number of iterations for the solution path. Default `50`.
- `lambda_threshold`: `float` -- the algorithm stops when the regularization parameter falls below this threshold. Default `0.01`.
- `tau`: `float` -- quantile level, must be in `(0, 1)`. Default `0.5` (median regression).

## Returns

`PrimalResult` with:

- `type`: `"Quantile Regression"`
- `data`: the `(n, d)` data matrix from the input
- `response`: the length `n` response vector from the input
- `beta`: a `(d, iterN)` matrix of coefficient estimates whose columns correspond to regularization parameters along the solution path
- `beta0`: `None` (not applicable for quantile regression)
- `df`: degrees of freedom (number of nonzero coefficients) along the solution path
- `value`: the sequence of optimal objective function values corresponding to each `lambda`
- `iterN`: the number of iterations performed
- `lambda_`: the sequence of regularization parameters obtained along the solution path
- `runtime`: elapsed time in seconds

## See also

[`dantzig_solver`](dantzig_solver.md),
[`sparse_svm_solver`](sparse_svm_solver.md),
[`compressed_sensing_solver`](compressed_sensing_solver.md)

## Examples

```python
import numpy as np
from pyprimal import quantile_regression_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
y = X @ np.array([1, -1, 2] + [0] * 17) + rng.standard_normal(n)

result = quantile_regression_solver(X, y, tau=0.5)
print(result)
```
