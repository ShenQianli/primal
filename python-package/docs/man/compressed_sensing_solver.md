# `compressed_sensing_solver`

## Usage

```python
compressed_sensing_solver(
    X,
    y,
    max_it=50,
    lambda_threshold=0.01,
) -> PrimalResult
```

## Description

Solve the compressed sensing problem using the parametric simplex method.
Given an underdetermined linear system `y = X @ beta`, the method recovers
the sparsest solution by tracing the L1-regularization path.

## Arguments

- `X`: `ndarray` of shape `(n, d)` -- measurement matrix.
- `y`: `ndarray` of shape `(n,)` -- observation vector.
- `max_it`: `int` -- maximum number of iterations for the solution path. Default `50`.
- `lambda_threshold`: `float` -- the algorithm stops when the regularization parameter falls below this threshold. Default `0.01`.

## Returns

`PrimalResult` with:

- `type`: `"Compressed sensing"`
- `data`: the `(n, d)` data matrix from the input
- `response`: the length `n` response vector from the input
- `beta`: a `(d, iterN)` matrix of coefficient estimates whose columns correspond to regularization parameters along the solution path
- `beta0`: `None` (not applicable for compressed sensing)
- `df`: degrees of freedom (number of nonzero coefficients) along the solution path
- `value`: the sequence of optimal objective function values corresponding to each `lambda`
- `iterN`: the number of iterations performed
- `lambda_`: the sequence of regularization parameters obtained along the solution path
- `runtime`: elapsed time in seconds

## See also

[`dantzig_solver`](dantzig_solver.md),
[`sparse_svm_solver`](sparse_svm_solver.md),
[`quantile_regression_solver`](quantile_regression_solver.md)

## Examples

```python
import numpy as np
from pyprimal import compressed_sensing_solver

rng = np.random.default_rng(42)
n, d = 50, 100
X = rng.standard_normal((n, d))
beta_true = np.zeros(d)
beta_true[:5] = [3, -2, 1, -1, 0.5]
y = X @ beta_true

result = compressed_sensing_solver(X, y)
print(result)
```
