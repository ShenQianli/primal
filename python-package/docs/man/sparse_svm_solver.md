# `sparse_svm_solver`

## Usage

```python
sparse_svm_solver(
    X,
    y,
    max_it=50,
    lambda_threshold=0.01,
) -> PrimalResult
```

## Description

Solve the sparse support vector machine (SVM) classification problem using
the parametric simplex method. The method finds a sparse linear classifier
that separates two classes while controlling the L1 norm of the coefficient
vector.

## Arguments

- `X`: `ndarray` of shape `(n, d)` -- data matrix.
- `y`: `ndarray` of shape `(n,)` -- response vector with values in {-1, 1}.
- `max_it`: `int` -- maximum number of iterations for the solution path. Default `50`.
- `lambda_threshold`: `float` -- the algorithm stops when the regularization parameter falls below this threshold. Default `0.01`.

## Returns

`PrimalResult` with:

- `type`: `"SparseSVM"`
- `data`: the `(n, d)` data matrix from the input
- `response`: the length `n` response vector from the input
- `beta`: a `(d, iterN)` matrix of coefficient estimates whose columns correspond to regularization parameters along the solution path
- `beta0`: a length `iterN` vector of intercept estimates corresponding to regularization parameters along the solution path
- `df`: degrees of freedom (number of nonzero coefficients) along the solution path
- `value`: the sequence of optimal objective function values corresponding to each `lambda`
- `iterN`: the number of iterations performed
- `lambda_`: the sequence of regularization parameters obtained along the solution path
- `runtime`: elapsed time in seconds

## See also

[`dantzig_solver`](dantzig_solver.md),
[`compressed_sensing_solver`](compressed_sensing_solver.md),
[`quantile_regression_solver`](quantile_regression_solver.md)

## Examples

```python
import numpy as np
from pyprimal import sparse_svm_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
y = np.sign(X @ np.array([1, -1, 2] + [0] * 17))

result = sparse_svm_solver(X, y)
print(result)
```
