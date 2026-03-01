"""Example: Dantzig selector for sparse linear regression."""

import numpy as np
from pyprimal import dantzig_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
beta_true = np.array([1, -1, 2, -2, 3] + [0] * 15)
y = X @ beta_true + 0.1 * rng.standard_normal(n)

result = dantzig_solver(X, y)
print(result)
print()
result.coef()
result.plot()
