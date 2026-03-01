"""Example: Sparse quantile regression."""

import numpy as np
from pyprimal import quantile_regression_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
y = X @ np.array([1, -1, 2] + [0] * 17) + rng.standard_normal(n)

# Median regression (tau=0.5)
result = quantile_regression_solver(X, y, tau=0.5)
print(result)
print()
result.coef()
result.plot()
