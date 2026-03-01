"""Example: Compressed sensing for sparse signal recovery."""

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
print()
result.coef()
result.plot()
