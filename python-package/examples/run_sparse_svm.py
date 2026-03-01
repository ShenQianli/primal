"""Example: Sparse SVM for classification."""

import numpy as np
from pyprimal import sparse_svm_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
y = np.sign(X @ np.array([1, -1, 2] + [0] * 17))

result = sparse_svm_solver(X, y)
print(result)
print()
result.coef()
result.plot()
