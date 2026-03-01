# Quick Start

## 0) Runtime check

```python
import pyprimal
pyprimal.test()
```

## 1) Dantzig selector

```python
import numpy as np
from pyprimal import dantzig_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
beta_true = np.array([1, -1, 2, -2, 3] + [0] * 15)
y = X @ beta_true + 0.1 * rng.standard_normal(n)

result = dantzig_solver(X, y)
print(result)
result.coef()
```

## 2) Sparse SVM

```python
import numpy as np
from pyprimal import sparse_svm_solver

rng = np.random.default_rng(42)
n, d = 100, 20
X = rng.standard_normal((n, d))
y = np.sign(X @ np.array([1, -1, 2] + [0] * 17))

result = sparse_svm_solver(X, y)
print(result)
result.coef()
```

## 3) Compressed sensing

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

## 4) Quantile regression

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

## 5) Visualization

```python
result.plot()       # all three plots
result.plot(1)      # coefficients vs lambda
result.plot(2)      # coefficients vs iteration
result.plot(3)      # lambda vs iteration
```

## 6) Coefficient extraction

```python
result.coef()       # last iteration (default)
result.coef(1)      # first iteration (1-based index)
result.coef(5)      # fifth iteration
```
