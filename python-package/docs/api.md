# API Reference

## Solvers

### `pyprimal.dantzig_solver`

```python
dantzig_solver(X, y, max_it=50, lambda_threshold=0.01) -> PrimalResult
```

Solve the Dantzig selector problem via parametric simplex method.

### `pyprimal.sparse_svm_solver`

```python
sparse_svm_solver(X, y, max_it=50, lambda_threshold=0.01) -> PrimalResult
```

Solve the sparse SVM classification problem via parametric simplex method.

### `pyprimal.compressed_sensing_solver`

```python
compressed_sensing_solver(X, y, max_it=50, lambda_threshold=0.01) -> PrimalResult
```

Solve the compressed sensing problem via parametric simplex method.

### `pyprimal.quantile_regression_solver`

```python
quantile_regression_solver(X, y, max_it=50, lambda_threshold=0.01, tau=0.5) -> PrimalResult
```

Solve the sparse quantile regression problem via parametric simplex method.

## Result object

### `pyprimal.PrimalResult`

Dataclass returned by all solvers with fields:

- `type`, `data`, `response`, `beta`, `beta0`, `df`, `value`, `iterN`, `lambda_`, `runtime`

### Methods

```python
result.summary() -> str
result.coef(n=None) -> dict
result.plot(n=None)
```

- `summary()`: formatted multi-line summary (also used by `print(result)`)
- `coef(n)`: extract coefficients at iteration `n` (1-based index)
- `plot(n)`: visualize the regularization path

## Runtime helper

### `pyprimal.test`

```python
pyprimal.test()
```

Verify package installation and C library availability.
