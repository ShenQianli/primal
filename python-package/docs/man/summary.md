# `PrimalResult.summary`

## Usage

```python
result.summary() -> str
```

## Description

Return a formatted multi-line summary string including the model type,
regularization parameter path, and degrees of freedom along the solution
path. This is the Python equivalent of R's `print.primal()`.

Calling `print(result)` or `str(result)` invokes `summary()` internally.

## Returns

`str` with columns:

- `lambda`: regularization parameters
- `df`: degrees of freedom
- `value`: objective function values

## See also

[`PrimalResult`](PrimalResult.md),
[`PrimalResult.coef`](coef.md)
