# `PrimalResult.plot`

## Usage

```python
result.plot(n=None)
```

## Description

Plot the regularization path and parameters obtained from the parametric
simplex method. This is the Python equivalent of R's `plot.primal()`.

Requires `matplotlib`.

## Arguments

- `n`: `int` or `None` -- if `None`, all three plots are shown together in a single figure. If `n` is a number (1, 2, or 3), only the corresponding plot is shown:
    - `1`: coefficient estimates vs. regularization parameter (lambda)
    - `2`: coefficient estimates vs. iteration index
    - `3`: regularization parameter (lambda) vs. iteration index

## Notes

- Requires `matplotlib >= 3.5`. Install with `pip install "pyprimal[viz]"`.

## See also

[`PrimalResult`](PrimalResult.md),
[`PrimalResult.coef`](coef.md),
[`dantzig_solver`](dantzig_solver.md),
[`sparse_svm_solver`](sparse_svm_solver.md)
