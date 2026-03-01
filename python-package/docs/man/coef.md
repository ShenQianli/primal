# `PrimalResult.coef`

## Usage

```python
result.coef(n=None) -> dict
```

## Description

Extract and display the estimated coefficients corresponding to a specific
regularization parameter. This is the Python equivalent of R's
`coef.primal()`.

## Arguments

- `n`: `int` or `None` -- the index of the regularization parameter along the solution path (1-based, matching R convention). If `None`, the last iteration is used.

## Returns

`dict` with keys:

- `"lambda"`: the regularization parameter at iteration `n`
- `"df"`: degrees of freedom at iteration `n`
- `"beta"`: coefficient vector at iteration `n`
- `"beta0"`: intercept value at iteration `n` (only for SparseSVM)

## Notes

- The index `n` is **1-based** to match R's convention.
- Raises `IndexError` if `n` is out of range.

## See also

[`PrimalResult`](PrimalResult.md),
[`PrimalResult.summary`](summary.md),
[`dantzig_solver`](dantzig_solver.md),
[`sparse_svm_solver`](sparse_svm_solver.md)
