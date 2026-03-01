# `test`

## Usage

```python
pyprimal.test()
```

## Description

Verify that the `pyprimal` package is installed and the C library (`libpsm`)
is loaded correctly.

## Returns

Prints version and library status to stdout.

## Notes

- If the C library cannot be found, `test()` will raise an `OSError`.
