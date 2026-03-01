# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PRIMAL (PaRametric sImplex Method for spArse Learning) implements parametric simplex method solvers for sparse learning problems: Dantzig selector, sparse SVM, compressed sensing, and quantile regression. The core algorithm is C++ with Eigen, exposed via R package (`PRIMAL`) and Python package (`pyprimal`).

## Build Commands

```bash
# Build shared library (platform-specific: .dylib/.so/.dll)
make clean && make dylib

# Build static + shared library
make all

# Python: build dylib, copy to package, install editable
make Pyinstall
# Or manually:
make clean && make dylib
cd python-package && pip install -e ".[viz,test]"

# R: build and install
make Rinstall
# R: build tarball only
make Rbuild
# R: check package
make Rcheck
```

## Testing

```bash
# Python tests (pytest)
cd python-package && pytest tests/ -v

# Quick Python verification
python -c "import pyprimal; pyprimal.test()"

# Run a single Python test
cd python-package && pytest tests/test_solvers.py::TestDantzig::test_basic -v

# R package check
make Rcheck
```

## Architecture

**Three-layer design:** C++ core → C API → language bindings (R/Python)

### C++ Core (`src/`, `include/PSM/`)
- `PSM.hpp` / `PSM.cpp` — PSM solver class using Eigen sparse matrices and LU factorization. Methods: `init()`, `solve(max_it, lambda_threshold)`, `lusolve_update_dxb()`, `lusolve_update_dzn()`.
- `api.hpp` / `api.cpp` — Four `extern "C"` functions (`Dantzig_api`, `SparseSVM_api`, `CompressedSensing_api`, `QuantileRegression_api`) that formulate each problem as an LP and call PSM.

### Python Package (`python-package/pyprimal/`)
- `core.py` — Solver functions (`dantzig_solver`, `sparse_svm_solver`, `compressed_sensing_solver`, `quantile_regression_solver`) that call the C API via `ctypes`. Returns `PrimalResult` dataclass with `.summary()`, `.coef(n)`, `.plot(n)`.
- `libpath.py` — Discovers platform-specific shared library from `pyprimal/lib/`.
- `setup.py` — Custom `build_py` that compiles C++ via `make dylib` if no pre-built library exists.

### R Package (`R-package/`)
- Solver functions in `R/` call compiled C code via `.C()` interface. Returns S3 `"primal"` objects with `print`, `coef`, `plot` methods.

## Compiler Configuration

- C++14 standard. macOS: clang/clang++. Linux: gcc/g++.
- Flags: `-std=c++14 -Wall -fPIC -I./include -I./include/eigen`
- Eigen headers are vendored at `include/eigen/` (not a system dependency).

## CI/CD

GitHub Actions (`.github/workflows/build-publish.yml`) builds wheels via `cibuildwheel` on Ubuntu and macOS for Python 3.9–3.13. Publishes to PyPI on version tags (`v*`).
