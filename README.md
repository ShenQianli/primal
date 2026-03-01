<h1 align="center">PRIMAL</h1>
<h4 align="center">Parametric Simplex Method for Sparse Learning (R + Python)</h4>

___PRIMAL___ (PaRametric sImplex Method for spArse Learning) implements a unified framework of parametric simplex method for a variety of sparse learning problems (e.g., Dantzig selector, sparse quantile regression, sparse support vector machines, and compressed sensing) combined with efficient hyper-parameter selection strategies. The core algorithm is implemented in C++ with Eigen support for portable high performance linear algebra.

This repository provides two package variants:

- R package: `PRIMAL` (native R interface, available on CRAN)
- Python package: `pyprimal` (native Python interface via ctypes)

Both variants share the same C++ backend and produce numerically identical results.

## Package Variants

- **R version (`PRIMAL`)**: see [R Package](#r-package-primal).
- **Python version (`pyprimal`)**: see [Python Package](#python-package-pyprimal).

## Table of contents

- [Introduction](#introduction)
- [Directory structure](#directory-structure)
- [R Package (`PRIMAL`)](#r-package-primal)
- [Python Package (`pyprimal`)](#python-package-pyprimal)
- [Performance](#performance)
- [References](#references)

## Introduction

Linear Programming (LP) based sparse learning methods, such as the Dantzig selector (for linear regression) [1], sparse quantile regression [2], sparse support vector machines [3], have been widely used in machine learning for high dimensional data analysis [4, 5]. Despite their popularity, their software implementations are quite limited. This library is proposed for the aforementioned LP-based sparse learning methods. It has the following two key features: 1) It provides a highly efficient optimization engine based on the parametric simplex method [6], which can efficiently solve large scale sparse learning problems; 2) Besides the estimation procedures, it provides additional functional modules such as data-dependent model selection and model visualization.

We also provide tutorials on the theoretical background and the code. Please see the `tutorials` folder.

## Directory structure

```
primal/
├── src/                    C++ implementation of the PSM algorithm
│   ├── api.cpp             C API as interface for R and Python
│   └── PSM.cpp             Core PSM solver
├── include/
│   ├── PSM/                C++ header files
│   └── eigen/              Eigen header files for linear algebra
├── R-package/              R package (PRIMAL)
├── python-package/         Python package (pyprimal)
│   ├── pyprimal/           Package source
│   ├── docs/               MkDocs documentation
│   ├── examples/           Runnable example scripts
│   └── tests/              Unit tests
├── tutorials/              Tutorials for R and Python
├── profiling/              Performance benchmarks
├── Makefile                Build configuration
└── CMakeLists.txt          CMake build configuration
```

## R Package (`PRIMAL`)

### Installation

Install from CRAN (recommended):

```R
install.packages("PRIMAL")
```

Install from source:

```bash
git clone https://github.com/Gatech-Flash/primal.git
cd primal
make Rinstall
```

### R User Interface

```R
set.seed(42)
library(PRIMAL)

## Generate the design matrix and coefficient vector
n <- 100; d <- 20; s <- 5
X <- matrix(rnorm(n*d), n, d)
beta <- c(1, -1, 2, -2, 3, rep(0, d-s))
Y <- X %*% beta + 0.1 * rnorm(n)

## Dantzig selector
fit <- Dantzig_solver(X, Y)
print(fit)
coef(fit)
plot(fit)
```

### R API

| Function | Problem |
|----------|---------|
| `Dantzig_solver(X, y)` | Dantzig selector |
| `SparseSVM_solver(X, y)` | Sparse support vector machine |
| `CompressedSensing_solver(X, y)` | Compressed sensing |
| `QuantileRegression_solver(X, y, tau)` | Quantile regression |

All solvers return an S3 object of class `"primal"` with methods `print()`, `coef()`, and `plot()`.

## Python Package (`pyprimal`)

### Python package location

- `python-package/README.md`
- `python-package/docs/`
- `python-package/examples/`

### Python installation

```bash
git clone https://github.com/Gatech-Flash/primal.git
cd primal

# Build the C++ library
make clean && make dylib

# Install the Python package
cd python-package
pip install -e ".[viz,test]"

# Verify
python -c "import pyprimal; pyprimal.test()"
```

Optional extras:

```bash
pip install -e ".[viz]"    # matplotlib for plotting
pip install -e ".[test]"   # pytest for testing
pip install -e ".[docs]"   # mkdocs for documentation
```

### Python documentation

Build the documentation locally:

```bash
cd python-package
pip install -e ".[docs]"
mkdocs build
mkdocs serve          # opens at http://127.0.0.1:8000
```

### Python User Interface

```python
import numpy as np
from pyprimal import dantzig_solver

## Generate the design matrix and coefficient vector
np.random.seed(42)
n, d, s = 100, 20, 5
X = np.random.randn(n, d)
beta = np.array([1, -1, 2, -2, 3] + [0] * (d - s))
y = X @ beta + 0.1 * np.random.randn(n)

## Dantzig selector
result = dantzig_solver(X, y)
print(result)
result.coef()
result.plot()
```

### Python API coverage

| Function | Problem |
|----------|---------|
| `dantzig_solver(X, y)` | Dantzig selector |
| `sparse_svm_solver(X, y)` | Sparse support vector machine |
| `compressed_sensing_solver(X, y)` | Compressed sensing |
| `quantile_regression_solver(X, y, tau)` | Quantile regression |

All solvers return a `PrimalResult` dataclass with `.summary()`, `.coef(n)`, and `.plot(n)` methods.

## Performance

### R benchmark

Compares `PRIMAL` against the R package `fastclime`. Sample size n = 200; dimension d varies from 1000 to 7000 (Dantzig) and 200 to 1000 (compressed sensing).

```bash
cd profiling
Rscript benchmark.R
```

- Dantzig selector: PRIMAL achieves similar optimization performance to fastclime but is 6 to 10 times faster.
- Compressed sensing: PRIMAL is 2 times faster than fastclime and achieves similar optimization.

### Python benchmark

Benchmarks `pyprimal` solvers across the same dimension ranges.

```bash
cd profiling
python benchmark.py
```

The script measures CPU time for Dantzig selector (d = 1000–7000) and compressed sensing (d = 200–1000), and saves a plot to `profiling/images/performance_Python.png`.

## References

[1] Dantzig G. Linear programming and extensions, 2016.

[2] Belloni A, Chernozhukov V. L1-penalized quantile regression in high-dimensional sparse models, 2011.

[3] Candes E, Tao T. The Dantzig selector: Statistical estimation when p is much larger than n, 2007.

[4] Wang L. The L1 penalized LAD estimator for high dimensional linear regression, 2013.

[5] Bandyopadhyay S, Mehta M, Kuo D, et al. Rewiring of genetic networks in response to DNA damage, 2010.

[6] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.
