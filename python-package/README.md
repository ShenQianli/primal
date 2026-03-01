# pyprimal

**pyprimal**: Parametric Simplex Method for Sparse Learning

A native Python package implementing the parametric simplex method for a
variety of sparse learning problems including Dantzig selector, sparse SVM,
compressed sensing, and quantile regression. The core algorithm is implemented
in C++ with Eigen support; Python wraps it via ctypes with no R dependency.

## Directory structure

```
python-package/
├── pyprimal/           Package source
│   ├── core.py         Solvers and PrimalResult dataclass
│   ├── libpath.py      Shared library loader
│   └── lib/            Compiled libpsm (.so / .dylib)
├── docs/               MkDocs documentation pages
├── examples/           Runnable example scripts
├── tests/              Unit tests (pytest)
├── mkdocs.yml          Documentation configuration
└── pyproject.toml      Build configuration
```

## Requirements

- Python >= 3.9
- NumPy >= 1.23
- Compiled `libpsm` shared library (Linux `.so` or macOS `.dylib`)
- Optional: matplotlib >= 3.5 for plotting

## Installation

From source (recommended):

```bash
git clone https://github.com/Gatech-Flash/primal.git
cd primal

# Build the C++ library
make clean && make dylib

# Install the Python package
cd python-package
pip install -e ".[viz,test]"
```

Optional extras:

```bash
pip install -e ".[viz]"    # matplotlib for plotting
pip install -e ".[test]"   # pytest for testing
pip install -e ".[docs]"   # mkdocs for documentation
```

Verify installation:

```bash
python -c "import pyprimal; pyprimal.test()"
```

## Usage

```python
import numpy as np
from pyprimal import dantzig_solver, sparse_svm_solver

# Dantzig selector
rng = np.random.default_rng(42)
X = rng.standard_normal((100, 20))
beta_true = np.array([1]*5 + [0]*15, dtype=float)
y = X @ beta_true + 0.1 * rng.standard_normal(100)

result = dantzig_solver(X, y)
print(result)           # summary
result.coef()           # coefficients at last iteration
result.coef(3)          # coefficients at iteration 3
result.plot()           # all three plots
result.plot(1)          # regularization path only

# Sparse SVM
y_svm = np.where(X[:, 0] + X[:, 1] > 0, 1.0, -1.0)
result_svm = sparse_svm_solver(X, y_svm)
print(result_svm)
```

### Available solvers

| Function | Problem |
|----------|---------|
| `dantzig_solver(X, y)` | Dantzig selector |
| `sparse_svm_solver(X, y)` | Sparse support vector machine |
| `compressed_sensing_solver(X, y)` | Compressed sensing |
| `quantile_regression_solver(X, y, tau=0.5)` | Quantile regression |

All solvers return a `PrimalResult` dataclass with:

- `.beta` -- coefficient matrix (d x iterN)
- `.lambda_` -- regularization parameter path
- `.value` -- objective function values
- `.df` -- degrees of freedom along the path
- `.summary()` / `print(result)` -- formatted summary
- `.coef(n)` -- extract coefficients at iteration n (1-based index)
- `.plot(n)` -- visualize the solution path

## Documentation

Build the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs build
mkdocs serve          # opens at http://127.0.0.1:8000
```

## Developer workflow

```bash
# Run tests
pip install -e ".[test]"
pytest tests/ -v

# Build docs
pip install -e ".[docs]"
mkdocs build --strict
```

## Citation

```bibtex
@article{li2018primal,
  title={The Parametric Simplex Method for Sparse Learning},
  author={Li, Zichong and Shen, Qianli and Zhao, Tuo},
  year={2018}
}
```

## License

GPL-3.0

## Authors

Zichong Li, Qianli Shen, Tuo Zhao

Maintainer: Tuo Zhao <tourzhao@gatech.edu>
