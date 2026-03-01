"""PYPRIMAL: Parametric Simplex Method for Sparse Learning"""

import os

from .core import (
    PrimalResult,
    dantzig_solver,
    sparse_svm_solver,
    compressed_sensing_solver,
    quantile_regression_solver,
)

__all__ = [
    "PrimalResult",
    "dantzig_solver",
    "sparse_svm_solver",
    "compressed_sensing_solver",
    "quantile_regression_solver",
]

_VERSION_PATH = os.path.join(os.path.dirname(__file__), "VERSION")
with open(_VERSION_PATH) as _f:
    __version__ = _f.read().strip()


def test():
    """Verify that the package is installed and the C library is loaded."""
    print(f"pyprimal {__version__} loaded successfully!")
    print(f"Available solvers: {', '.join(__all__[1:])}")
