"""PRIMAL: Parametric Simplex Method for Sparse Learning

Function-based API with PrimalResult dataclass, matching the R package.
"""

import time
from dataclasses import dataclass
from typing import Optional

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

from .libpath import find_lib_path

# ---------------------------------------------------------------------------
# Load shared library
# ---------------------------------------------------------------------------


def _load_lib():
    """Load the libpsm shared library."""
    lib_path = find_lib_path()
    return ctypes.cdll.LoadLibrary(lib_path[0])


_LIB = _load_lib()

# ctypes helpers
_CDouble = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
_CInt = ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")

# ---------------------------------------------------------------------------
# PrimalResult dataclass
# ---------------------------------------------------------------------------


@dataclass
class PrimalResult:
    """Result object returned by all PRIMAL solvers.

    Attributes
    ----------
    type : str
        Problem type (e.g., ``"Dantzig"``, ``"SparseSVM"``).
    data : ndarray of shape (n, d)
        Input data matrix.
    response : ndarray of shape (n,)
        Input response vector.
    beta : ndarray of shape (d, iterN)
        Coefficient estimates; each column corresponds to a regularization
        parameter.
    beta0 : ndarray of shape (iterN,) or None
        Intercept estimates (only for SparseSVM).
    df : ndarray of shape (iterN,)
        Degrees of freedom (number of nonzero coefficients) along the path.
    value : ndarray of shape (iterN,)
        Optimal objective function values along the path.
    iterN : int
        Number of iterations performed.
    lambda_ : ndarray of shape (iterN,)
        Regularization parameters along the solution path.
    runtime : float
        Elapsed time in seconds.
    """

    type: str
    data: np.ndarray
    response: np.ndarray
    beta: np.ndarray
    beta0: Optional[np.ndarray]
    df: np.ndarray
    value: np.ndarray
    iterN: int
    lambda_: np.ndarray
    runtime: float

    # ---- print / summary -------------------------------------------------

    def summary(self) -> str:
        """Print summary information including the model type, regularization
        parameter path, and degrees of freedom along the solution path.

        Returns
        -------
        str
            Multi-line summary of the result.
        """
        lines = [
            "",
            f"*******Parametric Simplex Method solving "
            f"{self.type} problem*********",
            f"iteration times =  {self.iterN}",
            "lambda list:",
            "  ".join(f"{v:.5g}" for v in self.lambda_),
            f"Degree of freedom: {self.df[0]} -----> {self.df[-1]}",
            f"Runtime: {self.runtime:.4f}  secs",
        ]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        return (
            f"PrimalResult(type={self.type!r}, iterN={self.iterN}, "
            f"d={self.beta.shape[0]}, runtime={self.runtime:.4f}s)"
        )

    # ---- coef ------------------------------------------------------------

    def coef(self, n: Optional[int] = None) -> dict:
        """Extract and display the estimated coefficients corresponding to
        a specific regularization parameter.

        Parameters
        ----------
        n : int or None
            The index of the regularization parameter along the solution
            path. If ``None``, the last iteration is used.

        Returns
        -------
        dict
            Dictionary with keys ``"lambda"``, ``"df"``, ``"beta"``,
            and optionally ``"beta0"`` (for SparseSVM).
        """
        if n is None:
            n = self.iterN
        if n < 1 or n > self.iterN:
            raise IndexError(f"n must be between 1 and {self.iterN}, got {n}")
        idx = n - 1  # convert to 0-based
        print(f"{self.type}  problem")
        print(f"index:  {n}")
        print(f"lambda:  {self.lambda_[idx]:.5g}")
        print(f"degree of freedom:  {self.df[idx]}")
        print("beta:")
        print(self.beta[:, idx])
        result = {
            "lambda": self.lambda_[idx],
            "df": self.df[idx],
            "beta": self.beta[:, idx].copy(),
        }
        if self.beta0 is not None:
            print(f"beta0:  {self.beta0[idx]:.5g}")
            result["beta0"] = self.beta0[idx]
        return result

    # ---- plot ------------------------------------------------------------

    def plot(self, n: Optional[int] = None) -> None:
        """Plot the regularization path and parameters obtained from the
        parametric simplex method.

        Parameters
        ----------
        n : int or None
            If ``None``, all three plots are shown together. If ``n`` is
            a number (1, 2, or 3), only the corresponding plot is shown:

            1. Coefficients vs regularization parameter
            2. Coefficients vs iteration
            3. Lambda vs iteration
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plotting. "
                "Install it with: pip install pyprimal[viz]"
            )

        iters = np.arange(1, self.iterN + 1)

        if n is None:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))
            # Plot 1: Coef vs Lambda
            axes[0].plot(self.lambda_, self.beta.T)
            axes[0].set_xlabel("Regularization Parameter")
            axes[0].set_ylabel("Coefficient")
            axes[0].set_title("Regularization Path")
            # Plot 2: Coef vs Iteration
            axes[1].plot(iters, self.beta.T)
            axes[1].set_xlabel("Iteration")
            axes[1].set_ylabel("Coefficient")
            axes[1].set_title("Regularization Path")
            # Plot 3: Lambda vs Iteration
            axes[2].plot(iters, self.lambda_)
            axes[2].set_xlabel("Iteration")
            axes[2].set_ylabel("Lambda")
            axes[2].set_title("Value of Lambda along the Path")
            plt.tight_layout()
            plt.show()
        elif n == 1:
            plt.figure()
            plt.plot(self.lambda_, self.beta.T)
            plt.xlabel("Regularization Parameter")
            plt.ylabel("Coefficient")
            plt.title("Regularization Path")
            plt.show()
        elif n == 2:
            plt.figure()
            plt.plot(iters, self.beta.T)
            plt.xlabel("Iteration")
            plt.ylabel("Coefficient")
            plt.title("Regularization Path")
            plt.show()
        elif n == 3:
            plt.figure()
            plt.plot(iters, self.lambda_)
            plt.xlabel("Iteration")
            plt.ylabel("Lambda")
            plt.title("Value of Lambda along the Path")
            plt.show()
        else:
            raise ValueError(
                f"n must be 1, 2, or 3 (or None for all), got {n}"
            )


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def _validate_inputs(X, y):
    """Validate and convert inputs to contiguous double arrays."""
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"X must be a 2D array, got {X.ndim}D")
    if y.ndim != 1:
        raise ValueError(f"y must be a 1D array, got {y.ndim}D")
    if X.shape[0] != y.shape[0]:
        raise ValueError(
            f"X and y must have the same number of rows: "
            f"X has {X.shape[0]}, y has {y.shape[0]}"
        )
    return X, y


# ---------------------------------------------------------------------------
# Solver functions
# ---------------------------------------------------------------------------


def dantzig_solver(X, y, max_it=50, lambda_threshold=0.01):
    """Solve a Dantzig selector problem using the parametric simplex method.

    Parameters
    ----------
    X : array-like of shape (n, d)
        An ``n`` by ``d`` data matrix.
    y : array-like of shape (n,)
        A length ``n`` response vector.
    max_it : int, default=50
        Maximum number of iterations for the solution path.
        The default is ``50``.
    lambda_threshold : float, default=0.01
        The algorithm stops when the regularization parameter falls below
        this threshold. The default value is ``0.01``.

    Returns
    -------
    PrimalResult
        An object containing the following attributes:

        - **data** -- The ``n`` by ``d`` data matrix from the input.
        - **response** -- The length ``n`` response vector from the input.
        - **beta** -- A matrix of coefficient estimates whose columns
          correspond to regularization parameters along the solution path.
        - **df** -- The degrees of freedom (number of nonzero coefficients)
          along the solution path.
        - **value** -- The sequence of optimal objective function values
          corresponding to each ``lambda``.
        - **iterN** -- The number of iterations performed.
        - **lambda_** -- The sequence of regularization parameters obtained
          along the solution path.
        - **type** -- The type of the problem, e.g., ``"Dantzig"``
          or ``"SparseSVM"``.
        - **runtime** -- Elapsed time in seconds.

    See Also
    --------
    sparse_svm_solver, compressed_sensing_solver,
    quantile_regression_solver

    Examples
    --------
    >>> import numpy as np
    >>> from pyprimal import dantzig_solver
    >>> rng = np.random.default_rng(42)
    >>> X = rng.standard_normal((100, 20))
    >>> beta_true = np.array([1]*5 + [0]*15, dtype=float)
    >>> y = X @ beta_true + 0.1 * rng.standard_normal(100)
    >>> result = dantzig_solver(X, y)
    >>> print(result)
    """
    X, y = _validate_inputs(X, y)
    n, d = X.shape

    _func = _LIB.Dantzig_api
    _func.argtypes = [
        _CInt, _CInt, _CDouble, _CDouble,
        _CInt, _CDouble,
        _CInt, _CDouble, _CDouble, _CDouble,
    ]

    iter_count = np.zeros(1, dtype=np.int32)
    lambda_list = np.zeros(max_it, dtype=np.float64)
    x_list = np.zeros(max_it * d, dtype=np.float64)
    y_list = np.zeros(max_it, dtype=np.float64)

    t0 = time.time()
    _func(
        np.array([n], dtype=np.int32),
        np.array([d], dtype=np.int32),
        X, y,
        np.array([max_it], dtype=np.int32),
        np.array([lambda_threshold], dtype=np.float64),
        iter_count, lambda_list, x_list, y_list,
    )
    runtime = time.time() - t0

    T = int(iter_count[0])
    beta = x_list[: T * d].reshape(T, d).T  # (d, T)
    df = np.count_nonzero(beta, axis=0)

    return PrimalResult(
        type="Dantzig",
        data=X,
        response=y,
        beta=beta,
        beta0=None,
        df=df,
        value=y_list[:T].copy(),
        iterN=T,
        lambda_=lambda_list[:T].copy(),
        runtime=runtime,
    )


def sparse_svm_solver(X, y, max_it=50, lambda_threshold=0.01):
    """Solve a sparse SVM problem using the parametric simplex method.

    Parameters
    ----------
    X : array-like of shape (n, d)
        An ``n`` by ``d`` data matrix.
    y : array-like of shape (n,)
        A length ``n`` response vector with values in {-1, 1}.
    max_it : int, default=50
        Maximum number of iterations for the solution path.
        The default is ``50``.
    lambda_threshold : float, default=0.01
        The algorithm stops when the regularization parameter falls below
        this threshold. The default value is ``0.01``.

    Returns
    -------
    PrimalResult
        An object containing the following attributes:

        - **data** -- The ``n`` by ``d`` data matrix from the input.
        - **response** -- The length ``n`` response vector from the input.
        - **beta** -- A matrix of coefficient estimates whose columns
          correspond to regularization parameters along the solution path.
        - **beta0** -- A vector of intercept estimates corresponding to
          regularization parameters along the solution path.
        - **df** -- The degrees of freedom (number of nonzero coefficients)
          along the solution path.
        - **value** -- The sequence of optimal objective function values
          corresponding to each ``lambda``.
        - **iterN** -- The number of iterations performed.
        - **lambda_** -- The sequence of regularization parameters obtained
          along the solution path.
        - **type** -- The type of the problem, e.g., ``"Dantzig"``
          or ``"SparseSVM"``.
        - **runtime** -- Elapsed time in seconds.

    See Also
    --------
    dantzig_solver, compressed_sensing_solver,
    quantile_regression_solver

    Examples
    --------
    >>> import numpy as np
    >>> from pyprimal import sparse_svm_solver
    >>> rng = np.random.default_rng(42)
    >>> X = rng.standard_normal((100, 20))
    >>> y = np.where(X[:, 0] + X[:, 1] > 0, 1.0, -1.0)
    >>> result = sparse_svm_solver(X, y)
    >>> print(result)
    """
    X, y = _validate_inputs(X, y)
    n, d = X.shape

    _func = _LIB.SparseSVM_api
    _func.argtypes = [
        _CInt, _CInt, _CDouble, _CDouble,
        _CInt, _CDouble,
        _CInt, _CDouble, _CDouble, _CDouble, _CDouble,
    ]

    iter_count = np.zeros(1, dtype=np.int32)
    lambda_list = np.zeros(max_it, dtype=np.float64)
    x_list = np.zeros(max_it * d, dtype=np.float64)
    y_list = np.zeros(max_it, dtype=np.float64)
    x0_list = np.zeros(max_it, dtype=np.float64)

    t0 = time.time()
    _func(
        np.array([n], dtype=np.int32),
        np.array([d], dtype=np.int32),
        X, y,
        np.array([max_it], dtype=np.int32),
        np.array([lambda_threshold], dtype=np.float64),
        iter_count, lambda_list, x_list, y_list, x0_list,
    )
    runtime = time.time() - t0

    T = int(iter_count[0])
    beta = x_list[: T * d].reshape(T, d).T  # (d, T)
    df = np.count_nonzero(beta, axis=0)

    return PrimalResult(
        type="SparseSVM",
        data=X,
        response=y,
        beta=beta,
        beta0=x0_list[:T].copy(),
        df=df,
        value=y_list[:T].copy(),
        iterN=T,
        lambda_=lambda_list[:T].copy(),
        runtime=runtime,
    )


def compressed_sensing_solver(X, y, max_it=50, lambda_threshold=0.01):
    """Solve a compressed sensing problem using the parametric simplex method.

    Parameters
    ----------
    X : array-like of shape (n, d)
        An ``n`` by ``d`` data matrix.
    y : array-like of shape (n,)
        A length ``n`` response vector.
    max_it : int, default=50
        Maximum number of iterations for the solution path.
        The default is ``50``.
    lambda_threshold : float, default=0.01
        The algorithm stops when the regularization parameter falls below
        this threshold. The default value is ``0.01``.

    Returns
    -------
    PrimalResult
        The returned value is similar to that of :func:`dantzig_solver`.

    See Also
    --------
    dantzig_solver, sparse_svm_solver, quantile_regression_solver

    Examples
    --------
    >>> import numpy as np
    >>> from pyprimal import compressed_sensing_solver
    >>> rng = np.random.default_rng(42)
    >>> X = rng.standard_normal((50, 100))
    >>> beta_true = np.zeros(100); beta_true[:5] = 1.0
    >>> y = X @ beta_true + 0.1 * rng.standard_normal(50)
    >>> result = compressed_sensing_solver(X, y)
    >>> print(result)
    """
    X, y = _validate_inputs(X, y)
    n, d = X.shape

    _func = _LIB.CompressedSensing_api
    _func.argtypes = [
        _CInt, _CInt, _CDouble, _CDouble,
        _CInt, _CDouble,
        _CInt, _CDouble, _CDouble, _CDouble,
    ]

    iter_count = np.zeros(1, dtype=np.int32)
    lambda_list = np.zeros(max_it, dtype=np.float64)
    x_list = np.zeros(max_it * d, dtype=np.float64)
    y_list = np.zeros(max_it, dtype=np.float64)

    t0 = time.time()
    _func(
        np.array([n], dtype=np.int32),
        np.array([d], dtype=np.int32),
        X, y,
        np.array([max_it], dtype=np.int32),
        np.array([lambda_threshold], dtype=np.float64),
        iter_count, lambda_list, x_list, y_list,
    )
    runtime = time.time() - t0

    T = int(iter_count[0])
    beta = x_list[: T * d].reshape(T, d).T  # (d, T)
    df = np.count_nonzero(beta, axis=0)

    return PrimalResult(
        type="Compressed sensing",
        data=X,
        response=y,
        beta=beta,
        beta0=None,
        df=df,
        value=y_list[:T].copy(),
        iterN=T,
        lambda_=lambda_list[:T].copy(),
        runtime=runtime,
    )


def quantile_regression_solver(X, y, max_it=50, lambda_threshold=0.01,
                               tau=0.5):
    """Solve a quantile regression problem using the parametric simplex method.

    Parameters
    ----------
    X : array-like of shape (n, d)
        An ``n`` by ``d`` data matrix.
    y : array-like of shape (n,)
        A length ``n`` response vector.
    max_it : int, default=50
        Maximum number of iterations for the solution path.
        The default is ``50``.
    lambda_threshold : float, default=0.01
        The algorithm stops when the regularization parameter falls below
        this threshold. The default value is ``0.01``.
    tau : float, default=0.5
        The quantile level for the regression, must be in (0, 1).
        The default is ``0.5``.

    Returns
    -------
    PrimalResult
        The returned value is similar to that of :func:`dantzig_solver`.

    See Also
    --------
    dantzig_solver, sparse_svm_solver, compressed_sensing_solver

    Examples
    --------
    >>> import numpy as np
    >>> from pyprimal import quantile_regression_solver
    >>> rng = np.random.default_rng(42)
    >>> X = rng.standard_normal((100, 20))
    >>> beta_true = np.array([1]*5 + [0]*15, dtype=float)
    >>> y = X @ beta_true + rng.standard_normal(100)
    >>> result = quantile_regression_solver(X, y, tau=0.5)
    >>> print(result)
    """
    if not 0 < tau < 1:
        raise ValueError(f"tau must be in (0, 1), got {tau}")
    X, y = _validate_inputs(X, y)
    n, d = X.shape

    _func = _LIB.QuantileRegression_api
    _func.argtypes = [
        _CInt, _CInt, _CDouble, _CDouble,
        _CDouble,
        _CInt, _CDouble,
        _CInt, _CDouble, _CDouble, _CDouble,
    ]

    iter_count = np.zeros(1, dtype=np.int32)
    lambda_list = np.zeros(max_it, dtype=np.float64)
    x_list = np.zeros(max_it * d, dtype=np.float64)
    y_list = np.zeros(max_it, dtype=np.float64)

    t0 = time.time()
    _func(
        np.array([n], dtype=np.int32),
        np.array([d], dtype=np.int32),
        X, y,
        np.array([tau], dtype=np.float64),
        np.array([max_it], dtype=np.int32),
        np.array([lambda_threshold], dtype=np.float64),
        iter_count, lambda_list, x_list, y_list,
    )
    runtime = time.time() - t0

    T = int(iter_count[0])
    beta = x_list[: T * d].reshape(T, d).T  # (d, T)
    df = np.count_nonzero(beta, axis=0)

    return PrimalResult(
        type="Quantile Regression",
        data=X,
        response=y,
        beta=beta,
        beta0=None,
        df=df,
        value=y_list[:T].copy(),
        iterN=T,
        lambda_=lambda_list[:T].copy(),
        runtime=runtime,
    )
