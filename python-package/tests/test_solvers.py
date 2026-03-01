"""Tests for pyprimal solver functions."""

import numpy as np
import pytest

from pyprimal import (
    PrimalResult,
    dantzig_solver,
    sparse_svm_solver,
    compressed_sensing_solver,
    quantile_regression_solver,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def regression_data():
    """Generate a simple regression dataset."""
    rng = np.random.default_rng(42)
    n, d = 80, 15
    X = rng.standard_normal((n, d))
    beta_true = np.zeros(d)
    beta_true[:4] = [2.0, -1.5, 1.0, -0.5]
    y = X @ beta_true + 0.1 * rng.standard_normal(n)
    return X, y


@pytest.fixture
def svm_data():
    """Generate a simple binary classification dataset."""
    rng = np.random.default_rng(42)
    n, d = 80, 15
    X = rng.standard_normal((n, d))
    y = np.where(X[:, 0] + X[:, 1] > 0, 1.0, -1.0)
    return X, y


# ---------------------------------------------------------------------------
# Structure tests
# ---------------------------------------------------------------------------


def _check_result_structure(result, X, y, expect_beta0=False):
    """Common assertions for all PrimalResult objects."""
    assert isinstance(result, PrimalResult)
    n, d = X.shape

    # Field types
    assert isinstance(result.type, str)
    assert isinstance(result.iterN, int)
    assert isinstance(result.runtime, float)

    # Shapes
    assert result.data.shape == (n, d)
    assert result.response.shape == (n,)
    assert result.beta.shape[0] == d
    assert result.beta.shape[1] == result.iterN
    assert result.df.shape == (result.iterN,)
    assert result.value.shape == (result.iterN,)
    assert result.lambda_.shape == (result.iterN,)

    # df consistency
    expected_df = np.count_nonzero(result.beta, axis=0)
    np.testing.assert_array_equal(result.df, expected_df)

    # beta0
    if expect_beta0:
        assert result.beta0 is not None
        assert result.beta0.shape == (result.iterN,)
    else:
        assert result.beta0 is None

    # iterN > 0
    assert result.iterN > 0
    assert result.runtime >= 0


class TestDantzig:
    def test_basic(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        _check_result_structure(result, X, y)
        assert result.type == "Dantzig"

    def test_custom_params(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y, max_it=20, lambda_threshold=0.1)
        _check_result_structure(result, X, y)
        assert result.iterN <= 20


class TestSparseSVM:
    def test_basic(self, svm_data):
        X, y = svm_data
        result = sparse_svm_solver(X, y)
        _check_result_structure(result, X, y, expect_beta0=True)
        assert result.type == "SparseSVM"


class TestCompressedSensing:
    def test_basic(self, regression_data):
        X, y = regression_data
        result = compressed_sensing_solver(X, y)
        _check_result_structure(result, X, y)
        assert result.type == "Compressed sensing"


class TestQuantileRegression:
    def test_basic(self, regression_data):
        X, y = regression_data
        result = quantile_regression_solver(X, y, tau=0.5)
        _check_result_structure(result, X, y)
        assert result.type == "Quantile Regression"

    def test_invalid_tau(self, regression_data):
        X, y = regression_data
        with pytest.raises(ValueError, match="tau must be in"):
            quantile_regression_solver(X, y, tau=0.0)
        with pytest.raises(ValueError, match="tau must be in"):
            quantile_regression_solver(X, y, tau=1.0)
        with pytest.raises(ValueError, match="tau must be in"):
            quantile_regression_solver(X, y, tau=-0.5)


# ---------------------------------------------------------------------------
# Input validation tests
# ---------------------------------------------------------------------------


class TestValidation:
    def test_x_not_2d(self):
        with pytest.raises(ValueError, match="2D"):
            dantzig_solver(np.ones(10), np.ones(10))

    def test_y_not_1d(self):
        with pytest.raises(ValueError, match="1D"):
            dantzig_solver(np.ones((5, 3)), np.ones((5, 1)))

    def test_shape_mismatch(self):
        with pytest.raises(ValueError, match="same number of rows"):
            dantzig_solver(np.ones((5, 3)), np.ones(4))


# ---------------------------------------------------------------------------
# summary / coef tests
# ---------------------------------------------------------------------------


class TestSummaryCoef:
    def test_summary(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        s = result.summary()
        assert "Dantzig" in s
        assert "iteration times" in s
        assert "lambda list" in s
        assert "Degree of freedom" in s
        assert "Runtime" in s

    def test_str(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        assert str(result) == result.summary()

    def test_coef_default(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        c = result.coef()
        assert "lambda" in c
        assert "df" in c
        assert "beta" in c
        assert c["beta"].shape == (regression_data[0].shape[1],)

    def test_coef_indexed(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        c = result.coef(1)
        np.testing.assert_array_equal(c["beta"], result.beta[:, 0])

    def test_coef_out_of_range(self, regression_data):
        X, y = regression_data
        result = dantzig_solver(X, y)
        with pytest.raises(IndexError):
            result.coef(0)
        with pytest.raises(IndexError):
            result.coef(result.iterN + 1)

    def test_coef_svm_has_beta0(self, svm_data):
        X, y = svm_data
        result = sparse_svm_solver(X, y)
        c = result.coef()
        assert "beta0" in c
