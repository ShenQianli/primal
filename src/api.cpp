#include "PSM/api.hpp"

PLP *LinearProgrammingWrapper(
	int *_m,	/*row num*/
	int *_n,	/*col num*/
	double *_A,	/* A of size m*n*/
	double *_b,	/* b of size m*/
	double *_b_bar,	/* b_bar of size m*/
	double *_c,	/* c of size n*/
	double *_c_bar	/* c_bar of size n*/
){
	int m = *_m;
	int n = *_n;
	MatrixXd A(m, n);
	VectorXd b(m);
	VectorXd b_bar(m);
	VectorXd c(n);
	VectorXd c_bar(n);
	
	for(int row = 0; row < m; ++row){
		for(int col = 0; col < n; ++col){
			A(row, col) = _A[row * n + col];
		}
	}
	
	for(int row = 0; row < m; ++row){
		b(row) = _b[row];
		b_bar(row) = _b_bar[row];
	}
	
	for(int col = 0; col < n; ++col){
		c(col) = _c[col];
		c_bar(col) = _c_bar[col];
	}
	
	return new PLP(A, b, b_bar, c, c_bar);
}

PLP *QuantileRegressionWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
){
	MatrixXd X(_n, _d);
	VectorXd y(_n);
	for(int row = 0; row < _n; ++row){
		for(int col = 0; col < _d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	return QuantileRegression(X, y);
}

PLP *SparseSVMWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
){
	MatrixXd X(_n, _d);
	VectorXd y(_n);
	for(int row = 0; row < _n; ++row){
		for(int col = 0; col < _d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	return SparseSVM(X, y);
}

PLP *DantzigWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
){
	MatrixXd X(_n, _d);
	VectorXd y(_n);
	for(int row = 0; row < _n; ++row){
		for(int col = 0; col < _d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	return Dantzig(X, y);
}

PLP *CompressedSensingWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
){
	MatrixXd X(_n, _d);
	VectorXd y(_n);
	for(int row = 0; row < _n; ++row){
		for(int col = 0; col < _d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	return CompressedSensing(X, y);
}

