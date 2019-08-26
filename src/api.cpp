#include "PSM/api.hpp"

extern "C" void QuantileRegression_api(
	int *_n,
	int *_d,
	double *_X,
	double *_y,
 	double *_tau,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
 ){
	int n = *_n;
	int d = *_d;
	MatrixXd X(n, d);
	VectorXd y(n);
	double tau = *_tau;
	VectorXd x;
	
	MatrixXd A(n, 2*d+2*n);
	VectorXd b(n);
	VectorXd b_bar(n);
	VectorXd c(2*d+2*n);
	VectorXd c_bar(2*d+2*n);
	
	for(int row = 0; row < n; ++row){
		for(int col = 0; col < d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	
	A.block(0, 0, n, d) = X;
	A.block(0, d, n, d) = -X;
	A.block(0, 2*d+n, n, n).setIdentity();
	A.block(0, 2*d, n, n) = -A.block(0, 2*d+n, n, n);
	
	b = y;
	b_bar.setZero();
	
	c.setZero();
	c.tail(2*n).setOnes();
	c.tail(2*n) = (1-tau) * c.tail(2*n);
	c.tail(n).setOnes();
	c.tail(n) = tau * c.tail(n);
	c = -c;
	
	c_bar.setZero();
	c_bar.head(2*d).setOnes();
	c_bar = -c_bar;
	
	PSM psm(A, b, b_bar, c, c_bar);
	
	PSMresult result = psm.solve(*max_it, *lambda_threshold);
	x.resize(result.d);
	*T = result.T;
	for(int t = 0; t < *T; ++t){
		lambda_list[t] = result.lambda_list[t];
		x = result.x_list.col(t);
		for(int i = 0; i < d; ++i){
			x_list[t*d+i] = x[i] - x[i+d];
		}
		y_list[t] = -result.y_list[t];
	}
}

extern "C" void SparseSVM_api(
	int *_n,
	int *_d,
	double *_X,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list,
	double *x0_list
 ){
	int n = *_n;
	int d = *_d;
	MatrixXd X(n, d);
	VectorXd y(n);
	VectorXd x;
	
	MatrixXd A(n+1, 2*n+2*d+3);
	VectorXd b(n+1);
	VectorXd b_bar(n+1);
	VectorXd c(2*n+2*d+3);
	VectorXd c_bar(2*n+2*d+3);
	
	for(int row = 0; row < n; ++row){
		for(int col = 0; col < d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	
	A.setZero();
	A.block(0, 0, 1, 2*d).setOnes();
	MatrixXd Z(n, d);
	for(int i = 0; i < n; ++i){
		Z.row(i) = y(i)*X.row(i);
	}
	A.block(1, 0, n, d) = Z;
	A.block(1, d, n, d) = -Z;
	A.block(1, 2*d, n, 1) = y;
	A.block(1, 2*d+1, n, 1) = -y;
	A.block(1, 2*d+2, n, n).setIdentity();
	A.block(1, 2*d+2, n, n) = -A.block(1, 2*d+2, n, n);
	A.block(0, 2*d+n+2, n+1, n+1).setIdentity();
	
	b.setOnes();
	b(0) = 0;
	
	b_bar.setZero();
	b_bar(0) = 1;
	
	c.setZero();
	c.tail(n).setOnes();
	c = -c;
	
	c_bar.setZero();
	
	PSM psm(A, b, b_bar, c, c_bar);
	PSMresult result = psm.solve(*max_it, *lambda_threshold);
	x.resize(result.d);
	*T = result.T;
	for(int t = 0; t < *T; ++t){
		lambda_list[t] = result.lambda_list[t];
		x = result.x_list.col(t);
		for(int i = 0; i < d; ++i){
			x_list[t*d+i] = x[i] - x[i+d];
		}
		y_list[t] = -result.y_list[t];
		x0_list[t] = x[2*d] - x[2*d+1];
	}
}

extern "C" void Dantzig_api(
	int *_n,
	int *_d,
	double *_X,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
 ){
	int n = *_n;
	int d = *_d;
	MatrixXd X(n, d);
	VectorXd y(n);
	VectorXd x;
	
	MatrixXd A(2*d, 4*d);
	VectorXd b(2*d);
	VectorXd b_bar(2*d);
	VectorXd c(4*d);
	VectorXd c_bar(4*d);
	
	for(int row = 0; row < n; ++row){
		for(int col = 0; col < d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	
	A.block(0, 0, d, d) = X.transpose() * X;
	A.block(0, d, d, d) = - A.block(0, 0, d, d);
	A.block(d, 0, d, d) = - A.block(0, 0, d, d);
	A.block(d, d, d, d) = A.block(0, 0, d, d);
	A.block(0, 2*d, 2*d, 2*d).setIdentity();
	
	b.head(d) = X.transpose() * y;
	b.tail(d) = -b.head(d);
	
	b_bar.setOnes();
	
	c.setZero();
	c.head(2*d).setOnes();
	c = -c;
	
	c_bar.setZero();
	PSM psm(A, b, b_bar, c, c_bar);
	PSMresult result = psm.solve(*max_it, *lambda_threshold);
	x.resize(result.d);
	*T = result.T;
	for(int t = 0; t < *T; ++t){
		lambda_list[t] = result.lambda_list[t];
		x = result.x_list.col(t);
		for(int i = 0; i < d; ++i){
			x_list[t*d+i] = x[i] - x[i+d];
		}
		y_list[t] = -result.y_list[t];
	}
}

extern "C" void CompressedSensing_api(
	int *_n,
	int *_d,
	double *_X,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
 ){
	int n = *_n;
	int d = *_d;
	MatrixXd X(n, d);
	VectorXd y(n);
	VectorXd x;
	
	MatrixXd A(2*n, 2*d+2*n);
	VectorXd b(2*n);
	VectorXd b_bar(2*n);
	VectorXd c(2*d+2*n);
	VectorXd c_bar(2*d+2*n);
	
	for(int row = 0; row < n; ++row){
		for(int col = 0; col < d; ++col){
			X(row, col) = _X[row * d + col];
		}
		y(row) = _y[row];
	}
	
	A.block(0, 0, n, d) = X;
	A.block(0, d, n, d) = - X;
	A.block(n, 0, n, d) = - X;
	A.block(n, d, n, d) = X;
	A.block(0, 2*d, 2*n, 2*n).setIdentity();
	
	b.head(n) = y;
	b.tail(n) = - y;
	
	b_bar.setOnes();
	
	c.setZero();
	c.head(2*d).setOnes();
	c = -c;
	
	c_bar.setZero();
	
	PSM psm(A, b, b_bar, c, c_bar);
	PSMresult result = psm.solve(*max_it, *lambda_threshold);
	x.resize(result.d);
	*T = result.T;
	for(int t = 0; t < *T; ++t){
		lambda_list[t] = result.lambda_list[t];
		x = result.x_list.col(t);
		for(int i = 0; i < d; ++i){
			x_list[t*d+i] = x[i] - x[i+d];
		}
		y_list[t] = -result.y_list[t];
	}
}

extern "C" void ParametricSimplexMethod_api(
	int *_M,
	int *_N,
	double *_A,
	double *_b,
	double *_b_bar,
	double *_c,
	double *_c_bar,
	int *B_init,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
 ){
	int M = *_M;
	int N = *_N;
	int d = N;
	MatrixXd A(M, N);
	VectorXd b(M);
	VectorXd b_bar(M);
	VectorXd c(N);
	VectorXd c_bar(N);
	VectorXd x;
	
	for(int i = 0; i < M; ++i){
		for(int j = 0; j < N; ++j){
			A(i, j) = _A[i*N+j];
		}
		b(i) = _b[i];
		b_bar(i) = _b_bar[i];
	}
	for(int i = 0; i < N; ++i){
		c(i) = _c[i];
		c_bar(i) = _c_bar[i];
	}
	
	PSM psm(A, b, b_bar, c, c_bar);
	PSMresult result = psm.solve(*max_it, *lambda_threshold);
	
	x.resize(result.d);
	*T = result.T;
	for(int t = 0; t < *T; ++t){
		lambda_list[t] = result.lambda_list[t];
		x = result.x_list.col(t);
		for(int i = 0; i < d; ++i){
			x_list[t*d+i] = x[i];
		}
		y_list[t] = result.y_list[t];
	}
}
