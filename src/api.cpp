#include "PSM/api.hpp"

static void load_input(int *_n, int *_d, double *_X, double *_y,
                       Eigen::MatrixXd& X, Eigen::VectorXd& y){
	int n = *_n;
	int d = *_d;
	if(n <= 0 || d <= 0) return;
	X.resize(n, d);
	y.resize(n);
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X_map(_X, n, d);
	X = X_map;
	Eigen::Map<Eigen::VectorXd> y_map(_y, n);
	y = y_map;
}

static void write_output(const PSMresult& result, int d, int *T,
                         double *lambda_list, double *x_list, double *y_list){
	Eigen::VectorXd x;
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

extern "C" void QuantileRegression_api
(	int *_n,
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
	if(n <= 0 || d <= 0) return;
	Eigen::MatrixXd X;
	Eigen::VectorXd y;
	double tau = *_tau;
	load_input(_n, _d, _X, _y, X, y);

	Eigen::MatrixXd A(n, 2*d+2*n);
	Eigen::VectorXd b(n);
	Eigen::VectorXd b_bar(n);
	Eigen::VectorXd c(2*d+2*n);
	Eigen::VectorXd c_bar(2*d+2*n);

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
	write_output(result, d, T, lambda_list, x_list, y_list);
}

extern "C" void SparseSVM_api
(	int *_n,
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
	if(n <= 0 || d <= 0) return;
	Eigen::MatrixXd X;
	Eigen::VectorXd y;
	Eigen::VectorXd x;
	load_input(_n, _d, _X, _y, X, y);

	Eigen::MatrixXd A(n+1, 2*n+2*d+3);
	Eigen::VectorXd b(n+1);
	Eigen::VectorXd b_bar(n+1);
	Eigen::VectorXd c(2*n+2*d+3);
	Eigen::VectorXd c_bar(2*n+2*d+3);

	A.setZero();
	A.block(0, 0, 1, 2*d).setOnes();
	Eigen::MatrixXd Z(n, d);
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

extern "C" void Dantzig_api
(	int *_n,
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
	if(n <= 0 || d <= 0) return;
	Eigen::MatrixXd X;
	Eigen::VectorXd y;
	load_input(_n, _d, _X, _y, X, y);

	Eigen::MatrixXd A(2*d, 4*d);
	Eigen::VectorXd b(2*d);
	Eigen::VectorXd b_bar(2*d);
	Eigen::VectorXd c(4*d);
	Eigen::VectorXd c_bar(4*d);

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
	write_output(result, d, T, lambda_list, x_list, y_list);
}

extern "C" void CompressedSensing_api
(
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
	if(n <= 0 || d <= 0) return;
	Eigen::MatrixXd X;
	Eigen::VectorXd y;
	load_input(_n, _d, _X, _y, X, y);

	Eigen::MatrixXd A(2*n, 2*d+2*n);
	Eigen::VectorXd b(2*n);
	Eigen::VectorXd b_bar(2*n);
	Eigen::VectorXd c(2*d+2*n);
	Eigen::VectorXd c_bar(2*d+2*n);

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
	write_output(result, d, T, lambda_list, x_list, y_list);
}
