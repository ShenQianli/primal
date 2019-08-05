#include "PSM/problems.hpp"

PLP *QuantileRegression(const MatrixXd& X, const VectorXd& y){
	int n = (int)X.rows();
	int d = (int)X.cols();
	MatrixXd A(n, 2*d+2*n);
	VectorXd b(n);
	VectorXd b_bar(n);
	VectorXd c(2*d+2*n);
	VectorXd c_bar(2*d+2*n);
	
	A.block(0, 0, n, d) = X;
	A.block(0, d, n, d) = -X;
	A.block(0, 2*d, n, n).setIdentity();
	A.block(0, 2*d+n, n, n) = -A.block(0, 2*d, n, n);
	
	b = y;
	b_bar.setZero();
	
	c.setZero();
	c.tail(n).setOnes();
	c = -c;
	
	c_bar.setZero();
	c_bar.tail(2*n).setOnes();
	c_bar.tail(n) = -c_bar.tail(n);
	c_bar = -c_bar;
	
	return new PLP(A, b, b_bar, c, c_bar);
}

PLP *SparseSVM(const MatrixXd& X, const VectorXd& y){
	int n = (int)X.rows();
	int d = (int)X.cols();
	MatrixXd A(n+1, 2*n+3*d+2);
	VectorXd b(n+1);
	VectorXd b_bar(n+1);
	VectorXd c(2*n+3*d+2);
	VectorXd c_bar(2*n+3*d+2);
	
	A.setZero();
	A.block(0, 0, n, n).setIdentity();
	A.block(0, n, n, n) = -A.block(0, 0, n, n);
	MatrixXd Z(n, d);
	for(int i = 0; i < n; ++i){
		Z.row(i) = y(i) * X.row(i);
	}
	A.block(0, 2*n, n, d) = Z;
	A.block(0, 2*n+d, n, d) = -Z;
	A.block(n, 2*n, 1, 2*d).setOnes();
	A.block(0, 2*n+2*d, n, 1) = y;
	A.block(0, 2*n+2*d+1, n, 1) = -y;
	A.block(n, 2*n+2*d+2, 1, d).setOnes();
	
	b.setZero();
	b.head(n).setOnes();
	
	b_bar.setZero();
	b(n) = 1;
	
	c.setZero();
	c.head(n).setOnes();
	c.head(n) = -c.head(n);
	
	c_bar.setZero();
	
	return new PLP(A, b, b_bar, c, c_bar);
}


PLP *Dantzig(const MatrixXd& X, const VectorXd& y){
	//	int n = (int)X.rows();
	int d = (int)X.cols();
	MatrixXd A(2*d, 2*d);
	VectorXd b(2*d);
	VectorXd b_bar(2*d);
	VectorXd c(2*d);
	VectorXd c_bar(2*d);
	
	A.block(0, 0, d, d) = X.transpose() * X;
	A.block(0, d, d, d) = - A.block(0, 0, d, d);
	A.block(d, 0, d, d) = - A.block(0, 0, d, d);
	A.block(d, d, d, d) = A.block(0, 0, d, d);
	
	b.head(d) = X.transpose() * y;
	b.tail(d) = -b.head(d);
	
	b_bar.setOnes();
	
	c.setOnes();
	c = -c;
	
	c_bar.setZero();
	
	return new PLP(A, b, b_bar, c, c_bar);
}

PLP *CompressedSensing(const MatrixXd& X, const VectorXd& y){
	int n = (int)X.rows();
	int d = (int)X.cols();
	MatrixXd A(2*n, 2*d);
	VectorXd b(2*n);
	VectorXd b_bar(2*n);
	VectorXd c(2*d);
	VectorXd c_bar(2*d);
	
	A.block(0, 0, n, d) = X;
	A.block(0, d, n, d) = - X;
	A.block(n, 0, n, d) = - X;
	A.block(n, d, n, d) = X;
	
	b.head(n) = y;
	b.tail(n) = - y;
	
	b_bar.setOnes();
	
	c.setOnes();
	c = -c;
	
	c_bar.setZero();
	
	return new PLP(A, b, b_bar, c, c_bar);
}
