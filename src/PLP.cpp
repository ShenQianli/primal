#include "PSM/PLP.hpp"

PLP::PLP(const MatrixXd& _A,	/*parameter*/
		 const VectorXd& _b,	/*parameter*/
		 const VectorXd& _b_bar,	/*parameter*/
		 const VectorXd& _c,	/*parameter*/
		 const VectorXd& _c_bar	/*parameter*/
){
	int m = (int)_A.rows();
	int n = (int)_A.cols();
	M = m;
	N = n + m;
	vector<T> v; /* vector to store matrix as triplets */
	A.resize(M, N);
	b.resize(M);
	b_bar.resize(M);
	c.resize(N);
	c_bar.resize(M);
	
	A.block(0, 0, m, n) = _A;
	A.block(0, n, m, m).setIdentity(); //slack
	
	for(int i = 0; i < m; ++i){
		b(i) = _b(i);
		b_bar(i) = _b_bar(i);
	}
	
	c.setZero();
	c_bar.setZero();
	for(int i = 0; i < n; ++i){
		c(i) = _c(i);
		c_bar(i) = _c_bar(i);
	}
}
