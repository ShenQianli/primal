#ifndef PLP_HPP
#define PLP_HPP

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef Triplet<double> T;

class PLP{
	/*
	 parametric linear programming class
	 
	 primal:
	 min_x (c^T)x s.t. Ax <= b, x >= 0
	 after slack:
	 min_x (c^T)x s.t. Ax = b, x >= 0
	 add pertubations:
	 min_x ((c + lambda * c_bar)^T)dx s.t. Ax = b + lambda * b_bar, x>= 0
	 */
public:
	/*problem paramters*/
	int M; /*M = m*/
	int N; /*N = n+m*/
	MatrixXd A; /* A of size M*N */
	VectorXd b; /* b of size M*/
	VectorXd b_bar; /* b_bar of size M*/
	VectorXd c; /* c of size N*/
	VectorXd c_bar; /* c_bar of size N*/
	
	PLP(const MatrixXd& _A,	/*parameter*/
		const VectorXd& _b,	/*parameter*/
		const VectorXd& _b_bar,	/*parameter*/
		const VectorXd& _c,	/*parameter*/
		const VectorXd& _c_bar	/*parameter*/
	)
};

#endif
