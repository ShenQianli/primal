#ifndef PSM_H
#define PSM_H

#include <iostream>
#include <cstring>
#include <cfloat>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

constexpr double EPS = 1e-10;

enum FLAGTYPE{NONE, PRIMAL, DUAL};

class PSMresult{
public:
	int T;
	int d;
	std::vector<double> lambda_list;
	Eigen::MatrixXd x_list;
	std::vector<double> y_list;
	PSMresult(int max_it, int _d);
	void update(double lambda, Eigen::VectorXd x, double y);
};

class PSM{
public:
	/*Parametric Linear Programming problem parameters*/
	int M; /*row num*/
	int N; /*col num*/
	int m; /*m = M*/
	int n; /*n = N-M*/
	const Eigen::MatrixXd& A; /* A of size M*N */
	const Eigen::VectorXd& b; /* b of size M*/
	const Eigen::VectorXd& b_bar; /* b_bar of size M*/
	const Eigen::VectorXd& c; /* c of size N*/
	const Eigen::VectorXd& c_bar; /* c_bar of size N*/

	std::vector<int> B;    /*basic indices size of M*/
	std::vector<int> NB;    /*non-basic indices size of N-M*/
	std::vector<int> inner_dict;

	Eigen::VectorXd E_d;
	Eigen::MatrixXd Eta;
	Eigen::MatrixXd A_N_t;

	/* LU solve state for dxb */
	int dxb_itern_;
	bool dxb_first_;
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> dxb_solver_;

	/* LU solve state for dzn */
	int dzn_itern_;
	bool dzn_first_;
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> dzn_solver_;

	PSM(const Eigen::MatrixXd& _A,
		const Eigen::VectorXd& _b,
		const Eigen::VectorXd& _b_bar,
		const Eigen::VectorXd& _c,
		const Eigen::VectorXd& _c_bar);
	void init();
	void reset_lu_state();
	Eigen::VectorXd lusolve_update_dxb(int col_in);
	Eigen::VectorXd lusolve_update_dzn(int col_out);
	PSMresult solve(int max_it, double lambda_threshold);
};

#endif
