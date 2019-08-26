#ifndef PSM_H
#define PSM_H

#include <iostream>
#include <cstring>  /* malloc calloc memset*/
#include <cfloat>	/* DBL_MAX DBL_MIN */
#include <eigen/Eigen/Sparse>
#include <eigen/Eigen/Dense>

#define EPS 1e-10

using namespace std;
using namespace Eigen;

enum FLAGTYPE{NONE, PRIMAL, DUAL};

class PSMresult{
public:
	int T;
	int d;
	double *lambda_list;
	MatrixXd x_list;
	double *y_list;
	PSMresult(int max_it, int _d);
	~PSMresult();
	void update(double lambda, VectorXd x, double y);
};

class PSM{
public:
	/*Parametric Linear Programming problem parameters*/
	int M; /*row num*/
	int N; /*col num*/
	int m; /*m = M*/
	int n; /*n = N-M*/
	const MatrixXd& A; /* A of size M*N */
	const VectorXd& b; /* b of size M*/
	const VectorXd& b_bar; /* b_bar of size M*/
	const VectorXd& c; /* c of size N*/
	const VectorXd& c_bar; /* c_bar of size N*/
	
	int *B;    /*basic indices size of M*/
	int *NB;    /*non-basic indices size of N-M*/
	int *inner_dict;
	
	VectorXd E_d;
	MatrixXd Eta;
	MatrixXd A_N_t;
	
	PSM(const MatrixXd& _A,
		const VectorXd& _b,
		const VectorXd& _b_bar,
		const VectorXd& _c,
		const VectorXd& _c_bar);
	~PSM();
	void init(int *B_init = NULL);
	VectorXd A_B_solve(VectorXd y, int itern);
	VectorXd A_B_t_solve(VectorXd y, int itern);
	VectorXd lusolve_update_dxb(int col_in);
	VectorXd lusolve_update_dzn(int col_out);
	PSMresult solve(int max_it, double lambda_threshold, int *B_init = NULL);
};

#endif
