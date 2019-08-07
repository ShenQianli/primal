#ifndef PSM_H
#define PSM_H

#include <cstring>  /* malloc calloc memset*/
#include <cfloat>	/* DBL_MAX DBL_MIN */
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define EPS 1.0e-100

using namespace std;
using namespace Eigen;

enum FLAGTYPE{NONE, BASIC, NONBASIC};

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
	void init();
	VectorXd lusolve_update_dxb(int col_in);
	VectorXd lusolve_update_dzn(int col_out);
	PSMresult solve(int max_it, double lambda_threshold);
};

#endif
