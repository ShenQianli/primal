#ifndef PSM_H
#define PSM_H

#include <cstring>  /* malloc calloc memset*/
#include <cfloat>	/* DBL_MAX DBL_MIN */
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "PSM/PLP.hpp"

#define EPS 1.0e-100

using namespace std;
using namespace Eigen;

enum FLAGTYPE{NONE, BASIC, NONBASIC};

class PSMresult{
public:
	int T;
	int d;
	double *lambda_list;
	double *x_list;
	double *y_list;
	PSMresult(int max_it, int _d);
	~PSMresult();
	void update(double lambda, VectorXd x, double y);
};

class PSM{
public:
	int m;
	int n;
	int M; /*M = m*/
	int N; /*N = n+m*/
	MatrixXd A; /* A of size M*N */
	VectorXd b; /* b of size M*/
	VectorXd b_bar; /* b_bar of size M*/
	VectorXd c; /* c of size N*/
	VectorXd c_bar; /* c_bar of size N*/
	
	int *B;    /*basic indices size of m*/
	int *NB;    /*non-basic indices size of n-m*/
	int *inner_dict;
	
	VectorXd E_d;
	MatrixXd Eta;
	MatrixXd A_N_t;
	
	PSM(PLP *pPLP);
	~PSM();
	void init();
	VectorXd lusolve_update_dxb(int col_in);
	VectorXd lusolve_update_dzn(int col_out);
	PSMresult solve(int max_it, double lambda_threshold);
};

#endif
