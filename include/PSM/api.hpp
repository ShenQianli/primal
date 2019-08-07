#ifndef API_HPP
#define API_HPP

#include <Eigen/Dense>

#include "PSM/PSM.hpp"


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
);

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
);

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
);

extern "C" void CompressedSensing_api
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
);

#endif
