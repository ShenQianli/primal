#ifndef API_HPP
#define API_HPP

#include <Eigen/Dense>

#include "PSM/PLP.hpp"
#include "PSM/PSM.hpp"
#include "PSM/problems.hpp"

extern "C" void QuantileRegression_api
(	int *n,	/*row num*/
	int *d,	/*col num*/
	double *X,
	double *y,
 	int max_it,
 	double lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

extern "C" void SparseSVM_api
(	int *n,	/*row num*/
	int *d,	/*col num*/
	double *X,
	double *y,
	int max_it,
	double lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list,
 	double *x0
);

extern "C" void Dantzig_api
(	int *n,	/*row num*/
	int *d,	/*col num*/
	double *X,
	double *y,
	int max_it,
	double lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

extern "C" void CompressedSensing_api
(	int *n,	/*row num*/
	int *d,	/*col num*/
	double *X,
	double *y,
	int max_it,
	double lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

#endif
