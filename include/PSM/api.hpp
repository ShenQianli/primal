#ifndef API_HPP
#define API_HPP

#include <Eigen/Dense>

#include "PSM/PLP.hpp"
#include "PSM/problems.hpp"

PLP *LinearProgrammingWrapper(
	int *_m,	/*row num*/
	int *_n,	/*col num*/
	double *_A,	/* A of size m*n*/
	double *_b,	/* b of size m*/
	double *_b_bar,	/* b_bar of size m*/
	double *_c,	/* c of size n*/
	double *_c_bar	/* c_bar of size n*/
);

PLP *QuantileRegressionWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
);

PLP *SparseSVMWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
);

PLP *DantzigWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
);

PLP *CompressedSensingWrapper(
	int *_n,	/*row num*/
	int *_d,	/*col num*/
	double *_X,
	double *_y
);

#endif
