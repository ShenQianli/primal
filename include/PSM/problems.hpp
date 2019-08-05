#ifndef PROBLEMS_HPP
#define PROBLEMS_HPP

#include <Eigen/Dense>

#include "PSM/PLP.hpp"

using namespace Eigen;

PLP *QuantileRegression(const MatrixXd& X, const VectorXd& y);

PLP *SparseSVM(const MatrixXd& X, const VectorXd& y);

PLP *Dantzig(const MatrixXd& X, const VectorXd& y);

PLP *CompressedSensing(const MatrixXd& X, const VectorXd& y);

#endif
