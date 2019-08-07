#include <iostream>
#include <Eigen/Dense>

#include "PSM/api.hpp"
#include <iostream>

using namespace std;

int main()
{
	const int n = 20;
	const int d = 100;
	int nn = 20;
	int dd = 100;
	MatrixXd X(n, d);
	VectorXd y(n);
	X.setRandom();
	y.setRandom();
	double _X[n*d];
	double _y[n];
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < d; ++j){
			_X[d*i+j] = X(i, j);
		}
		_y[i] = y(i);
	}
	
	int max_it = 100;
	double lambda_threshold = 0.001;
	int T;
	double lambda_list[100];
	double x_list[100*d];
	double y_list[100];
	double x0_list[100];
	
	Dantzig_api(&nn, &dd, _X, _y, &max_it, &lambda_threshold, &T, lambda_list, x_list, y_list);

	for(int i = 0; i < T; ++i){
		cout << i << " " << lambda_list[i] << " " << y_list[i]<<endl;
	}
}
