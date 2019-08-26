#include "PSM/PSM.hpp"

PSMresult::PSMresult(int max_it, int _d){
	T = 0;
	d = _d;
	lambda_list = (double *)malloc(max_it * sizeof(double));
	x_list.resize(d, max_it);
	y_list = (double *)malloc(max_it * sizeof(double));
}

PSMresult::~PSMresult(){
	free(lambda_list);
	free(y_list);
}

void PSMresult::update(double lambda, VectorXd x, double y){
	lambda_list[T] = lambda;
	x_list.col(T) = x;
	y_list[T] = y;
	T++;
}

PSM::PSM(const MatrixXd& _A,
		 const VectorXd& _b,
		 const VectorXd& _b_bar,
		 const VectorXd& _c,
		 const VectorXd& _c_bar
		 ):A(_A),b(_b),b_bar(_b_bar),c(_c),c_bar(_c_bar){
	
	M = A.rows();
	N = A.cols();
	m = M;
	n = N - m;
	inner_dict = (int *)malloc(N*sizeof(int));
	B = (int *)malloc(m*sizeof(int));
	NB = (int *)malloc(n*sizeof(int));
	E_d.resize(m);
	Eta.resize(m,m);
	A_N_t.resize(n,m);
}

PSM::~PSM(){
	free(B);
	free(NB);
	free(inner_dict);
}

void PSM::init(int *B_init){
	memset(B, 0, m*sizeof(int));
	memset(NB, 0, n*sizeof(int));
	if(B_init == NULL){
		for(int i = 0; i < m; ++i){
			B[i] = i + n;
			inner_dict[i + n] = i;
		}
		for(int i = 0; i < n; ++i){
			NB[i] = i;
			inner_dict[i] = i;
		}
	}
	else{
		int p = 0, q = 0;
		for(int i = 0; i < N; ++i){
			if(i == B_init[p]){
				B[p] = i;
				inner_dict[i] = p++;
			}
			else{
				NB[q] = i;
				inner_dict[i] = q++;
			}
		}
	}
	E_d.setZero();
	Eta.setZero();
	for (int i=0;i<n;i++)
	{
		A_N_t.row(i)=A.col(NB[i]);
	}
}

VectorXd PSM::A_B_solve(VectorXd y, int itern){
	static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
	MatrixXd A_B(m,m);
	VectorXd x(m);
	SparseMatrix<double> A_B_sparse(m,m);
	int i;
	if(itern>=m/2||itern==0){
		for (i=0;i<m;i++)
		{
			A_B.col(i)=A.col(B[i]);
		}
		A_B_sparse = A_B.sparseView();
		solver.analyzePattern(A_B_sparse);
		solver.factorize(A_B_sparse);
	}
	x = solver.solve(y);
	return x;
}

VectorXd PSM::A_B_t_solve(VectorXd y, int itern){
	static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
	MatrixXd A_B(m,m);
	VectorXd x(m);
	SparseMatrix<double> A_B_sparse(m,m);
	SparseMatrix<double> A_B_sparse_t(m,m);
	int i;
	if(itern>=m/2||itern==0){
		for (i=0;i<m;i++)
		{
			A_B.col(i)=A.col(B[i]);
		}
		A_B_sparse = A_B.sparseView();
		A_B_sparse_t = A_B_sparse.transpose();
		solver.analyzePattern(A_B_sparse_t);
		solver.factorize(A_B_sparse_t);
	}
	x = solver.solve(y);
	return x;
}


VectorXd PSM::lusolve_update_dxb(int col_in){
	/*input i return A_B^(-1)*A_N*ei*/
	static int itern = 0;
	if(col_in == -1){
		itern = 0;
		return MatrixXd::Zero(1,1);
	}
	int i,j=0;
	int col;
	double temp;
	VectorXd x(m);
	VectorXd y(m);
	y = A.col(col_in);
	if(itern>=m/2||itern==0){
		itern=0;
	}
	x = A_B_solve(y, itern);
	
	if (itern <= 0) {
		Eta.setZero();
		Eta.col(itern)=x;
		itern++;
		return x;
	}
	for (j=0; j<itern; j++) {
		col = E_d(j);
		
		temp = x(col)/Eta(col,j);
		if (temp != 0.0) {
			for (i=0;i<m;i++) {
				x(i) -= Eta(i,j) * temp;
			}
			x(col) = temp;
		}
	}
	//save the updated data
	Eta.col(itern)=x;
	itern++;
	return x;
}

VectorXd PSM::lusolve_update_dzn(int col_out){
	/*input j returnA_N_t*A_B^(-1)^T*ej*/
	static int itern = 0;
	if(col_out == -1){
		itern = 0;
		return MatrixXd::Zero(1,1);
	}
	int j = 0,k;
	int col;
	double temp;
	VectorXd x(n);
	SparseVector<double> y(m);
	static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
	y.insert(inner_dict[col_out])=-1;
	x.setZero();
	if(itern>=m/2||itern==0){
		itern=0;
	}
	if (itern == 0){
		x = A_B_t_solve(y,itern);
		itern++;
		E_d.setZero();
		E_d(0) = inner_dict[col_out];
		x = A_N_t*x;
		return x;
	}
	for (j=itern-1; j>=0; j--) {
		temp=0;
		col = E_d(j);
		for(SparseVector<double>::InnerIterator it(y); it; ++it){
			if(it.index()>=0){
				if(it.index()!=col){
					k=it.index();
					temp -= Eta(k,j)*it.value();
				} else temp += it.value();
			}
		}
		y.coeffRef(col) = temp/Eta(col,j);
	}
	x = A_B_t_solve(y, itern);
	E_d(itern)=inner_dict[col_out];
	itern++;
	x = A_N_t*x.sparseView();
	return x;
}

PSMresult PSM::solve(int max_it, double lambda_threshold, int *B_init){
	/*local variables*/
	FLAGTYPE flag = NONE;
	int T = 0;
	int col_in = -1;
	int col_out = -1;
	VectorXd xB_star(m);
	VectorXd xB_bar(m);
	VectorXd dxB(m);
	VectorXd zN_star(n);
	VectorXd zN_bar(n);
	VectorXd dzN(n);
	VectorXd xB(m);
	VectorXd x_output(N);
	double t;
	double t_bar;
	double s;
	double s_bar;
	double lambda_star;
	double y_output;
	PSMresult result(max_it, N);
	
	/*initialize*/
	init(B_init);
//	xB_star = b;
//	xB_bar = b_bar;
//	zN_star = A_N_t * c.tail(m) - c.head(n);
//	zN_bar = A_N_t * c_bar.tail(m) - c_bar.head(n);
	xB_star = A_B_solve(b, 0);
	xB_bar = A_B_solve(b_bar, 0);
	VectorXd cB(m);
	VectorXd cN(n);
	VectorXd cB_bar(m);
	VectorXd cN_bar(n);
	for(int i = 0; i < m; ++i){
		cB(i) = c(B[i]);
		cB_bar(i) = c_bar(B[i]);
	}
	for(int i = 0; i < n; ++i){
		cN(i) = c(NB[i]);
		cN_bar(i) = c_bar(NB[i]);
	}
	zN_star = A_N_t * A_B_t_solve(cB, 0) - cN;
	zN_bar = A_N_t * A_B_t_solve(cB_bar, 0) - cN_bar;
	
	/*set initial Basic*/
	if(zN_bar.isZero()){
		flag = PRIMAL;
	}
	else if(xB_bar.isZero()){
		flag = DUAL;
	}
	else{
		return result;
	}
	
	/*Simplex Method*/
	while(1){
		/*break when primal(dual) optimal if flag is PRIMAL(DUAL)*/
		col_in = -1;
		col_out = -1;
		if(flag == PRIMAL){
			double min = -EPS;
			for (int i = 0; i < n; ++i)
			{
				if(zN_star(i) < min){
					min = zN_star(i);
					col_in = NB[i];
				}
			}
			if(col_in == -1){
				/*optimal*/
				break;
			}
			dxB = lusolve_update_dxb(col_in);
			
			min = DBL_MAX;
			for (int i = 0; i < m; ++i)
			{
				if(xB_bar(i) == 0 && dxB(i) > EPS){
					double tmp = xB_star(i) / dxB(i);
					if(tmp < min){
						min = tmp;
						col_out = B[i];
					}
				}
			}
			if(col_out == -1){
				break;
			}
			dzN = lusolve_update_dzn(col_out);
		}
		else if(flag == DUAL){
			double min = -EPS;
			for (int i = 0; i < m; ++i)
			{
				if(xB_star(i) < min){
					min = xB_star(i);
					col_out = B[i];
				}
			}
			if(col_out == -1){
				/*optimal*/
				break;
			}
			dzN = lusolve_update_dzn(col_out);
			min = DBL_MAX;
			for (int i = 0; i < n; ++i)
			{
				if(zN_bar(i) == 0 && dzN(i) > EPS){
					double tmp = zN_star(i) / dzN(i);
					if(tmp < min){
						min = tmp;
						col_in = NB[i];
					}
				}
			}
			if(col_in == -1){
				break;
			}
			dxB = lusolve_update_dxb(col_in);
		}
		
		/*Compute the dual and primal step lengths \
		 for both variables and perturbations*/
		t = xB_star[inner_dict[col_out]] / dxB[inner_dict[col_out]];
		t_bar = xB_bar[inner_dict[col_out]] / dxB[inner_dict[col_out]];
		s = zN_star[inner_dict[col_in]] / dzN[inner_dict[col_in]];
		s_bar = zN_bar[inner_dict[col_in]] / dzN[inner_dict[col_in]];
		
		/*Update the primal and dual solutions*/
		xB_star = xB_star - t * dxB;
		xB_bar = xB_bar - t_bar * dxB;
		zN_star = zN_star - s * dzN;
		zN_bar = zN_bar - s_bar * dzN;
		
		zN_star(inner_dict[col_in]) = s;
		zN_bar(inner_dict[col_in]) = s_bar;
		xB_star(inner_dict[col_out]) = t;
		xB_bar(inner_dict[col_out]) = t_bar;
		
		/*Update the basic and nonbasic index sets*/
		A_N_t.row(inner_dict[col_in])=A.col(col_out);
		B[inner_dict[col_out]] = col_in;
		NB[inner_dict[col_in]] = col_out;
		swap(inner_dict[col_in], inner_dict[col_out]);
	}
	while(T < max_it){
		T++;
		lambda_star = 0;
		col_in = -1;
		col_out = -1;
		if(flag == DUAL){
			/*compute lambda_star and select entering variable*/
			for (int i = 0; i < n; ++i)
			{
				if(zN_bar(i) > EPS){
					double tmp = - zN_star(i) / zN_bar(i);
					if(tmp > lambda_star){
						lambda_star = tmp;
						col_in = NB[i];
					}
				}
			}
			if(col_in == -1){
				break;
			}
			/*Compute primal step direction*/
			dxB = lusolve_update_dxb(col_in);
			
			/*Select leaving variable*/
			double min = DBL_MAX;
			for (int i = 0; i < m; ++i)
			{
				if(dxB(i) > EPS){
					double tmp;
					tmp = (xB_star(i) + lambda_star * xB_bar(i)) / dxB(i);
					if(tmp < min){
						min = tmp;
						col_out = B[i];
					}
				}
			}
			if(col_out == -1){
				break;
			}
			/*Compute dual step direction*/
			dzN = lusolve_update_dzn(col_out);
		}
		else if(flag == PRIMAL){
			/*compute lambda_star and select leaving variable*/
			for (int i = 0; i < m; ++i)
			{
				if(xB_bar(i) > EPS){
					double tmp = - xB_star(i) / xB_bar(i);
					if(tmp > lambda_star){
						lambda_star = tmp;
						col_out = B[i];
					}
				}
			}
			if(col_out == -1){
				break;
			}
			/*Compute dual step direction*/
			dzN = lusolve_update_dzn(col_out);
			
			/*Select entering variable*/
			double min = DBL_MAX;
			for (int i = 0; i < n; ++i)
			{
				if(dzN(i) > EPS){
					double tmp;
					tmp = (zN_star(i) + lambda_star * zN_bar(i)) / dzN(i);
					if(tmp < min){
						min = tmp;
						col_in = NB[i];
					}
				}
			}
			if(col_in == -1){
				break;
			}
			/*Compute primal step direction*/
			dxB = lusolve_update_dxb(col_in);
		}
		
		/*Compute the dual and primal step lengths \
		 for both variables and perturbations*/
		t = xB_star[inner_dict[col_out]] / dxB[inner_dict[col_out]];
		t_bar = xB_bar[inner_dict[col_out]] / dxB[inner_dict[col_out]];
		s = zN_star[inner_dict[col_in]] / dzN[inner_dict[col_in]];
		s_bar = zN_bar[inner_dict[col_in]] / dzN[inner_dict[col_in]];
		
		/*Update the primal and dual solutions*/
		xB_star = xB_star - t * dxB;
		xB_bar = xB_bar - t_bar * dxB;
		zN_star = zN_star - s * dzN;
		zN_bar = zN_bar - s_bar * dzN;
		
		zN_star(inner_dict[col_in]) = s;
		zN_bar(inner_dict[col_in]) = s_bar;
		xB_star(inner_dict[col_out]) = t;
		xB_bar(inner_dict[col_out]) = t_bar;
		
		/*Update the basic and nonbasic index sets*/
		A_N_t.row(inner_dict[col_in])=A.col(col_out);
		B[inner_dict[col_out]] = col_in;
		NB[inner_dict[col_in]] = col_out;
		swap(inner_dict[col_in], inner_dict[col_out]);
		
		/*Update the output result*/
		xB = xB_star + lambda_star * xB_bar;
		x_output.setZero();
		for(int i = 0; i < m; ++i){
			x_output[B[i]] = xB[i];
		}
		y_output = x_output.transpose()*(c + lambda_star * c_bar);
		result.update(lambda_star , x_output, y_output);
		/*check threshold*/
		if(lambda_star <= lambda_threshold){
			break;
		}
	}
	
	/*reset*/
	lusolve_update_dxb(-1);
	lusolve_update_dzn(-1);
	
	return result;
}
