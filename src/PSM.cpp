#include "PSM/PSM.hpp"

bool iszero(double x){
	return x > -EPS && x < EPS;
}

PSMresult::PSMresult(int max_it, int _d){
	T = 0;
	d = _d;
	lambda_list = (double *)calloc(max_it, sizeof(double));
	x_list = (double *)calloc(max_it * d, sizeof(double));
	y_list = (double *)calloc(max_it, sizeof(double));
}

PSMresult::~PSMresult(){
	free(lambda_list);
	free(x_list);
	free(y_list);
}

void PSMresult::update(double lambda, VectorXd x, double y){
	lambda_list[T] = lambda;
	for(int i = 0; i < d; ++i)
		x_list[T*d+i] = x(i);
	y_list[T] = y;
	T++;
}

PSM::PSM(PLP *pPLP){
	m = pPLP->M;
	n = pPLP->N - m;
	M = m;
	N = pPLP->N;
	A = pPLP->A;
	b = pPLP->b;
	b_bar = pPLP->b_bar;
	c = pPLP->c;
	c_bar = pPLP->c_bar;
	inner_dict = (int *)malloc(N*sizeof(int));
	B = (int *)malloc(M*sizeof(int));
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

void PSM::init(){
	memset(B, 0, M*sizeof(int));
	memset(NB, 0, n*sizeof(int));
	for(int i = 0; i < M; ++i){
		B[i] = i + n;
		inner_dict[i + n] = i;
	}
	for(int i = 0; i < M; ++i){
		B[i] = i + n;
		inner_dict[i + n] = i;
	}
	E_d.setZero();
	Eta.setZero();
	for (int i=0;i<n;i++)
	{
		A_N_t.row(i)=A.col(NB[i]);
	}
}

VectorXd PSM::lusolve_update_dxb(int col_in){
	/*input i return A_B^(-1)*A_N*ei*/
	static int itern = 0;
	static bool isfirsttime = true;
	if(col_in == -1){
		itern = 0;
		isfirsttime = true;
		return MatrixXd::Zero(1,1);
	}
	int i,j=0;
	int col;
	double temp;
	VectorXd x(m);
	VectorXd y(m);
	MatrixXd A_B(m,m);
	SparseMatrix<double> A_B_sparse(m,m);
	static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
	y = A.col(col_in);
	if(itern>=m/2) isfirsttime=false;
	if((itern>=m/2||itern==0)&&isfirsttime==false){
		for (i=0;i<m;i++)
		{
			A_B.col(i)=A.col(B[i]);
		}
		A_B_sparse = A_B.sparseView();
		solver.analyzePattern(A_B_sparse);
		solver.factorize(A_B_sparse);
		itern=0;
		
	}
	if(isfirsttime==true){
		x = y;
	}else x = solver.solve(y);
	
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
	static bool isfirsttime = true;
	if(col_out == -1){
		itern = 0;
		isfirsttime = true;
		return MatrixXd::Zero(1,1);
	}
	int i,j = 0,k;
	int col;
	double temp;
	VectorXd x(n);
	MatrixXd A_B(m,m);
	SparseMatrix<double> A_B_sparse(m,m);
	SparseMatrix<double> A_B_sparse_t(m,m);
	SparseVector<double> y(m);
	static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
	y.insert(inner_dict[col_out])=-1;
	x.setZero();
	if(itern>=m/2) isfirsttime=false;
	if((itern>=m/2||itern==0)&&isfirsttime==false){
		for (i=0;i<m;i++)
		{
			A_B.col(i)=A.col(B[i]);
		}
		A_B_sparse = A_B.sparseView();
		A_B_sparse_t = A_B_sparse.transpose();
		solver.analyzePattern(A_B_sparse_t);
		solver.factorize(A_B_sparse_t);
		itern=0;
	}
	if (itern == 0){
		if(isfirsttime==true){
			x = y;
		}else x = solver.solve(y);
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
	if(isfirsttime==true){
		x = y;
	}else x = solver.solve(y);
	E_d(itern)=inner_dict[col_out];
	itern++;
	x = A_N_t*x.sparseView();
	return x;
	}

PSMresult PSM::solve(int max_it, double lambda_threshold){
	/*local variables*/
	FLAGTYPE flag = NONE;
	int T = 0;
	int col_in = -1;
	int col_out = -1;
	VectorXd xB_star(M);
	VectorXd xB_bar(M);
	VectorXd dxB(M);
	VectorXd zN_star(n);
	VectorXd zN_bar(n);
	VectorXd dzN(n);
	double t;
	double t_bar;
	double s;
	double s_bar;
	double lambda_star;
	PSMresult result(max_it, m);
	
	/*initialize*/
	init();
	xB_star = b;
	xB_bar = b_bar;
	zN_star = -c.head(n);
	zN_bar = -c_bar.head(n);
	
	while(T < max_it){
		T++;
		
		/*compute lambda_star*/
		lambda_star = DBL_MIN;
		col_in = -1;
		col_out = -1;
		flag = NONE;
		for (int i = 0; i < n; ++i)
		{
			if(zN_bar(i) > EPS){
				double tmp = - zN_star(i) / zN_bar(i);
				if(tmp > lambda_star){
					lambda_star = tmp;
					col_in = NB[i];
					flag = NONBASIC;
				}
			}
		}
		for (int i = 0; i < M; ++i)
		{
			if(xB_bar(i) > EPS){
				double tmp = - xB_star(i) / xB_bar(i);
				if(tmp > lambda_star){
					lambda_star = tmp;
					col_out = B[i];
					flag = BASIC;
				}
			}
		}
		
		/*check threshold*/
		if(lambda_star <= lambda_threshold){
			break;
		}
		switch (flag) {
			case NONBASIC:
			{
				/*Compute primal step direction*/
				dxB = lusolve_update_dxb(col_in);
				
				/*Select leaving variable*/
				double max = DBL_MIN;
				for (int i = 0; i < M; ++i)
				{
					if(!iszero(xB_star[i] + lambda_star * xB_bar[i])){
						double tmp = dxB[i] / \
						(xB_star[i] + lambda_star * xB_bar[i]);
						if(tmp > max){
							max = tmp;
							col_out = B[i];
						}
					}
				}
				
				/*Compute dual step direction*/
				dzN = lusolve_update_dzn(col_out);
				
				break;
			}
			case BASIC:
			{
				/*Compute dual step direction*/
				dzN = lusolve_update_dzn(col_out);
				
				/*Select leaving variable*/
				double max = DBL_MIN;
				for (int i = 0; i < n; ++i)
				{
					if(!iszero(zN_star[i] + lambda_star * zN_bar[i])){
						double tmp = dzN[i] / \
						(zN_star[i] + lambda_star * zN_bar[i]);
						if(tmp > max){
							max = tmp;
							col_in = NB[i];
						}
					}
				}
				
				/*Compute primal step direction*/
				dxB = lusolve_update_dxb(col_in);
				break;
			}
			default:
			{
				break;
			}
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
	/*reset*/
	lusolve_update_dxb(-1);
	lusolve_update_dzn(-1);
	
	return result;
}
