#include "PSM/PSM.hpp"

PSMresult::PSMresult(int max_it, int _d){
    T = 0;
    d = _d;
    lambda_list.resize(max_it, 0.0);
    x_list.resize(d, max_it);
    y_list.resize(max_it, 0.0);
}

void PSMresult::update(double lambda, Eigen::VectorXd x, double y){
    lambda_list[T] = lambda;
    x_list.col(T) = x;
    y_list[T] = y;
    T++;
}

PSM::PSM(const Eigen::MatrixXd& _A,
		 const Eigen::VectorXd& _b,
		 const Eigen::VectorXd& _b_bar,
		 const Eigen::VectorXd& _c,
		 const Eigen::VectorXd& _c_bar
		 ):A(_A),b(_b),b_bar(_b_bar),c(_c),c_bar(_c_bar){

    M = A.rows();
	N = A.cols();
	m = M;
	n = N - m;
    inner_dict.resize(N, 0);
    B.resize(m, 0);
    NB.resize(n, 0);
    E_d.resize(m);
    Eta.resize(m,m);
    A_N_t.resize(n,m);

    dxb_itern_ = 0;
    dxb_first_ = true;
    dzn_itern_ = 0;
    dzn_first_ = true;
}

void PSM::reset_lu_state(){
    dxb_itern_ = 0;
    dxb_first_ = true;
    dzn_itern_ = 0;
    dzn_first_ = true;
}

void PSM::init(){
    std::memset(B.data(), 0, m*sizeof(int));
    std::memset(NB.data(), 0, n*sizeof(int));
    for(int i = 0; i < m; ++i){
        B[i] = i + n;
        inner_dict[i + n] = i;
    }
    for(int i = 0; i < n; ++i){
        NB[i] = i;
        inner_dict[i] = i;
    }
    E_d.setZero();
    Eta.setZero();
    for (int i = 0; i < n; i++)
    {
        A_N_t.row(i) = A.col(NB[i]);
    }
}

Eigen::VectorXd PSM::lusolve_update_dxb(int col_in){
    /*input i return A_B^(-1)*A_N*ei*/
    int i = 0, j = 0;
    int col = 0;
    double temp = 0.0;
    Eigen::VectorXd x(m);
    Eigen::VectorXd y(m);
    Eigen::MatrixXd A_B(m,m);
    Eigen::SparseMatrix<double> A_B_sparse(m,m);
    y = A.col(col_in);
    if(dxb_itern_ >= m/2) dxb_first_ = false;
    if((dxb_itern_ >= m/2 || dxb_itern_ == 0) && dxb_first_ == false){
        for (i = 0; i < m; i++)
        {
            A_B.col(i) = A.col(B[i]);
        }
        A_B_sparse = A_B.sparseView();
        dxb_solver_.analyzePattern(A_B_sparse);
        dxb_solver_.factorize(A_B_sparse);
        dxb_itern_ = 0;
    }
    if(dxb_first_ == true){
        x = y;
    }else x = dxb_solver_.solve(y);

    if (dxb_itern_ <= 0) {
        Eta.setZero();
        Eta.col(dxb_itern_) = x;
        dxb_itern_++;
        return x;
    }
    for (j = 0; j < dxb_itern_; j++) {
        col = E_d(j);

        temp = x(col)/Eta(col,j);
        if (temp != 0.0) {
            for (i = 0; i < m; i++) {
                x(i) -= Eta(i,j) * temp;
            }
            x(col) = temp;
        }
    }
    //save the updated data
    Eta.col(dxb_itern_) = x;
    dxb_itern_++;
    return x;
}

Eigen::VectorXd PSM::lusolve_update_dzn(int col_out){
    /*input j return A_N_t*A_B^(-1)^T*ej*/
    int i = 0, j = 0, k = 0;
    int col = 0;
    double temp = 0.0;
    Eigen::VectorXd x(n);
    Eigen::MatrixXd A_B(m,m);
    Eigen::SparseMatrix<double> A_B_sparse(m,m);
    Eigen::SparseMatrix<double> A_B_sparse_t(m,m);
    Eigen::SparseVector<double> y(m);
    y.insert(inner_dict[col_out]) = -1;
    x.setZero();
    if(dzn_itern_ >= m/2) dzn_first_ = false;
    if((dzn_itern_ >= m/2 || dzn_itern_ == 0) && dzn_first_ == false){
        for (i = 0; i < m; i++)
        {
            A_B.col(i) = A.col(B[i]);
        }
        A_B_sparse = A_B.sparseView();
        A_B_sparse_t = A_B_sparse.transpose();
        dzn_solver_.analyzePattern(A_B_sparse_t);
        dzn_solver_.factorize(A_B_sparse_t);
        dzn_itern_ = 0;
    }
    if (dzn_itern_ == 0){
        if(dzn_first_ == true){
            x = y;
        }else x = dzn_solver_.solve(y);
        dzn_itern_++;
        E_d.setZero();
        E_d(0) = inner_dict[col_out];
        x = A_N_t * x;
        return x;
    }
    for (j = dzn_itern_ - 1; j >= 0; j--) {
        temp = 0;
        col = E_d(j);
        for(Eigen::SparseVector<double>::InnerIterator it(y); it; ++it){
            if(it.index() >= 0){
                if(it.index() != col){
                    k = it.index();
                    temp -= Eta(k,j) * it.value();
                } else temp += it.value();
            }
        }
        y.coeffRef(col) = temp / Eta(col,j);
    }
    if(dzn_first_ == true){
        x = y;
    }else x = dzn_solver_.solve(y);
    E_d(dzn_itern_) = inner_dict[col_out];
    dzn_itern_++;
    x = A_N_t * x.sparseView();
    return x;
}

PSMresult PSM::solve(int max_it, double lambda_threshold){
    /*local variables*/
    FLAGTYPE flag = NONE;
    int T = 0;
    int col_in = -1;
    int col_out = -1;
    Eigen::VectorXd xB_star(m);
    Eigen::VectorXd xB_bar(m);
    Eigen::VectorXd dxB(m);
    Eigen::VectorXd zN_star(n);
    Eigen::VectorXd zN_bar(n);
    Eigen::VectorXd dzN(n);
	Eigen::VectorXd xB(m);
    Eigen::VectorXd x_output(N);
    double t = 0.0;
    double t_bar = 0.0;
    double s = 0.0;
    double s_bar = 0.0;
    double lambda_star = 0.0;
    double y_output = 0.0;
    PSMresult result(max_it, N);

    /*initialize*/
    init();
    reset_lu_state();
	xB_star = b;
	xB_bar = b_bar;
	zN_star = A_N_t * c.tail(m) - c.head(n);
	zN_bar = A_N_t * c_bar.tail(m) - c_bar.head(n);

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
		A_N_t.row(inner_dict[col_in]) = A.col(col_out);
		B[inner_dict[col_out]] = col_in;
		NB[inner_dict[col_in]] = col_out;
		std::swap(inner_dict[col_in], inner_dict[col_out]);
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
			if(col_in != -1){
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
				if(col_out != -1){
					/*Compute dual step direction*/
					dzN = lusolve_update_dzn(col_out);
				}
			}
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
			if(col_out != -1){
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
				if(col_in != -1){
					/*Compute primal step direction*/
					dxB = lusolve_update_dxb(col_in);
				}
			}
		}

		/*Update the output result*/
		xB = xB_star + lambda_star * xB_bar;
		x_output.setZero();
		for(int i = 0; i < m; ++i){
			x_output[B[i]] = xB[i];
		}
		y_output = x_output.transpose()*(c + lambda_star * c_bar);
		result.update(lambda_star , x_output, y_output);

		/*check threshold*/
		if(lambda_star <= lambda_threshold || col_in == -1 || col_out == -1){
			break;
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
		A_N_t.row(inner_dict[col_in]) = A.col(col_out);
		B[inner_dict[col_out]] = col_in;
		NB[inner_dict[col_in]] = col_out;
		std::swap(inner_dict[col_in], inner_dict[col_out]);
	}

    /*reset*/
    reset_lu_state();

    return result;
}
