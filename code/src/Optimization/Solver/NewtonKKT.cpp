#include "NewtonKKT.h"

#include "Quad.h"
#include "line_search.h"
#include "Logger.h"

#include "igl/Timer.h"

#include "igl/cat.h"
#include "igl/matlab_format.h"

using namespace std;
//m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
                //m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
                //m_solver.iparm[20] = 1;
double NewtonKKT::solve_constrained(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x,
        double convergence_threshold) {

    if (lambda.rows()!= constraints.getConstNum()) {
        lambda.resize(constraints.getConstNum());
        lambda.setZero();
        //int wait; std::cout << "init lagrange mult" << std::endl; std::cin >> wait;
    }

    auto prev_x = x0; double current_merit_p = merit_p; x = x0;
    current_merit_p = max(0.05,lambda.cwiseAbs().maxCoeff()*1.1); // (18.32) in Nocedal
    double ret = f.obj(x0); int iter = 0;
	write_log(4) << "convergence_threshold = " << convergence_threshold << endl;
    //if (const_dev > infeasability_filter) { x = prev_x; current_merit_p *=2;} prev_x = x;
    //while ((const_dev > infeasability_epsilon) && iter < max_newton_iters) {
    double error = 10000*kkt_error(x, f, constraints);
    //error = convergence_threshold+1e-2; // hack to get at least one iteration. Problem is that the gradient condition is too tight, maybe.
    // this might be due to the modified jacobian. We should probably use the normal non modified one for the termination condition..
    while ( (error > infeasability_epsilon) && (iter < max_newton_iters) ) {
        double ret = one_iter(x,f,constraints,x,current_merit_p);
        double const_dev = constraints.Vals(x).norm();
        iter++;
		write_log(4) << "KKT error = " << error << ", const_dev = " << const_dev << " after " << iter  << " iters" << std::endl;
        error = kkt_error(x, f, constraints);
        //if (const_dev > infeasability_filter) { x = prev_x; current_merit_p *=2; iter--;/*do another iter*/} prev_x = x;
        //current_merit_p *=2;
    }
	write_log(3) << "finished optimization after " << iter << " / " << max_newton_iters << " iterations" << endl;
	return ret;
}

void NewtonKKT::build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, int var_n,
                                     const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n) {
    // build an IJV which contains the hessian, diag matrix along the hessian, jacobian, jacobian transpose,
    //  and diagonal matrix along the J to make sure the whole diagonal is nnz (needed for Pardiso solver)
    int ijv_idx = 0;
    if (kkt_IJV.size() == 0) {
        int kkt_size = hessian_IJV.size() + 2*jacobian_IJV.size() + var_n + const_n;
        kkt_IJV.resize(kkt_size); int ijv_idx = 0;
    }

    // add hessian+id*eps_id
    double eps = 1e-8;
    for (int i = 0; i < hessian_IJV.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(hessian_IJV[i].row(),hessian_IJV[i].col(),-hessian_IJV[i].value());
    }
    for (int i = 0; i < var_n; i++) { kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,-eps);}

    // Add both J and J transpose
    for (int i = 0; i < jacobian_IJV.size(); i++) {
        int row = jacobian_IJV[i].row(), col = jacobian_IJV[i].col(); double val = jacobian_IJV[i].value();
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(col,row+var_n, val); // Jt columed offseted at var_n
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(row+var_n, col, val); // J offseted in var_n
    }
    // add zeros along var_n+i,var_n_i
    for (int i = var_n; i < var_n+const_n; i++) kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,0);

    if ( A.rows() == 0) {
      A =  Eigen::SparseMatrix<double>(var_n+const_n,var_n+const_n);
      igl::sparse_cached_precompute(kkt_IJV, cached_ijv_data, A);
    } else {
      igl::sparse_cached(kkt_IJV, cached_ijv_data, A);
    }
}
void NewtonKKT::build_kkt_system(const Eigen::SparseMatrix<double>& hessian,
                        const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& KKT) {
    Eigen::SparseMatrix<double> id (hessian.rows(),hessian.cols()); id.setIdentity();
    auto eps_id = 1e-8*id;

    Eigen::SparseMatrix<double> H = -hessian - eps_id;
    
    Eigen::SparseMatrix<double> Jt = J.transpose();

    Eigen::SparseMatrix<double> H_jt; igl::cat(2,H,Jt, H_jt);
    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    igl::cat(1, H_jt, J_0, KKT);

    //A.makeCompressed();
    Eigen::SparseMatrix<double> id_KKT(KKT.rows(),KKT.rows()); id_KKT.setIdentity(); id_KKT = 0*id_KKT;
    KKT = KKT + id_KKT; // todo: stupid but Paradiso wants to add zeros explicitly
}

double NewtonKKT::one_iter(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x,
        double current_merit) {
    x = x0;
    int vnum = x.rows()/3;
    double new_e;

    igl::Timer timer; auto init_time = timer.getElapsedTime(); auto t = init_time;
    // Get Hessian
    //auto hessian = f.hessian(x); 
    auto hessian_ijv = f.update_and_get_hessian_ijv(x);
    auto hessian_time = timer.getElapsedTime()-t;

    // Get Jacobian
    t = timer.getElapsedTime();
    //auto jacobian = constraints.Jacobian(x);
    auto jacobian_ijv = constraints.update_and_get_jacobian_ijv(x);
    auto jacobian_time = timer.getElapsedTime()-t;

    t = timer.getElapsedTime();
    //Eigen::SparseMatrix<double> A; build_kkt_system(hessian,jacobian,A);
    int var_n = x.rows(); int const_n = constraints.getConstNum();
    build_kkt_system_from_ijv(hessian_ijv, var_n,
                        jacobian_ijv, const_n);

    auto kkt_time = timer.getElapsedTime()-t;  

    t = timer.getElapsedTime();

    //energy->check_grad(x);
    double old_e = f.obj(x);
    Eigen::VectorXd g(f.grad(x)); Eigen::VectorXd neg_g = -1*g;
    Eigen::VectorXd d(g.rows()); Eigen::VectorXd lambda_d(lambda.rows());
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    //cout << "analayzing pattern" << endl;
    //solver.analyzePattern(A);
    // cout << "NewtonKKT A.norm() = " << A.norm() << endl;
    if (first_solve) {
        m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
        m_solver.set_pattern();
        m_solver.analyze_pattern();
        //m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
        //m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
        first_solve = false;
    } else {
        m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    auto analyze_pattern_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    m_solver.factorize();
    /*solver.factorize(A);
    if(solver.info()!=Eigen::Success) {
        cout << "Eigen Failure!" << endl;
        exit(1);
    }*/
    auto factorize_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    
    //m_solver.factorize();

    Eigen::VectorXd constraints_deviation = -1*constraints.Vals(x);
    Eigen::VectorXd g_const; igl::cat(1, g, constraints_deviation, g_const);
    //g_const = -1*g_const;
    
    Eigen::VectorXd res;
    //cout << "solving!" << endl;
    //res = solver.solve(g_const);
    m_solver.solve(g_const,res);

    //cout << "Checking precision: (A*res-g_const).norm() = " << (A*res-g_const).norm() << endl;
    
    for (int d_i = 0; d_i < g.rows(); d_i++) {d[d_i] = res[d_i];}
    for (int d_i = g.rows(); d_i < res.rows(); d_i++) {lambda_d[d_i-g.rows()] = res[d_i];}
    auto solve_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    double step_size = 1;
    //new_e = line_search(x,d,init_t,f);
    new_e = exact_l2_merit_linesearch(x,d,step_size,f,constraints,current_merit);
    lambda = lambda + step_size*lambda_d;
    auto linesearch_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();


    double total_time = timer.getElapsedTime()-init_time;
    /*
    cout << endl << endl << "total_kkt_system_time  = " << total_time << endl;
    cout << "Hessian_compute_time  = " << hessian_time << endl;
    cout << "Jacobian_compute_time  = " << jacobian_time << endl;
    cout << "kkt_system_build_time  = " << kkt_time << endl;
    
    cout << "hessian analyze_pattern_time = " << analyze_pattern_time << endl;
    cout << "factorize_time = " << factorize_time << endl;
    cout << "solve_time  = " << solve_time << endl;
    cout << "linesearch_time  = " << linesearch_time << endl;
    */
    //std::ofstream out_file(std::string("KKT_mat"));
    //out_file << igl::matlab_format(A,"A");
    //out_file << igl::matlab_format(g_const,"g_const");
    //std::ofstream out_file(std::string("null_space.m"));
    //out_file << igl::matlab_format(Jt,"Jt");
    
    //std::cout << "NewtonKKT: old_e = " << old_e << " new_e = " << new_e << std::endl;
    
    return new_e;
}

double NewtonKKT::kkt_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints) {
    auto jacobian = eq_constraints.Jacobian(x);
    
    auto g = obj.grad(x);
    // TODO this is wrong, add the lambda (even if we don't use it for optimization)
    double grad_error = 0.0001*g.norm();//0.0001*(g-jacobian.transpose()*lambda).norm();
    double eq_const_error = eq_constraints.Vals(x).norm();
    
    //std::cout << "\tgrad_error = " << grad_error << std::endl;
    //std::cout << "\teq_const_error = " << eq_const_error << std::endl;
    return (std::max)(grad_error, eq_const_error);
    //return eq_const_error;
}