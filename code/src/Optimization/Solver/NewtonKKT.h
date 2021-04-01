#pragma once

#include "Solver.h"
#include "PardisoSolver.h"

// This just used the linearized constraints but the lagrangian but solves [H,J^t;J,0] 
//  meaning H is just the obj Hessian without the second order parts of the constraints
class NewtonKKT : public ConstrainedSolver {
  
public:
	NewtonKKT(const double& infeasability_epsilon, const double& infeasability_filter, const int& max_newton_iters, const double& merit_p) :
			 infeasability_epsilon(infeasability_epsilon), infeasability_filter(infeasability_filter), 
			 max_newton_iters(max_newton_iters) ,merit_p(merit_p), m_solver(ai,aj,K) {m_solver.set_type(-2);}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x,
			double convergence_threshold = 1e-4);

private:
	double one_iter(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x, double current_merit);
	void build_kkt_system(const Eigen::SparseMatrix<double>& hessian, const Eigen::SparseMatrix<double>& Jacobian,
						Eigen::SparseMatrix<double>& KKT);

	void build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, int var_n,
                                     const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n);

	double kkt_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints);

	const double& infeasability_epsilon; 
	const double& infeasability_filter;
	const int& max_newton_iters;
	const double& merit_p;

	std::vector<Eigen::Triplet<double> > kkt_IJV; Eigen::VectorXi cached_ijv_data;
	Eigen::SparseMatrix<double> A;

	PardisoSolver m_solver;

	bool first_solve = true;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
	Eigen::VectorXd lambda;
};