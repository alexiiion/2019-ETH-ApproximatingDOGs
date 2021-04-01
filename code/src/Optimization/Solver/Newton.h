#pragma once

#include <igl/Timer.h>

#include "Solver.h"
#include "PardisoSolver.h"
#include "line_search.h"

#include "Logger.h"

class Newton : public Solver 
{

	//TODO add version for using Eigen Solver, for easier access and publishing code!

public:
	Newton() : m_solver(ai, aj, K) {
		m_solver.set_type(2);
	}


	bool first_solve = true;
	PardisoSolver m_solver;
	Eigen::VectorXd lambda;
	Eigen::VectorXi ai, aj;
	Eigen::VectorXd K;


	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) 
	{
        x = x0;
        
		//igl::Timer timer;
		//double t = timer.getElapsedTime();

		//Eigen::SparseMatrix<double> id(x.rows(),x.rows()); id.setIdentity();
		//Eigen::SparseMatrix<double> H = obj.hessian(x) + 1e-2*id;
		////Eigen::SparseMatrix<double> H = 1e-2*id;
        Eigen::SparseMatrix<double> H = obj.hessian(x);
		/*
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(H);*/
		if (first_solve) {
			m_solver.set_system_matrix(H.triangularView<Eigen::Upper>());
			m_solver.set_pattern();
			m_solver.analyze_pattern();
			first_solve = false;
		}
		else {
			m_solver.update_system_matrix(H.triangularView<Eigen::Upper>());
		}
		m_solver.factorize();

		//write_log(4) << timer.getElapsedTime() - t << " --> duration for solver.compute(H)" << std::endl;
		//t = timer.getElapsedTime();
		/*
        if(solver.info() != Eigen::Success) 
		{
            std::cout << "Eigen Failure! error code: " << solver.info() << std::endl;
            exit(1);
        }
		*/
        
		//auto d = -1*solver.solve(obj.grad(x));
		Eigen::VectorXd rhs = -1 * obj.grad(x);
		Eigen::VectorXd d; m_solver.solve(rhs, d);
        double init_step_size = 1;

		//write_log(4) << timer.getElapsedTime() - t << " --> duration for solver.solve(obj.grad(x))" << std::endl;
		//t = timer.getElapsedTime();

        double new_e = line_search(x, d, init_step_size, obj);

		//write_log(4) << timer.getElapsedTime() - t << " --> duration for line_search" << std::endl << std::endl;
		//t = timer.getElapsedTime();

        return new_e;
    }
};