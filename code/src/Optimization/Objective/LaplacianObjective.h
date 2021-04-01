#pragma once

#include <igl/repdiag.h>

#include "Objective.h"

#include "DataModel.h"
#include "Logger.h"


class LaplacianObjective : public Objective
{
public:
	LaplacianObjective(const Eigen::SparseMatrix<double>& L_i)
	{
		igl::repdiag(L_i, 3, L);
		//IJV.resize(L.nonZeros());

		IJV = to_triplets(L);

		////Eigen::SparseMatrix<double> mat(rows, cols);
		//int ijv_cnt = 0;
		//for (int k = 0; k < L.outerSize(); ++k) 
		//{
		//	for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
		//		IJV[ijv_cnt++] = Eigen::Triplet<double>(it.row(), it.col(), it.value());
		//	}
		//}
	};

	virtual LaplacianObjective* clone() const { return new LaplacianObjective(*this); }

	virtual double obj(const Eigen::VectorXd& x) const
	{
		double objective = x.transpose() * L * x;

		write_log(DataModel::log_level_optimization) << "    LaplacianObjective obj: " << objective << std::endl;

		return objective;
	}

	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const
	{
		Eigen::VectorXd gradient = 2 * L*x;
		write_log(DataModel::log_level_optimization) << "    LaplacianObjective grad.norm(): " << gradient.norm() << std::endl;

		return gradient;
	}


private:
	Eigen::SparseMatrix<double> L;

	virtual void updateHessianIJV(const Eigen::VectorXd& x)
	{
		/* empty on purpose */
	}
};