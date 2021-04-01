#pragma once

#include "Objective.h"

#include "DataModel.h"
#include "Logger.h"

class ProximityObjective : public Objective
{
public:
	ProximityObjective(const Eigen::VectorXd& x_target) : x_target(x_target)
	{
		IJV.resize(x_target.rows());
	};

	virtual ProximityObjective* clone() const { return new ProximityObjective(*this); }

	virtual double obj(const Eigen::VectorXd& x) const
	{
		double objective = (x - x_target).squaredNorm();
		write_log(DataModel::log_level_optimization) << "    ProximityObjective obj: " << objective << std::endl;

		return objective;
	}

	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const
	{
		Eigen::VectorXd gradient = 2 * (x - x_target);
		return gradient;
	}


private:
	Eigen::VectorXd x_target;

	virtual void updateHessianIJV(const Eigen::VectorXd& x)
	{
		/* ??? empty on purpose ??? */
	}
};