#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"

#include "DataModel.h"
#include "Logger.h"


class PointPairConstraints : public Constraints 
{
public:
	//pairs are of size #V*dim, 2
	PointPairConstraints(const Eigen::MatrixXi& pairs) : pairs(pairs)
	{
		// Here the pairs are assumed to be defined on the vector pairs (so pairing up points means having 3 index pairs in the c'tor here)
		const_n = pairs.rows();
		IJV.resize(2*const_n);
	};

	virtual PointPairConstraints* clone() const 
	{ 
		return new PointPairConstraints(*this);
	}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const 
	{
		Eigen::VectorXd vals(pairs.rows());
		for (int i = 0; i < pairs.rows(); i++)
			vals(i) = x(pairs(i,0)) - x(pairs(i, 1));

		write_log(DataModel::log_level_optimization) << "    PointPairConstraints norm: " << vals.norm() << std::endl;

		return vals;
	}

	virtual void updateJacobianIJV(const Eigen::VectorXd& x) 
	{
		int const_n = 0; 
		int ijv_cnt = 0;
		
		for (int b_i = 0; b_i < pairs.rows(); b_i++ )
		{
			int var_const_idx1 = pairs(b_i, 0); 
			int var_const_idx2 = pairs(b_i, 1);

			// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_n, var_const_idx1, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_n, var_const_idx2, -1);
			const_n++;
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) 
	{
		// Linear constraints have zero second derivative. Empty on purpose
	};

	Eigen::MatrixXi pairs;
};