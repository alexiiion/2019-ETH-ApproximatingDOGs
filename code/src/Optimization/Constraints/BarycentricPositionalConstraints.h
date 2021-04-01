#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"
#include "Quad.h"

class BarycentricPositionalConstraints : public Constraints {
public:
	//BarycentricPositionalConstraints(const Eigen::MatrixXi& barycentric_indices_wrapper, const Eigen::MatrixXd& barycentric_weights_wrapper, const Eigen::VectorXd& mapped_points_target)
	//	: barycentric_indices_wrapper(barycentric_indices_wrapper), barycentric_weights_wrapper(barycentric_weights_wrapper), bc(mapped_points_target) {
	//		const_n = bc.rows(); approx_nnz = 3 * const_n;
	//	}
	BarycentricPositionalConstraints(const Eigen::MatrixXi& barycentric_indices_wrapper, const Eigen::MatrixXd& barycentric_weights_wrapper, const Eigen::VectorXd& mapped_points_target)
		: barycentric_indices_wrapper(barycentric_indices_wrapper), barycentric_weights_wrapper(barycentric_weights_wrapper), bc(mapped_points_target) {
		const_n = barycentric_indices_wrapper.rows();
		IJV.resize(const_n * 3);
	};

	BarycentricPositionalConstraints() {const_n = 0; bc.resize(0);} // empty set of constraints c'tor

	virtual BarycentricPositionalConstraints* clone() const {return new BarycentricPositionalConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;

	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) 
	{
		// Linear constraints have zero second derivative. Empty on purpose
	};

	void update_coords(const Eigen::VectorXd& mapped_points_target) 
	{ 
		//std::cout << "(bc-mapped_points_target).norm() = " << (bc - mapped_points_target).norm() << std::endl;
		//int wait; std::cin >> wait;
		bc = mapped_points_target; 
		//std::cout << "(bc-mapped_points_target).norm() = " << (bc - mapped_points_target).norm() << std::endl;
		//std::cin >> wait;
	}

	//Eigen::VectorXi getPositionIndices() const { return b; }
	//Eigen::VectorXd getPositionVals() const { return bc; }

	//std::vector<EdgePoint> getEdgePoints() const {return edgePoints;}
	//Eigen::VectorXd getEdgePointConstraints() const {return bc;}

private:

	Eigen::VectorXd bc;

	Eigen::MatrixXi barycentric_indices_wrapper;
	Eigen::MatrixXd barycentric_weights_wrapper;
};