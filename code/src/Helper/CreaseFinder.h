#pragma once

#include <Eigen/Core>
#include <vector>

namespace crease 
{

	void creases_ring_neighbors(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<std::vector<int>>& adjacency_VF, int number_rings, Eigen::VectorXd& out_creases);
	void creases_geodesic_neighbors(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double geodesic_distance, Eigen::VectorXd& out_creases);

}