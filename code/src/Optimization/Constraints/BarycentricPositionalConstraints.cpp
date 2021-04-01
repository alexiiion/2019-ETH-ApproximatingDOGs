#include "BarycentricPositionalConstraints.h"

#include "DataModel.h"
#include "Logger.h"

Eigen::VectorXd BarycentricPositionalConstraints::Vals(const Eigen::VectorXd& x_wrapper) const 
{
	//Eigen::VectorXd edgeCoords(EdgePoint::getPositionInMesh(edgePoints, x));
	Eigen::VectorXd current_barycenters_location(const_n);
	for (int i = 0; i < const_n; i++)
	{
		double point_wrapper = 0;

		for (int bary_i = 0; bary_i < 3; bary_i++)
		{
			double weight = barycentric_weights_wrapper(i, bary_i);
			int vertex_index = barycentric_indices_wrapper(i, bary_i);
			double vertex = x_wrapper(vertex_index);

			point_wrapper += weight * vertex;
		}
		current_barycenters_location(i) = point_wrapper;
	}

	Eigen::VectorXd distances = current_barycenters_location - bc;
	
	DataModel::output_fitting_energy = distances.norm();
	write_log(DataModel::log_level_optimization) << "    BarycentricPositionalConstraints norm: " << DataModel::output_fitting_energy << std::endl;

	return distances;
}

void BarycentricPositionalConstraints::updateJacobianIJV(const Eigen::VectorXd& x)
{
	int vn = x.rows() / 3;
	int const_cnt = 0; int ijv_idx = 0;
	for (int i = 0; i < const_n; i++) {

		for (int bary_i = 0; bary_i < 3; bary_i++) {
			double weight = barycentric_weights_wrapper(i, bary_i);
			int vertex_index = barycentric_indices_wrapper(i, bary_i);
			//IJV.push_back(Eigen::Triplet<double>(const_cnt, vertex_index, weight));
			IJV[ijv_idx++] = Eigen::Triplet<double>(const_cnt, vertex_index, weight);
		}
		const_cnt++;
	}
}

//virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
//	int const_n = 0;
//	for (int b_i = 0; b_i < b.rows(); b_i++) {
//		int var_const_idx = b(b_i);
//		// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
//		IJV[b_i] = Eigen::Triplet<double>(const_n++, var_const_idx, 1);
//	}
//}

