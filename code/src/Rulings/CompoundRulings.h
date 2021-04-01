#pragma once

#include <Eigen/Core>
#include "AnalyticRulings.h"

class CompoundRulings : public AnalyticRulings
{
public:
	CompoundRulings(Mesh& target, Eigen::MatrixXd& geodesic);
	~CompoundRulings();

	//virtual void compute() override;
	virtual void print_mathematica_data() override;

protected:
	double sample_step;

	Eigen::MatrixXd analytic_ruling_directions;
	Eigen::MatrixXd analytic_ruling_lengths;

	//Eigen::MatrixXd surface_normals;
	Eigen::MatrixXd discrete_tangents, discrete_normals, discrete_binormals;
	Eigen::VectorXd discrete_curvature;
	std::vector<std::pair<int, int>> filter_bins;

	virtual void compute_curve() override;
	virtual Eigen::MatrixXd compute_ruling_directions() override;
	//virtual Eigen::VectorXd compute_ruling_lengths(Eigen::MatrixXd& rulings) override;

	//filter: false = invalid vertex, analytic ruling should *not* be used
	std::vector<bool> build_discrete_curvature_filter();
	std::vector<bool> build_ruling_lengths_filter();
	std::vector<bool> build_tangent_proximity_filter();

	Eigen::MatrixXd get_discrete_rulings();
	Eigen::VectorXd get_discrete_curvature();
	Eigen::MatrixXd get_interpolated_rulings(const int start_index, const int end_index, const Eigen::MatrixXd& analytic_ruling_directions);
	Eigen::RowVector3d get_discrete_ruling_at(const int curve_index, const Eigen::MatrixXd& analytic_ruling_directions);

	void compute_discrete_descriptors(const Eigen::MatrixXd& sampled_curve, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals);
	void compute_discrete_descriptors();
	Eigen::MatrixXd get_discrete_binormals();
	//std::vector<std::pair<int, int>> flat_parts;
	//std::vector<std::pair<int, int>> get_flat_parts(const Eigen::VectorXd& discrete_curvature);

	//Eigen::MatrixXd get_interpolated_rulings(const int number_rulings, const Eigen::RowVector3d& start_ruling, const Eigen::RowVector3d& end_ruling);
	//Eigen::RowVector3d get_discrete_binormal_at(int curve_index);
	//Eigen::MatrixXd get_binormals_as_rulings(int start, int end);
};
