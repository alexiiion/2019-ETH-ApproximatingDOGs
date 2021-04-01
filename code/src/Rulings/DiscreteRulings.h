#pragma once

#include <Eigen/Core>

#include "Rulings.h"

class DiscreteRulings : public Rulings
{
public:
	DiscreteRulings(Mesh& target, Eigen::MatrixXd& geodesic);
	~DiscreteRulings();

	virtual void print_mathematica_data() override;

protected:
	double sample_step;
	bool is_flat = false;

	virtual void compute_curve() override;
	virtual Eigen::MatrixXd compute_ruling_directions() override;

	void compute_frames(const Eigen::MatrixXd& sampled_curve, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals);
	Eigen::MatrixXd compute_surface_normals(const Eigen::MatrixXd& sampled_curve);
	Eigen::MatrixXd get_binormals_as_rulings();
};
