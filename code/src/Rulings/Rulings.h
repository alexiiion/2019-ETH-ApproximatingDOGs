#pragma once

#include "MeshModel.h"

enum RulingsType
{
	Analytic,
	Discrete,
	Compound
};

class Rulings
{

public:
	Rulings(Mesh& target, Eigen::MatrixXd& geodesic) : target(target), geodesic(geodesic) {};
	~Rulings() {};

	Eigen::MatrixXd sampled_curve;
	Eigen::MatrixXd ruling_directions;
	Eigen::VectorXd ruling_lengths;

	virtual void compute();
	virtual void print_mathematica_data();

protected:
	Mesh& target;
	Eigen::MatrixXd& geodesic;

	Eigen::VectorXd curvature;
	Eigen::VectorXd torsion;

	Eigen::MatrixXd tangents;
	Eigen::MatrixXd principal_normals;
	Eigen::MatrixXd binormals;
	Eigen::Vector3d direction_normalizer;

	virtual void compute_curve() = 0;
	virtual Eigen::MatrixXd compute_ruling_directions() = 0;
	virtual Eigen::VectorXd compute_ruling_lengths(Eigen::MatrixXd& rulings);

	Eigen::Vector3d normalize_rulings_side(Eigen::MatrixXd& rulings);

	double compute_ruling_length(const Eigen::RowVector3d& tangent, const Eigen::RowVector3d& ruling_current, const Eigen::RowVector3d& ruling_next, const Eigen::RowVector3d& ruling_direction, const double step);
	double compute_ruling_length(const int index, const int other, const Eigen::MatrixXd& rulings, const Eigen::RowVector3d& tangent, const double step);

};
