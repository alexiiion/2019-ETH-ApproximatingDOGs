#pragma once

#include "MeshModel.h"
#include "Rulings.h"

//#include "Spline.h"
/*
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
	Eigen::VectorXd curvature;
	Eigen::VectorXd torsion;

	Eigen::MatrixXd tangents;
	Eigen::MatrixXd principal_normals;
	Eigen::MatrixXd binormals;

	virtual void compute_curve() = 0;	
	virtual Eigen::MatrixXd compute_ruling_directions() = 0;

	virtual void print_mathematica_data() {};

protected:
	Mesh& target;
	Eigen::MatrixXd& geodesic;
};
*/
class RuledDevelopableSurface
{
public:
	Mesh developable;

	RuledDevelopableSurface(Mesh& target, Eigen::MatrixXd& geodesic, const int at_vertex, const RulingsType& rulings_type);
	~RuledDevelopableSurface();

	Mesh create(const double width = 5.0);
	Eigen::MatrixXd generator_curve();
	int get_width();
	void print_mathematica_data();

private:
	Eigen::MatrixXd& geodesic;
	int vertex;
	double width;

	RulingsType rulings_type;
	Rulings* rulings_controller;

	//Eigen::MatrixXd rulings;
	//Eigen::VectorXd ruling_lengths;
	//Eigen::Vector3d direction_normalizer;

	Mesh create_mesh(const Eigen::MatrixXd& rulings, const Eigen::VectorXd& ruling_lengths, const double width);

	//Eigen::Vector3d normalize_rulings_side(Eigen::MatrixXd& rulings, int max_k_index);
	//Eigen::VectorXd compute_ruling_lengths(Eigen::MatrixXd& rulings);
	
	////double compute_ruling_length(Eigen::RowVector3d tangent, Eigen::RowVector3d& ruling_current, Eigen::RowVector3d& ruling_next, const Eigen::RowVector3d& ruling_first, const double step);
	//double compute_ruling_length(const Eigen::RowVector3d& tangent, const Eigen::RowVector3d& ruling_current, const Eigen::RowVector3d& ruling_next, const Eigen::RowVector3d& ruling_direction, const double step);

	//double compute_ruling_length(const int index, const int other, const Eigen::MatrixXd& rulings, const Eigen::RowVector3d& tangent, const double step);
	//double compute_ruling_length_bidirectional(const int index, const int previous, const int next, const Eigen::MatrixXd& rulings, const double step);
};
