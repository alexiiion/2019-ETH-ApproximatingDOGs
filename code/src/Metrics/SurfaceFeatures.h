#pragma once

#include <Eigen/Dense>

class Mesh;

class SurfaceFeatures
{
public:
	//Eigen::VectorXd weighted_curvature;
	//Eigen::VectorXi weighted_curvature_lookup;

	//gaussian curvature
	Eigen::VectorXd K;
	Eigen::VectorXd K_abs;
	Eigen::VectorXd K_abs_normalized;
	Eigen::VectorXi K_abs_lookup;

	//principal curvature
	Eigen::VectorXd principal_k_min;
	Eigen::VectorXd principal_k_max;
	Eigen::MatrixXd principal_k_min_direction;
	Eigen::MatrixXd principal_k_max_direction;
	Eigen::VectorXd Kp; //Gauss curvature computed form principal curvatures
	Eigen::VectorXi Kp_abs_lookup;

	//mean curvature
	Eigen::VectorXd H;
	Eigen::VectorXd H_normalized;
	Eigen::VectorXi H_lookup;

	//dihedral angle (creases)
	Eigen::VectorXd crease_vertices;
	Eigen::VectorXd crease_vertices_normalized;
	Eigen::VectorXi crease_vertices_lookup;

	bool are_built = false;
	void build(Mesh& mesh);
	void reset();

private:
	//void compute_features(Mesh& target);
	//void compute_weighted_features(Mesh& target, double weight_K, double weight_crease, double weight_H);

	void compute_principal_k(Mesh& mesh);
	void compute_K(Mesh& mesh);
	void compute_H(Mesh& mesh);
	void compute_creases(Mesh& mesh);

	void get_gaussian_curvature(Mesh& mesh, Eigen::VectorXd& out_K);
	void get_mean_curvature(Mesh& mesh, Eigen::VectorXd& out_H);
	void normalize_vector_entries(const Eigen::VectorXd& curvature, Eigen::VectorXd& out_normalized_curvature);
	void dihedral_triangle_angles(const Mesh& mesh, Eigen::VectorXd& out_creases);
};
