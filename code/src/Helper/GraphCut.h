#pragma once

#include <vector>
#include "MeshModel.h"

namespace graphcut
{
	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness = 50);
	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness, std::vector<double>& out_resulting_distances);

	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, const std::vector<Eigen::MatrixXd>& patch_vertex_normals, double smoothness, std::vector<double>& out_resulting_distances);

	std::vector<int> compute_graph_cut_face_labeling(const Eigen::MatrixXd& face_midpoints, const Eigen::MatrixXi& connectivity, const std::vector<Eigen::VectorXd>& patch_face_distances, double smoothness, std::vector<double>& out_resulting_distances);


	/*
	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness = 50);
	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness, std::vector<double>& out_resulting_distances);
	
	std::vector<int> compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, const std::vector<Eigen::VectorXd>& patch_vertex_normals, double smoothness, std::vector<double>& out_resulting_distances);

	std::vector<int> compute_graph_cut_face_labeling(const Eigen::MatrixXd& face_midpoints, const Eigen::MatrixXi& connectivity, const std::vector<Eigen::VectorXd>& patch_face_distances, double smoothness, std::vector<double>& out_resulting_distances);

	void data_to_matrix(const std::vector<Eigen::VectorXd>& data, Eigen::MatrixXd& matrix);
	void normalize_data(Eigen::MatrixXd& data, const double& max);
	void transform_data_exponentially(Eigen::MatrixXd& data, const double& multiplier, const double& cutoff);
	*/
}
