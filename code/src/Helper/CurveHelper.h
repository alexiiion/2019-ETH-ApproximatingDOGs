#pragma once

#include <Eigen/Dense>
#include "MeshModel.h"

namespace curve
{
	Eigen::MatrixXd remove_short_segments(const Eigen::MatrixXd& path, const double min_length = 0.01);
	
	//removes straight & short segments
	Eigen::MatrixXd remove_straight_segments(const Eigen::MatrixXd& path, const double epsilon_angle = 0.996, const double min_length = 0.01, const double max_length = 1.5);
	
	double compute_length(const Eigen::MatrixXd& path);

	bool get_frame_at(const int i, const Eigen::MatrixXd& path, Eigen::RowVector3d& out_T, Eigen::RowVector3d& out_N, Eigen::RowVector3d& out_B);
	Eigen::Matrix3d get_frame_at(int i, const Eigen::MatrixXd& path);
	void compute_frames(const Eigen::MatrixXd& path, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals);
	void compute_surface_frames(const Eigen::MatrixXd& sampled_curve, const Mesh& mesh, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals);

	double compute_curvature(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
	Eigen::VectorXd compute_curvature(const Eigen::MatrixXd& path);
	std::vector<double> compute_normal_deviations(const Eigen::MatrixXd& path, const Mesh& target);

	bool is_flat(const Eigen::MatrixXd& path);
	bool is_flat(const Eigen::MatrixXd& path, const double flat_threshold);
	bool is_flat(std::vector<double>& path_curvature);
	bool is_flat(std::vector<double>& path_curvature, const double flat_threshold);

	Eigen::MatrixXd smooth_curve(const Eigen::MatrixXd& path, const int iterations = 20, const double smooth_rate = 0.5, const double inflate_rate = -0.55);
	Eigen::MatrixXd smooth_normals(const Eigen::MatrixXd& normals, const int iterations = 20, const double smooth_rate = 0.5, const double inflate_rate = -0.55);
	
	void resample_uniformly(const Eigen::MatrixXd& original_path, const double segment_length, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
	void resample_uniformly(const Eigen::MatrixXd& original_path, Eigen::MatrixXd& out_sampled, double& out_sampled_length);

	void resample_matching_curve(const Eigen::MatrixXd& to_sample, const int center_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);

	void sample_forward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
	void sample_backward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);

	//TODO remove! it's legacy
	double compute_curvature_sum(const Eigen::MatrixXd& path);
	double compute_gauss_curvature_sum(Mesh& target, const Eigen::MatrixXd& path);
}