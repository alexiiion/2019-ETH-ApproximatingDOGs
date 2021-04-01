#pragma once

#include "MeshModel.h"
#include "OptimizationSettings.h"

struct BlockConstraintsConfig
{
	int index;

	bool is_in_use;
	bool has_constant_pairing;
	bool is_removing_outliers;
	bool is_removing_normal_outliers;
	bool is_limiting_face_constraints;

	int pairing_direction;
	int pairing_target;

	Eigen::MatrixXd points_from;
	Eigen::MatrixXd points_to;

	Eigen::MatrixXi barycentric_indices;
	Eigen::MatrixXd barycentric_weights;

	std::vector<bool> filtered_constraints;
};


class PositionConstraints
{

public:
	PositionConstraints(Mesh& target, WrapperMesh& wrapper, OptimizationSettings& optimization_settings);

	Mesh& target;
	WrapperMesh& wrapper;
	OptimizationSettings& optimization_settings;

	Eigen::MatrixXd paired_points_from;
	Eigen::MatrixXd paired_points_to;
	Eigen::MatrixXi constrained_barycentric_indices;
	Eigen::MatrixXd constrained_barycentric_weights;

	Eigen::MatrixXd* pair_from_V;
	Eigen::MatrixXi* pair_from_F;
	Eigen::MatrixXd* pair_to_V;
	Eigen::MatrixXi* pair_to_F;	


	//updating settings
	//float* outlier_threshold;
	bool is_keeping_pairs_constant = true;
	bool should_remove_outliers = true;
	bool should_remove_normal_outliers = false;
	bool should_limit_face_constraints = true;
	bool should_foster_stretch = false;

	bool use_previous_constraints = false; 
	bool use_all_constraints = false; // not tested

	int pairing_direction;
	int pairing_target;


	//void initialize(Mesh* target, WrapperMesh* wrapper, Eigen::MatrixXd& from, Eigen::MatrixXd& to);
	//void initialize(Mesh* target, WrapperMesh* wrapper, OptimizationSettings* optimization_settings);
	void update_positions();

	void add(Eigen::MatrixXd& from);
	void set_current_constraints(Eigen::MatrixXd& from);
	void remove_all_constraints();

	void get_flattend_constraints(Eigen::MatrixXi& out_wrapper_barycentric_indices_flat, Eigen::MatrixXd& out_wrapper_barycentric_weights_flat, Eigen::VectorXd& out_target_points_flat);

	void set_pairing_direction(int pairing_direction);
	void invert_pairing_direction();

private:

	int number_all_constraints = 0;
	int current_block = -1;
	std::vector<BlockConstraintsConfig> blocks; 
	

	void add(Eigen::MatrixXd& from, Eigen::MatrixXd& to);
	void update_block(BlockConstraintsConfig& block);
	void construct_optimized_constraints();

	void label_outliers(BlockConstraintsConfig& block);
	void label_normal_outliers(BlockConstraintsConfig& block);
	void label_exceeding_face_constraints(BlockConstraintsConfig& block);

	void add_block_settings(int block_size);

	void find_closest_points(const Eigen::MatrixXd& from, Eigen::MatrixXd& out_to);
	void find_closest_points(const Eigen::MatrixXd& from, Eigen::MatrixXd& out_to, Eigen::VectorXi& out_closest_faces);
};
