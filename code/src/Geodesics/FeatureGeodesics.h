#pragma once

class DataModel;

#include "MeshModel.h"
#include "GeodesicFeatureConfig.h"


class FeatureGeodesics
{
public:
	FeatureGeodesics() {};

	//feature vectors for all geodesics
	Eigen::VectorXd geodesics_merit; //unsorted (i.e. index is geodesic path index), column entries are properties
	Eigen::MatrixXd geodesic_features;
	//std::vector<int> used_geodesics;

	//all geodesics to select from
	std::vector<Eigen::MatrixXd> paths;
	//Eigen::VectorXd coverage;

	//geodesic features
	Eigen::VectorXd lengths;
	Eigen::VectorXi lengths_lookup; //in descending order

	Eigen::VectorXd gauss_curvatures;
	Eigen::VectorXi gauss_curvatures_lookup; //in descending order

	Eigen::VectorXd curvatures;
	Eigen::VectorXi curvatures_lookup; //in descending order

	int number_geodesics;
	int current_geodesic_index;
	int current_merit;

//	void get_all_geodesics(DataModel& data_model, std::vector<int> sources, std::vector<std::vector<int>> targets);
//	void get_all_geodesics(DataModel& data_model, std::vector<std::vector<int>>& feature_connectivity_graph);
//
//	void find_geodesics_on_target(const TargetMesh& target, std::vector<int>& sources, std::vector<int>& targets, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths);
//
//private:
//
//	//void get_all_geodesics(DataModel& data_model, std::vector<int>& segment_indices, int source_index, FeatureGeodesics& segment_geodesics);
//	//void find_geodesics_on_target(DataModel& data_model, std::vector<int>& segment_indices, int& source_index, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths);
//	void preprocess_target_geodesics(const TargetMesh& target);
//
//	void compute_constant_features(const GeodesicFeatureConfig& feature_settings);
//	void add_normalized_constant_feature(const int current_feature_index, const Eigen::VectorXd& feature, Eigen::MatrixXd& out_features);
//
//	void remove_straight_segments(const Eigen::VectorXd& geodesic_lengths, std::vector<Eigen::MatrixXd>& out_paths);
//	void compute_gauss_curvatures(const TargetMesh& target, const std::vector<Eigen::MatrixXd>& paths, const Eigen::VectorXd& lengths, Eigen::VectorXd& out_gaussian_curvatures);
//	void compute_curvatures(const std::vector<Eigen::MatrixXd>& paths, Eigen::VectorXd & out_curvatures);
//	void compute_curvature(const Eigen::MatrixXd& path, double& out_curvature);
};
