//#pragma once
//
//#include "DataModel.h"
//#include "GeodesicsController.h"
////#include "GeodesicCandidates.h"
//
//class FeatureGeodesicsController : public GeodesicsController
//{
//public:
//	FeatureGeodesicsController(DataModel& data_model) : GeodesicsController(data_model) {};
//
//	virtual bool initialize() override;
//	virtual bool get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve) override;
//	
//	void update_geodesics_merit(int index);
//
//
//private:
//	bool is_feature_covered = false;
//
//	std::vector<int> used_feature_neighborhoods;
//
//	void build_feature_graph();
//	void build_geodesics();
//
//	bool select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve);
//	void compute_dynamic_features(int index);
//	void compute_constant_features(const GeodesicFeatureConfig& feature_settings, GeodesicCandidates& geodesic_candidates);
//	void add_normalized_constant_feature(const int current_feature_index, const Eigen::VectorXd& feature, Eigen::MatrixXd& out_features);
//
//	GeodesicCandidates get_all_geodesics(DataModel& data_model, std::vector<int> sources, std::vector<std::vector<int>> targets);
//	void find_geodesics_on_target(const TargetMesh& target, std::vector<int>& sources, std::vector<int>& targets, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths);
//	void preprocess_target_geodesics(const TargetMesh& target, GeodesicCandidates& geodesic_candidates);
//
//	void resample_target_geodesic(int index, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
//	void resample_wrapper_geodesic(const Eigen::MatrixXd& to_sample, const int center_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
//	void sample_forward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
//	void sample_backward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length);
//};
