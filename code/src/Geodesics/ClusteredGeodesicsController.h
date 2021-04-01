//#pragma once
//
//#include "DataModel.h"
//#include "GeodesicsController.h"
//
//class ClusteredGeodesicsController : public GeodesicsController
//{
//public:
//	ClusteredGeodesicsController(DataModel& data_model) : GeodesicsController(data_model) {};
//
//	virtual bool initialize() override;
//	virtual bool get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve) override;
//
//
//private:
//	std::vector<int> random_geodesics;
//
//	void build_geodesics();
//	void cluster_geodesics();
//
//	void get_random_ruled_developables(int n);
//	void get_ruled_developable(const Eigen::MatrixXd& geodesic);
//	std::vector<int> get_random_geodesics(int n);
//
//	bool select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve);
//	void get_cluster_at(const int& start_index, std::vector<bool>& visited, std::vector<int>& out_cluster_vertices);
//
//	void compute_geodesics_descriptors();
//	void distances_between_geodesics(const Eigen::MatrixXd& path, const Eigen::MatrixXd& path_neighbor, const Eigen::RowVector3d& edge_offset, double& out_distance_mean, double& out_distance_variance, double& out_distance_sum);
//	void distances_between_geodesics(const Eigen::MatrixXd& path, const Eigen::MatrixXd& path_neighbor, const Eigen::RowVector3d& edge_offset,
//		Eigen::MatrixXd& out_distance_vectors, Eigen::VectorXd& out_shortest_distances, Eigen::VectorXd& out_distances_1order,
//		double& out_distance_mean, double& out_distance_variance, double& out_distance_sum);
//};
