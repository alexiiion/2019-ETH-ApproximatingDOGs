#pragma once

#include <vector>
#include <Eigen/Dense>
#include <igl/PI.h>

#include "GlobalSettings.h"
//#include "MeanCluster.h"
#include "GeodesicsStopping.h"
#include "MeshModel.h"

/*
enum WalkingStop
{
	NormalsDeviation,
	FrameDeviation,
	MaxLength,
	AdaptiveNormalsDeviation,
};
*/


class GeodesicWalker
{
public:
	GeodesicWalker(Mesh& target, GeodesicsStopping& stopping_criteria) : target(target), stopping_criteria(stopping_criteria)
	{ };

	/*
	GeodesicWalker(Mesh& target) : target(target) 
	{
		stopping_mode = WalkingStop::AdaptiveNormalsDeviation;
	};

	GeodesicWalker(Mesh& target, float* max_frame_error, bool* use_weighted_frame_error) : target(target), max_frame_error(max_frame_error), use_weighted_frame_error(use_weighted_frame_error) 
	{
		stopping_mode = WalkingStop::FrameDeviation;
	};
	GeodesicWalker(Mesh& target, float max_normals_deviation_deg) : target(target)
	{
		stopping_mode = WalkingStop::NormalsDeviation;
		max_normals_deviation = max_normals_deviation_deg / 180 * igl::PI;
	}; 
	*/

	~GeodesicWalker() {};

	GeodesicWalker& operator = (const GeodesicWalker &t)
	{
		return *this;
	}

	void get_geodesic_kmax_at(const int vertex_index, Eigen::MatrixXd& out_geodesic, double& out_length);
	void get_geodesic_kmin_at(const int vertex_index, Eigen::MatrixXd& out_geodesic, double& out_length);
	void get_geodesic_at(const int vertex_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_geodesic, double& out_length);
	void get_geodesic_at(const int vertex_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_geodesic, double& out_length, Eigen::MatrixXi& out_path_adjacent_vertices, Eigen::VectorXi& out_path_faces);


private:
	Mesh& target;
	int source_vertex;

	GeodesicsStopping& stopping_criteria;
	int added_segments = 0;
	

	double begin_crease_threshold = GlobalSettings::crease_angle_threshold;

	/*
	WalkingStop stopping_mode;
	//parameters for stopping criteria

	//for WalkingStop::FrameDeviation
	float* max_frame_error;
	bool* use_weighted_frame_error;

	//for WalkingStop::MaxLength
	float* max_length;
	double current_length = 0.0;

	//for WalkingStop::NormalsDeviation
	double max_normals_deviation = 1.0 / 180 * igl::PI;
	double previous_normals_deviation = 0.0;
	double current_normals_deviation = 0.0;
	
	//for WalkingStop::AdaptiveNormalsDeviation
	MeanCluster cluster;
	double last_deviation = 0.0;
	int fid_previous = -1;
	int fid_current = -1;

	//errors
	double error_previous = 0.0;
	double error_derivative = 0.0;
	double error_integral = 0.0;


	bool is_done_walking();
	void update_stopping_criteria(const Eigen::MatrixXd& path);
	void reset_stopping_criteria_data();
	*/


	void trace_geodesic(Eigen::Vector3d& direction, Eigen::MatrixXd& path, double& length, Eigen::MatrixXi& path_adjacent_vertices, Eigen::VectorXi& path_faces);
	
	double frame_error(const int vertex, const int vertex_compare, bool use_weighted_error = false);
	int get_source_face(const Eigen::Vector3d& direction, const Eigen::Vector3d& e0, const Eigen::Vector3d& e1);

	Eigen::MatrixXd merge_parts(const Eigen::MatrixXd& part_positive, const Eigen::MatrixXd& part_negative);
	Eigen::MatrixXi merge_parts(const Eigen::MatrixXi& part_positive, const Eigen::MatrixXi& part_negative);
	Eigen::VectorXi merge_parts(const Eigen::VectorXi& part_positive, const Eigen::VectorXi& part_negative);
};
