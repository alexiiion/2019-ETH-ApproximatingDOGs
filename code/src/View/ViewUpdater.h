#pragma once

#include "DataModel.h"
#include "ViewModel.h"
#include "InteractiveSegmentSelector.h"

class ViewUpdater
{
public:
	ViewUpdater(DataModel& data_model, ViewModel& view_model) : data_model(data_model), view_model(view_model) {};

	void update_target_view(MeshViewSettings& mesh_settings, const Mesh& mesh);
	void update_developables_parts_view();
	void update_patches_view();
	void update_gauss_map_view();
	void update_geodesics_init_view();

	void update_surface_properties_view();
	void update_feature_view(InteractiveSegmentSelector& segment_selection);
	void update_geodesics_view();
	void update_debug_view();
	void update_debug_geodesics_view();
	void update_result_view();
	//void update_test_view();
	
	void update_target_fill();
	void update_mesh_fill(MeshViewSettings& mesh_settings, Mesh& mesh);
	void update_mesh_color_assignment_fill(const Mesh& mesh, const std::vector<int>& label_assignment, const int number_labels);

private:
	DataModel& data_model;
	ViewModel& view_model;

	//void update_mesh_fill(const Mesh& mesh, const TargetVisualization& visualization, const std::vector<int>& label_assignment = {});

	void update_coverage_view();
	void update_mesh_view(const MeshViewSettings& mesh_settings);
	void update_single_patch_view(PatchViewSettings& view_settings, Patch& patch);

	void visualize_curvature(const Eigen::VectorXd& curvature, const Eigen::VectorXi& curvature_lookup);
	void show_filtered_curvature_points(const Eigen::VectorXd& feature, const Eigen::VectorXi& feature_lookup);
	void set_standard_fill();

	//void update_crease_view();
	//void update_cluster_view();
	//void update_graph_geodsics();
	//void show_geodesics_at_selected_vertex();
	//void update_all_geodesics_from_vertex(int vertex_index);
	//void show_specific_geodesics();
	//void update_specific_geodesics();


	//void update_patch_fill(const Patch& patch, MeshViewSettings& mesh_settings);

	/*
	void render_highlighted_curves(const std::vector<Eigen::MatrixXd>& curves, const int selected_index, const Eigen::RowVector3d& color, const bool show_points);
	void render_highlighted_curves(const std::vector< std::vector<Eigen::MatrixXd>>& curves, const int selected_index, const Eigen::RowVector3d& color, const bool show_points);

	void show_curve(const Eigen::MatrixXd& curve, const Eigen::RowVector3d& color, bool show_points = false);
	void show_curve(const std::vector<int>& curve_indices, const Eigen::MatrixXd& lookup, const Eigen::RowVector3d& color, bool show_points = false);

	void label_quad_faces(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F_quad);
	*/
};
