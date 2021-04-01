#pragma once

#include <Eigen/Core>
#include <vector>

#include "View/Model/ViewModes.h"
#include "View/Model/MeshVisualization.h"
#include "View/Model/MeshViewSettings.h"

#include "MeshModel.h"
#include "RuledDevelopableSurface.h"

#include <igl/Timer.h>

//enum InteractionMode
//{
//	Navigate,
//	Select,
//};
//
//enum SelectionMode
//{
//	NotSelecting,
//	DefineFeatures,
//	GeodesicsToAll,
//	SpecificGeodesic,
//	WalkingGeodesic
//	//MeshPicking
//};

//enum MeshVisualization
//{
//	Standard,
//	InitializationAssignment,
//	//ResultAssignment,
//	PrincipalCurvature,
//	GaussianCurvature,
//	MeanCurvature,
//	//WeightedCurvature
//};

//struct MeshViewSettings
//{
//	int view_index = -1;
//	bool has_changed = true;
//
//	bool show_mesh = true;
//	bool show_faces = true;
//	bool show_wireframe = true;
//	bool label_faces;
//	bool label_vertices;
//
//	Eigen::RowVectorXd dimensions;
//
//	MeshVisualization visualization;
//	Eigen::MatrixXd fill;
//};
//
//struct PatchViewSettings : MeshViewSettings
//{
//	bool label_quad_faces;
//	bool show_position_constraints = true;
//	bool show_target_constraints = true;
//};


struct ViewModel
{
	ViewModel()
	{
		//current_position_constraints.color_start = Colors::BLUE;
	};


	int debug_view_index = 0;


	/* --- modes -- */

	InteractionMode current_interaction_mode;
	SelectionMode current_selection_mode;
	bool is_mouse_down = false;

	//int selected_vertex = -1;
	//bool do_update_selected_vertex = false;
	//bool do_show_vertex = false;

	Eigen::MatrixXd selected_points;
	bool do_update_points = false;

	MeshVisualization overlay_visualization = MeshVisualization::Standard; //standard == none
	OverlaySettings overlay_settings;
	bool has_overlay_changed = false;

	///* --- curvatures --- */

	//bool do_update_surface_labels = true;

	float visualized_min = -0.005;
	float visualized_max = 0.005;

	bool show_filtered_curvature_points = false;
	bool show_filtered_curvature_values = false;
	float filter_curvature = 1.0;


	/* --- target --- */

	MeshViewSettings target_view_settings;
	//TargetVisualization current_target_visualization = TargetVisualization::Standard;
	//double target_scale = 1.0;
	int current_coverage_choice = -1;
	bool show_coverage = true;
	bool show_coverage_points = false;
	bool show_coverage_holes = false;
	//bool show_coverage_coloring = false;

	MeshViewSettings target_developable_view_settings;
	MeshViewSettings developable_parts_view_settings;
	int selected_developable_part = -1;

	/* --- patches --- */

	PatchViewSettings general_patch_view_settings;
	std::vector<PatchViewSettings> patch_view_settings;
	Eigen::RowVector3d patch_color = Colors::YELLOW; //Colors::GRAY_MID;
	bool use_general_patch_settings = false;
	bool show_only_current_patch = false;
	bool do_update_all_patches = false;
	int selected_patch_index = -1;

	
	/* --- gauss map --- */

	Mesh gauss_map_sphere;
	double gauss_map_sphere_radius = 5;
	Mesh* gauss_mapped_mesh;
	int gauss_map_view_index = -1;
	bool is_gauss_map_visible;
	bool do_update_gauss_map;
	int choice_gauss_map;


	/* --- creases --- */

	bool show_creases = false;
	bool label_crease_value = false;
	bool label_crease_vertex = false;
	float angle_crease_threshold = 0.0;
	Eigen::VectorXd crease_vertices;
	float normalized_creases_threshold = 1.0;
	Eigen::VectorXd crease_vertices_normalized;


	/* --- points of interest --- */
	bool do_render_POI = true;


	/* --- geodesics --- */
	   
	bool has_geodesic_selection_changed = false;
	bool do_update_geodesics = false;

	//walking geodesics
	bool do_trace_kmin = false;
	bool do_trace_kmax = true;
	std::vector<int> walking_geodesics;
	//bool show_debug_walking_geodesics_neighbors = false;
	//bool do_render_geodesics_clusters = true;
	//int render_geodesics_cluster_index = -1;

	bool is_selecting_debug_geodesic = false;
	bool do_render_debug_walking_geodesics = false;
	std::vector<Eigen::MatrixXd> debug_walking_geodesics;

	//geodesics from selected source to all other vertices
	int merit_max_index;

	int selected_geodesic_source_vertex = -1;
	bool do_render_graph_geodesic_paths = false;
	bool do_render_graph_geodesic_distances = false;
	std::vector<Eigen::MatrixXd> all_geodesic_paths_from_source;
	Eigen::VectorXd all_geodesic_lengths_from_source;

	//geodesics between two selected vertices
	bool do_render_specific_geodesics = false;
	std::vector<Eigen::MatrixXd> selectd_geodesic_paths;
	std::vector<double> selectd_geodesic_lengths;
	std::vector<int> paired_geodesic_sources;
	std::vector<int> paired_geodesic_targets;





	/* --- ruled developables --- */
	//MeshViewSettings ruled_developables_view_settings;
	int selected_ruled_index = -1;

	bool show_all_geodesics = false;
	bool show_random_geodesics = false;
	bool show_selected_geodesics = true;

	bool show_points_geodesics = false;
	bool show_curves_geodesics = true;
	bool show_surfaces_geodesics = true;

	bool do_update_init_geodesics = false;
	bool are_init_surfaces_visible = false;


	bool show_debug_ruled_developables = false;
	bool show_selected_ruled_developables = false;
	bool print_developables_data = true;
	//std::vector<int> ruled_view_indices;
	//float ruled_width = 5.0;

	//int number_random_points = 100;
	//std::vector<int> random_vertices;

	////spline
	//float spline_smoothness = 1.0;
	//float spline_sampling = 0.01;



	/* --- debug items --- */

	bool show_debug_items = true;

	int ruled_debug_view_index = -1;
	bool is_debugging_ruled_surfaces = false;
	bool do_update_debug_ruled_surfaces = false;
	int use_debug_rulings_type = 0; //0=auto, 1=analytic, 2=discrete, 3=compound
	int selected_debug_ruled_surface_index = -1;
	RuledDevelopableSurface* debug_ruled_surface = NULL;
	//bool has_debug_ruled_surface_changed = false;

	int geodesics_view_index = -1;
	int debug_geodesics_view_index = -1;
	int init_geodesics_view_index = -1;
	bool show_all_geodesics_candidates = false;
	bool show_geodesics_candidates = false;
	bool show_geodesics_candidates_labels = false;

	bool show_feature_cluster_graph = false;
	int selected_feature_cluster = -1;
	std::vector<std::vector<Eigen::MatrixXd>> graph_geodesic_paths;
	std::vector<std::vector<double>> graph_geodesic_lengths;
	std::vector<std::vector<int>> graph_geodesic_sources;
	std::vector<std::vector<int>> graph_geodesic_targets;




	//weights for weighted curvature (currently unused)
	float w_K = 0.9;
	float w_H = 1.0 - w_K;



	//OverlayData current_position_constraints;
	////current_position_constraints.color_start = Colors::BLUE;
	//
	//OverlayData current_target_curve;
	//OverlayData all_target_curves;
	//OverlayData coverage;
	//
	//OverlayData all_clustered_features;
	//
	//OverlayData geodesic_paths;
	//OverlayData preview_filtered_feature_points;
	//OverlayData preview_selected_feature_points;#

	//variables for animated mesh rotation
	Eigen::RowVector3d mid;
	igl::Timer timer;
	bool is_rotating = false;
	int rotation_axis = 0;
	int animation_speed_flag = 4; // 0 to 4 for different speeds
};

//struct OverlayData
//{
//	Eigen::MatrixXd start_points;
//	Eigen::MatrixXd end_points;
//
//	Eigen::Vector3d color_start;
//	Eigen::Vector3d color_end;
//
//	bool show_points;
//	bool show_edges;
//	bool interpolate_colors;
//
//	bool is_valid();
//};
