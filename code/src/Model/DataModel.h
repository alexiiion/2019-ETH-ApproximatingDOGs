#pragma once

#include <igl/serialize.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include "OptimizationSettings.h"
#include "GeodesicFeatureConfig.h"

#include "GeodesicCandidates.h"
#include "GeodesicsDirection.h"
#include "RuledDevelopableSurface.h"
#include "MeshModel.h"

#include "DevelopableModel.h"
#include "DevelopableOptimizationSettings.h"

//#include "GeodesicsModel.h"
//#include "RuledSurfacesModel.h"
#include "NormalDeviationStopping.h"

#include "GlobalSettings.h"
#include "Coverage.h"
#include "ApproximationState.h"

//forward declarations
class Patch; 


class DataModel
{

public:
	//DataModel()
	//	: geodesics(GeodesicsModel(target)), ruled_surfaces(RuledSurfacesModel(target, geodesics))
	//{
	//	update_flat_threshold();
	//};


	/*
	VIEW only
	*/
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;

	static int log_level_optimization;

	//
	static double output_smoothness;
	static double output_bending_energy;
	static double output_isometry_energy;
	static double output_regularizer_energy;
	static double output_fitting_energy;
	static double output_laplacian_energy;
	static double output_bilaplacian_energy;

	static double output_fitting_objective;
	static double output_DOG_objective;
	static bool is_local_minimum;
	static double current_constraints_step;

	Eigen::VectorXi face_labels; // #target.F
	Eigen::VectorXi vertex_labels; // #target.V


	/* --- editor settings  -- */

	int pre_selected_vertex = -1;

	bool is_pausing_on_adding_patch = false;
	bool is_pausing_on_changing_constraint_strategy = false;
	bool is_rendering_patch_coverage = false;


	/* --- optimization -- */

	OptimizationSettings optimization_settings;

	//optimization status
	bool is_optimizing = false;
	bool is_optimizing_DOG = false;
	bool is_optimizing_mesh = false;

	int iterations_since_converge = 0;


	ApproximationState state;


	/* --- target -- */

	std::string models_folder;
	std::string target_filename;

	float target_max_dimension = 20;
	float target_diagonal;
	Eigen::RowVectorXd target_dimensions;

	Mesh target;
	Mesh target_original;

	Coverage* patch_coverage = NULL;
	Coverage* result_coverage = NULL;
	Coverage* initialization_coverage = NULL;


	/* --- result -- */

	Developable developable_model;
	//DevelopableOptimizationSettings developable_optimization_settings;


	/* --- DOGs -- */

	std::vector<Patch*> patches;
	std::vector<int> wrapper_view_indices;
	Eigen::MatrixXd concatenated_patches_V;
	Eigen::MatrixXi concatenated_patches_F;


	/* --- holes in coverage --- */

	std::vector<Mesh> hole_meshes;
	std::vector<std::vector<int>> hole_vertex_indices;
	std::vector<int> hole_average_vertex;


	/* --- geodesics -- */
	

	std::vector<int> geodesics_directions; //in degree, where k_max = 0° and k_min = 90°
	GeodesicCandidates geodesics_candidates; //debug: store longest of geodesics_directions (so 1 geodesic per vertex)
	int current_candidates_index = -1;

	//int stopping_window_size = 5;
	//float stopping_outlier_factor = 2.0;
	NormalDeviationStopping geodesics_stopping;
	bool use_geodesics_flat_detetction;

	/* --- ruled developable surfaces -- */

	std::vector<int> ruled_vertex_indices;
	std::vector<int> selected_ruled_vertex_indices; 
	std::vector<RuledDevelopableSurface*> ruled_developables; //sparse vector of size #V, only has entries at ruled_vertex_indices 
	std::vector<int> vertex_ruled_assignment; // #target.V
	
	int ruled_view_index;
	std::vector<int> selected_ruled_view_indices;

	//float spline_smoothness = 0.05;
	float ruled_width = 15;
	int number_random_points = 100;

	int label_selection_smoothness = 50;



	/* --- REFACTORING -- */
	
	//TargetModel target;

	//ResultModel result; (Coverage, K)

	//PatchesModel patches; (Coverage, K)
	//HoleModel holes;

	//RuledSurfacesModel ruled_surfaces; (Coverage, K)
	//GeodesicsModel geodesics;


	//GeodesicsModel& geodesics;
	//RuledSurfacesModel& ruled_surfaces;

	//GeodesicsModel* geodesics;
	//RuledSurfacesModel* ruled_surfaces;




	void update_flat_threshold()
	{
		GlobalSettings::flat_curvature_threshold = 1.0 / (target_max_dimension * 0.5);
	}

};

