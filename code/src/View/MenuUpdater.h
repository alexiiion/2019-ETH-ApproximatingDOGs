#pragma once

#include "DataModel.h"
#include "ViewModel.h"
#include "InteractiveSegmentSelector.h"

class OptimizationController;
class RuledGeodesicsController;

class MenuUpdater
{
public:
	MenuUpdater(DataModel& data_model, ViewModel& view_model) : data_model(data_model), view_model(view_model) {};

	void update_editor_menu();
	void update_mode_menu();

	void update_wrapper_mesh_menu();
	void update_target_mesh_menu(OptimizationController& optimization);
	void update_optimization_menu(OptimizationController& optimization);
	void update_ruled_geodesics_menu(RuledGeodesicsController& ruled_geodesics);

	void update_preparataion_menu();
	void update_surface_analysis_menu();
	void update_features_menu(InteractiveSegmentSelector& segment_selection);
	void update_geodesics_menu();

	void mesh_view_options(MeshViewSettings& settings);
	void patch_view_options(PatchViewSettings& settings);

	void update_debug_window_menu(OptimizationController& optimization);

	void activate_elements();
	void deactivate_elements();

private:
	DataModel& data_model;
	ViewModel& view_model;

};

std::string get_curvature_range_info(const Eigen::VectorXd& curvature);