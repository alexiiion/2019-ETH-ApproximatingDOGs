#pragma once

#include "DataModel.h"
#include "ViewModel.h"

#include "ViewUpdater.h"
#include "MenuUpdater.h"

#include "InteractiveSegmentSelector.h"
#include "OptimizationController.h"


class View
{

	//TODO add hashing to check if this was added already https://ideone.com/tieHbd

public:
	View(DataModel& data_model, OptimizationController& optimization) : data_model(data_model), optimization(optimization), view_updater(ViewUpdater(data_model, view_model)), menu_updater(MenuUpdater(data_model, view_model)), segment_selection(InteractiveSegmentSelector(data_model)) {};
	void initialize();
	
	ViewModel view_model;
	InteractiveSegmentSelector segment_selection;

	bool update_view();
	bool pre_draw(igl::opengl::glfw::Viewer& viewer);
	void update_menu();
	void update_debug_menu();
	
	bool callback_key_down(unsigned int key);
	bool callback_key_up(unsigned int key);
	bool callback_mouse_down(int button);
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up(int button);


private:
	DataModel& data_model;
	OptimizationController& optimization;

	ViewUpdater view_updater;
	MenuUpdater menu_updater;

	bool is_initialized = false;
	bool has_changed = false;

	int get_vertex_from_screen(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	void style_viewer();


	//MeshRenderer
	//SurfaceRenderer
	//SegmentRenderer
	//OptimizationSettingsRenderer
	
	//not view
	//SurfaceAnalyzer 
	//SegmentSelector
};
