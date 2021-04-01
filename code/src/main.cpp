//#define EIGEN_USE_MKL_ALL

//#include <igl/readOBJ.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/upsample.h>

#include <algorithm>
#include <iostream>
#include <fstream>

#include "View.h"
#include "DataModel.h"
#include "OptimizationController.h"
#include "MeshController.h"

#include "Serializer.h"
#include "Logger.h"
#include "Utils.h"

#include "GlobalSettings.h"
#include "Helper/ConfigParser.h"

//#include "GeodesicsController.h"
//#include "CoordinateConverter.h"
#include "GeodesicWalker.h"
//


//#include "GeodesicsModel.h"
//#include "RuledSurfacesModel.h"
//#include "Initialization/GeodesicSurfacesController.h"


#define show_coordinate_indicator 0


using namespace std;

/* static variables */
int LOG_LEVEL = 4;

int DataModel::log_level_optimization = 6;

double DataModel::output_smoothness = 0.0;

double DataModel::output_bending_energy = 0.0;
double DataModel::output_isometry_energy = 0.0;
double DataModel::output_regularizer_energy = 0.0;
double DataModel::output_fitting_energy = 0.0;
double DataModel::output_laplacian_energy = 0.0;
double DataModel::output_bilaplacian_energy = 0.0;

double DataModel::output_fitting_objective = 0.0;
double DataModel::output_DOG_objective = 0.0;
bool   DataModel::is_local_minimum = false;
double DataModel::current_constraints_step = 0.0;

double GlobalSettings::flat_curvature_threshold = 0.1; // max curvature (i.e. 1 / osculating circle radius) to consider something flat is dependent on the mesh size 
double GlobalSettings::flat_angle_threshold = 0.035;   // angle; corresponds to  2.0°
double GlobalSettings::crease_angle_threshold = 0.7;   // angle; corresponds to 40.1°
double GlobalSettings::max_average_edge = 1.2;   // angle; corresponds to 40.1°



double GlobalSettings::rulings_spline_smoothness = 0.1;
double GlobalSettings::rulings_spline_sample_factor = 3; //number of samples per avg_edge_length
int    GlobalSettings::rulings_spline_degree = 5;
double GlobalSettings::min_tangent_ruling_deviation = 45.0 / 180.0 * M_PI;
double GlobalSettings::min_ruling_length = 1;



//bool is_loading_scene = false;
//bool is_loading_target = false;
//bool is_starting_optimization = false;
//
//string scene_file;




/* controllers and models */
DataModel data_model;

//GeodesicsController geodesics_controller;
//OptimizationController optimization(data_model, geodesics_controller);
OptimizationController optimization(data_model);
View view(data_model, optimization);



/* Pretend this is the ApplicationController */
//GeodesicSurfacesController surface_initialization(data_model.target, data_model.geodesics, data_model.ruled_surfaces);





void start_menu(float width, float offset_x, float offset_y, string menu_name)
{
	float menu_width = width; // *data_model.menu.menu_scaling();
	ImGui::SetNextWindowPos(ImVec2(offset_x, offset_y), ImGuiSetCond_FirstUseEver);
	//ImGui::SetNextWindowSize(ImVec2(menu_width, 0.0f), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
	bool _viewer_menu_visible = true;

	ImGui::Begin(menu_name.c_str(), &_viewer_menu_visible, ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
}

void end_menu()
{
	ImGui::PopItemWidth();
	ImGui::End();
}

bool callback_update_view(igl::opengl::glfw::Viewer& viewer)
{
	optimization.run_optimiztaion();
	return view.update_view();
}

void callback_update_menu()
{
	start_menu(400.f, 0.0f, 0.0f, "Viewer");

	view.update_menu();

	end_menu();
}

void callback_update_debug_menu()
{
	start_menu(300.f, 400.0f, 0.0f, "Debug");

	view.update_debug_menu();

	end_menu();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers)
{
	return view.callback_key_down(key);
}

bool callback_key_up(igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers)
{
	return view.callback_key_up(key);
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	return view.callback_mouse_down(button);
}

bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y)
{
	return view.callback_mouse_move(mouse_x, mouse_y);
}

bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
	return view.callback_mouse_up(button);
}

void callback_draw_custom_window()
{
	view.update_debug_menu();
}

void initialize_view()
{
	// Attach a menu plugin
	data_model.viewer.plugins.push_back(&data_model.menu);

	//data_model.viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 1);
	data_model.viewer.core().background_color << 1.0f, 1.0f, 1.0f, 1.0f;
	data_model.viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

	//data_model.viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
	//data_model.viewer.core().rotation_type = igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL;
	//data_model.viewer.data().point_size = 10.0;
	//data_model.viewer.data().line_width = 1.0;
	data_model.viewer.resize(1600, 1400);
	//data_model.viewer.resize(1920*0.65, 1080*0.9); //for recording

#if show_coordinate_indicator
	const auto coordinate_indicator = Eigen::MatrixXd::Identity(3, 3);
	data_model.viewer.data().add_edges(Eigen::MatrixXd::Zero(3, 3), coordinate_indicator*0.2, coordinate_indicator);
#endif
	view.view_model.debug_view_index = data_model.viewer.selected_data_index;

	view.initialize();
}

void viewer_register_callbacks()
{
	data_model.viewer.callback_pre_draw = callback_update_view; // calls at each frame
	//data_model.menu.callback_draw_viewer_menu = callback_update_menu;
	data_model.menu.callback_draw_viewer_window = callback_update_menu;
	data_model.menu.callback_draw_custom_window = callback_update_debug_menu;

	data_model.viewer.callback_key_down = callback_key_down;
	data_model.viewer.callback_key_up = callback_key_up;
	data_model.viewer.callback_mouse_down = callback_mouse_down;
	data_model.viewer.callback_mouse_move = callback_mouse_move;
	data_model.viewer.callback_mouse_up = callback_mouse_up;

	data_model.viewer.core().is_animating = true;
	data_model.viewer.core().animation_max_fps = 30.;
}

int main(int argc, char* argv[])
{
	//by default, load mesh as is
	GlobalSettings::max_average_edge = 1e6; //prevent upsampling 
	data_model.target_max_dimension = 0;	//prevent scaling

	StartConfig config;
	config.read_startup_config(get_folder_path(__FILE__) + "config.txt", data_model);

	view.initialize();


	if (config.is_loading_scene)
	{
		deserialize(config.scene_file, data_model, view.view_model, true);
		config.set_mesh_dependent_settings(data_model);

		if (config.is_starting_optimization)
			optimization.initialize();
	}
	else if (config.is_loading_target)
	{
		meshhelper::add_target(data_model.target_filename, data_model, view.view_model.target_view_settings.view_index, view.view_model.target_developable_view_settings.view_index);
		config.set_mesh_dependent_settings(data_model);

		if (config.is_starting_optimization)
			optimization.initialize();
	}
	else //start empty, load via UI
	{
		data_model.models_folder = "";
		view.view_model.target_view_settings.standard_fill_color = Eigen::RowVector3d(179.0 / 255.0, 176.0 / 255.0, 174.0 / 255.0);
	}

	viewer_register_callbacks();
	data_model.viewer.launch();
}