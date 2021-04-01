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



bool is_loading_scene = false;
bool is_loading_target = false;
bool is_starting_optimization = false;

string scene_file;




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

void string_to_vector(const string& input, vector<int>& out_list, const char deliminator = ' ')
{
	std::stringstream iss(input);

	int number;
	out_list.clear();
	while (iss >> number)
		out_list.push_back(number);
}

void read_config()
{
	write_log(2) << endl << "READING config settings: " << endl << endl;

	//get config file path
	string file = __FILE__;
	auto index = file.find_last_of("/\\");
	auto path = file.substr(0, index + 1);
	auto config_filepath = path + "config.txt";

	// Open the File
	std::ifstream in(config_filepath.c_str());

	// Check if object is valid
	if (!in)
	{
		std::cerr << "Cannot open the config file : " << config_filepath << std::endl;
		exit(1);
	}

	std::string line;
	// Read the next line from File untill it reaches the end.
	while (std::getline(in, line))
	{
		if (line.find('#') == 0)
			continue;
		if (line.empty())
			continue;

		auto colon_index = line.find_first_of(":");
		string qualifier = line.substr(0, colon_index);
		string content = line.substr(colon_index + 2, line.length() - colon_index + 2);
		
		if (content.empty())
			continue;
		
		if (qualifier == "start_empty")
		{
			is_loading_scene = false;
			is_loading_target = false;

			write_log(3) << "starting empty " << endl;
			write_log(3) << "-- stopping reading config here." << endl;
			break;
		}
		else if (qualifier == "models_folder")
		{
			data_model.models_folder = content.c_str();
			write_log(3) << "folder: " << data_model.models_folder << endl;
		}
		else if (qualifier == "load_scene")
		{
			scene_file = content;
			is_loading_scene = true;

			write_log(3) << "load_scene: " << scene_file << endl;
			write_log(3) << "-- stopping reading config here." << endl;
			break;
		}
		else if (qualifier == "initialize_optimization")
		{
			is_starting_optimization = atoi(content.c_str());
			write_log(3) << "initialize_optimization: " << boolalpha << is_starting_optimization << endl;
		}
		else if (qualifier == "target")
		{
			is_loading_target = true;
			data_model.target_filename = content;
			write_log(3) << "target_filename: " << data_model.target_filename << endl;
		}
		else if (qualifier == "scale_target")
		{
			data_model.target_max_dimension = atoi(content.c_str());
			data_model.update_flat_threshold();
			write_log(3) << "scale_target: " << data_model.target_max_dimension << endl;
		}
	}

	//Close The File
	in.close();
	write_log(2) << endl << "DONE reading config settings. " << endl << endl;
}

//bool is_nan(double d) { return std::isnan(d); }

int main(int argc, char* argv[])
{
	read_config();
	view.initialize();


	if (is_loading_scene)
	{
		GlobalSettings::max_average_edge = 1e6; //prevent upsampling

		deserialize(scene_file, data_model, view.view_model, true);
		//deserialize(scene_file, data_model, view.view_model);
		//optimization.geodesics_controller->initialize();

		
		//------ REVISION EXPERIMENTS -------

		////Random:: Bumpy side 

		//double target_size = data_model.target_dimensions.maxCoeff();
		////data_model.number_random_points = data_model.target.V.rows() * .20;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.ruled_width = target_size * 0.5;
		//data_model.optimization_settings.coverage_threshold = target_size * 0.035;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.25;
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;

			
		////Random:: face mask-aimshape
		////size = 40
		//data_model.label_selection_smoothness = 100;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.number_random_points = data_model.target.V.rows() * .1;
		//data_model.ruled_width = data_model.target_max_dimension * 0.3;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.035;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.1; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;
		//data_model.optimization_settings.weight_bending_energy = 10.0;

		
		
		
		
		if (is_starting_optimization)
			optimization.initialize();

		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.is_optimizing_DOG = false;
		//data_model.is_optimizing_mesh = true;









		//data_model.number_random_points = data_model.target.V.rows() * 1;

		////meshhelper::update_target_developable_optimized(data_model);
		
		//view.view_model.use_general_patch_settings = true;
		//view.view_model.general_patch_view_settings.show_mesh = false;
		//view.view_model.target_view_settings.show_mesh = false;
		//view.view_model.target_developable_view_settings.show_mesh = false;


		////for screencastst:
		//data_model.is_optimizing_DOG = true;
		//data_model.is_optimizing_mesh = false;

		//view.view_model.use_general_patch_settings = false;
		//view.view_model.target_developable_view_settings.show_mesh = false;
		//view.view_model.show_coverage = false;
		//view.view_model.show_coverage_points = false;

		////GlobalSettings::flat_curvature_threshold = 0.05;


		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;

		//data_model.geodesics_directions = { 0 };

		////bumpy
		////data_model.ruled_vertex_indices = { 3017 };
		////data_model.ruled_vertex_indices = { 44, 3281, 3576, 3017 };
		////data_model.ruled_vertex_indices = { 3017, 317, 3926, 34, 3576, 3604, 1702, 4296, 178, 1796, 1815 };
		////data_model.ruled_vertex_indices = { 1473, 2280, 2368, 269, 4391, 3567, 3134, 187, 3916, 400  };
		//data_model.ruled_vertex_indices = { 1473, 2368, 269, 4391, 3567, 3134, 187, 3916, 400  };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//
		//data_model.number_random_points = data_model.target.V.rows() * .2;

		//data_model.optimization_settings.outlier_threshold = 0.45 - data_model.target.average_edge_length;




		//// lilium  edge length: 0.41339
		//data_model.ruled_vertex_indices = { 670, 1496, 3038, 2739, 745 };
		////data_model.ruled_vertex_indices = { 1779, 562, 2453, 2365, 2739, 2858, 1379, 2158, 1562  };
		////data_model.ruled_vertex_indices = { 1779, 16, 514, 446, 1882, 3282, 2684, 2403, 802, 1409, 2369, 2672, 2739, 1914, 3037, 2514, 3091, 3145, 199, 614, 1442, 728, 710, 718, 3006, 2223, 1562, 1775, 320, 3025 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;
		//data_model.optimization_settings.outlier_threshold = 0.2; //~ edgelength * 0.5
		////data_model.optimization_settings.coverage_threshold = 0.41; //~ edgelength * 1.0
		//data_model.optimization_settings.coverage_threshold = 0.8; //~ edgelength * 2.0
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_directions = { 0 };
		////data_model.geodesics_stopping.window_size = 10;
		////data_model.geodesics_stopping.bounds_min_absolute = 0.05;

		////data_model.is_pausing_on_adding_patch = true;
		////data_model.is_pausing_on_changing_constraint_strategy = true;



		//// bumpy-side-irregular  edge length: 0.599774
		//data_model.ruled_vertex_indices = { 803, 222, 964, 1290, 2303 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;
		//data_model.optimization_settings.outlier_threshold = 0.2; //~ edgelength * 0.33
		//data_model.optimization_settings.coverage_threshold = 0.6; //~ edgelength * 1.0
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_directions = { 0 };
		//data_model.geodesics_stopping.window_size = 10;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;

		////data_model.is_pausing_on_adding_patch = true;
		////data_model.is_pausing_on_changing_constraint_strategy = true;

	}
	else if (is_loading_target)
	{
		GlobalSettings::max_average_edge = 1e6; //prevent upsampling
		//GlobalSettings::max_average_edge = 2; 
		//GlobalSettings::max_average_edge = 0.1; 
		

		//data_model.target_max_dimension = 0; //prevent scaling

		meshhelper::add_target(data_model.target_filename, data_model, view.view_model.target_view_settings.view_index, view.view_model.target_developable_view_settings.view_index);
		//meshhelper::add_new_target(data_model.target_filename, data_model, view.view_model.target_view_settings.view_index, view.view_model.target_developable_view_settings.view_index);




		//view.view_model.overlay_visualization = MeshVisualization::GaussianCurvature;
		//
		//int min_index;
		//int max_index;
		//double K_min = data_model.target.surface_features().K.minCoeff(&min_index);
		//double K_max = data_model.target.surface_features().K.maxCoeff(&max_index);
		//write_log(0) << "K range:  " << K_min << " - " << K_max << " at " << min_index << ", " << max_index << endl;
		//
		//double P_min_min = data_model.target.surface_features().principal_k_min(min_index);
		//double P_min_max = data_model.target.surface_features().principal_k_max(min_index);
		//double PK_min = P_min_min * P_min_max;
		////write_log(0) << "P at min: min = " << P_min_min << ", max = " << P_min_max << " --> K = " << P_min_min * P_min_max << endl;
		//
		//
		//double Kp_min = data_model.target.surface_features().Kp.minCoeff(&min_index);
		//double Kp_max = data_model.target.surface_features().Kp.maxCoeff(&max_index);
		//write_log(0) << "Kp range: " << Kp_min << " - " << Kp_max << " at " << min_index << ", " << max_index << endl;








		data_model.optimization_settings.weight_bending_energy = 6.0;
		data_model.optimization_settings.convergence_developablity_threshold = 0.01;
		
		data_model.ruled_width = 10;
		data_model.number_random_points = data_model.target.V.rows() * .05;
		data_model.geodesics_directions = {0, 90};
		//data_model.geodesics_directions = {0, 45, 90};


		data_model.label_selection_smoothness = 500;
		data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.0;
		data_model.optimization_settings.coverage_threshold = 0.8;
		


		//------ REVISION EXPERIMENTS -------

		////Random:: Bumpy side 

		//double target_size = data_model.target_dimensions.maxCoeff();
		//data_model.number_random_points = data_model.target.V.rows() * .20;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.ruled_width = target_size * 0.5;
		//data_model.optimization_settings.coverage_threshold = target_size * 0.035;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.25;
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		////Random:: lilium --> doesn't work great, 
		////size = 20
		//data_model.label_selection_smoothness = 500;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.number_random_points = data_model.target.V.rows() * .2;
		//data_model.ruled_width = data_model.target_max_dimension * 0.75;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.1;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		////Random:: cone_7-30_uniform2590
		////size = 15
		//data_model.label_selection_smoothness = 2000;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.number_random_points = data_model.target.V.rows() * .2;
		//data_model.ruled_width = data_model.target_max_dimension * 0.75;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.1;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		////Random:: face mask-aimshape
		////size = 40
		//data_model.label_selection_smoothness = 100;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.number_random_points = data_model.target.V.rows() * .1;
		//data_model.ruled_width = data_model.target_max_dimension * 0.3;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.035;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.1; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		//Random:: masonry
		//size = 30
		data_model.label_selection_smoothness = 250;
		data_model.use_geodesics_flat_detetction = false;
		data_model.geodesics_stopping.outlier_factor = 5.0;
		data_model.geodesics_stopping.window_size = 10.0;
		data_model.geodesics_stopping.winding_factor = 2.0;
		data_model.geodesics_directions = { 0 };
		data_model.number_random_points = data_model.target.V.rows() * .1;
		data_model.ruled_width = data_model.target_max_dimension * 0.3;
		data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.05;
		data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		data_model.optimization_settings.outlier_normal_threshold = 1.0;



		//------ PAPER RESULTS -------


		////Bunny *low*
		////...selected indices : 458, 292, 2129, 550, 1104, 985, 2485, 535, 157, 1265, 904, 307, 977, 2175, 1932
		////add ears: 1590, 795; add back feet: 1838, 1262; (add front feet: 1268, 1456)
		////data_model.ruled_vertex_indices = { 997, 296, 183, 55, 983 }; //testing normal deviation, pretty random, delete
		//data_model.geodesics_directions = {0};
		////data_model.ruled_vertex_indices = { 458, 292, 2129, 550, 1104, 985, 2485, 535, 157, 1265, 904, 307, 977, 2175, 1932, 1590, 795, 1838, 1262 };
		//data_model.ruled_vertex_indices = { 458, 292, 2129, 550, 1104, 985, 2485, 535, 157, 1265, 904, 307, 977, 2175, 1932, 1590, 795, 1838, 1262, 1268, 1456 };
		//data_model.label_selection_smoothness = 10;
		//data_model.use_geodesics_flat_detetction = true;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1;
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;
		
		////Bunny low again
		////...selected indices : 458, 1160, 630, 656, 1585, 2378, 775, 114, 2356, 1161, 119, 1714, 905, 144, 1812, 2429, 2081, 1723, 1375, 847, 2261, 1229, 1215, 2407
		////data_model.ruled_vertex_indices = { 458, 292, 2129, 550, 1104, 985, 2485, 535, 157, 1265, 904, 307, 977, 2175, 1932, 1590, 795, 1838, 1262, 1268, 1456 };
		////...selected indices : 384, 713, 2146, 1786, 1624, 1518, 735, 1242, 85, 161, 224, 658, 153, 990, 1325, 1864, 1451, 2221, 2274
		////add: tail 894, head 374, ears 795 1597
		//data_model.ruled_vertex_indices = { 384, 713, 2146, 1786, 1624, 1518, 735, 1242, 85, 161, 224, 658, 153, 990, 1325, 1864, 1451, 2221, 2274, 894, 374, 795, 1597 };
		//data_model.label_selection_smoothness = 300;
		//data_model.geodesics_directions = { 0 };
		//data_model.use_geodesics_flat_detetction = true;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.geodesics_stopping.window_size = 10;
		//data_model.ruled_width = data_model.target_max_dimension * 0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1;
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;
		


		////Bunny *high*
		////...selected indices : 458, 292, 2129, 550, 1104, 985, 2485, 535, 157, 1265, 904, 307, 977, 2175, 1932
		////add ears: 1590, 795; add back feet: 1838, 1262; (add front feet: 1268, 1456)
		//data_model.geodesics_directions = {0};
		//data_model.use_geodesics_flat_detetction = true;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.target_max_dimension = 45;
		//data_model.label_selection_smoothness = 350;
		//data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.025;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 0.5;
		//data_model.optimization_settings.outlier_normal_threshold = 1.6;


		
		////Fandisk
		////...selected indices : 426, 594, 2311, 2344, 590, 1254, 940, 2227, 2042, 1356, 372, 1268, 300, 1613, 2058, 744, 2046, 907, 808, 1128
		////data_model.ruled_vertex_indices = { 426, 594, 2311, 2344, 590, 1254, 940, 2227, 2042, 1356, 372, 1268, 300, 1613, 2058, 744, 2046, 907, 1128 }; //removed 808
		//data_model.ruled_vertex_indices = { 426, 594, 2164, 2344, 590, 1254, 940, 2227, 2042, 1356, 372, 1268, 300, 1613, 2058, 744, 2046, 907, 1128 }; //removed 808
		//data_model.use_geodesics_flat_detetction = true;
		//data_model.geodesics_stopping.outlier_factor = 6.0;
		//data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.02;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.001;
		//data_model.optimization_settings.outlier_normal_threshold = 0.7;
		////GlobalSettings::flat_angle_threshold = 0.002;


		//////Pavilion
		//////data_model.ruled_vertex_indices = { 1744 }; //pavilion
		//////data_model.ruled_width = data_model.target_max_dimension * 1.5;
		//////data_model.ruled_vertex_indices = { 1783 }; //pavilion
		//////data_model.ruled_vertex_indices = { 455 }; //pavilion, short legs
		//////data_model.ruled_vertex_indices = { 138, 544 }; //pavilion, short legs, fixing rulings
		//////data_model.ruled_vertex_indices = { 544 }; //pavilion, short legs, fixing rulings
		//////data_model.ruled_vertex_indices = { 647 }; 
		////new pavilion
		////data_model.ruled_vertex_indices = { 120 }; 
		////data_model.ruled_vertex_indices = { 83 }; 
		//data_model.ruled_vertex_indices = { 932 }; 
		//data_model.ruled_width = data_model.target_max_dimension * 1.25;
		//data_model.label_selection_smoothness = 50000;
		//data_model.geodesics_stopping.outlier_factor = 6.0;
		//data_model.optimization_settings.weight_bending_energy = 4.0;
		//data_model.optimization_settings.coverage_threshold = 1e3;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 2.5;
		//data_model.optimization_settings.outlier_normal_threshold = 4;//1.8;



		//// CUBE (works with making patch larger)
		//data_model.use_geodesics_flat_detetction = true;
		//data_model.geodesics_stopping.outlier_factor = 6.0;
		//data_model.ruled_width = data_model.target_max_dimension * 1.2;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 2.5;
		//data_model.optimization_settings.outlier_normal_threshold = 0.8;
		////GlobalSettings::flat_angle_threshold = 0.002;

		
		////Bumpy side--irregular 5 patches
		////...selected indices : 1113, 851, 2224, 1474
		//data_model.ruled_vertex_indices = { 1113, 851, 2224, 1474 };
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.035;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.25;
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		////Bumpy side--regular V2477 11 patches
		////...selected indices : 22, 1303, 2309, 338, 480, 1708, 1449, 1828, 1076, 464
		//data_model.ruled_vertex_indices = { 22, 1303, 2309, 338, 480, 1708, 1449, 1828, 1076, 464 };
		//data_model.target_max_dimension = 20;
		//data_model.target_max_dimension = 40; //for higher resolution
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.label_selection_smoothness = 350;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.005;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 0.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.6;
		
		
		////Bumpy side--regular V2477 3 patches
		//////...selected indices : 2008, 1097, 1218
		////...selected indices : 2447, 1655, 1374
		////data_model.ruled_vertex_indices = { 2008, 1097, 1218 };
		////data_model.ruled_vertex_indices = { 1320,1031,742 };
		//data_model.target_max_dimension = 20;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.label_selection_smoothness = 10000;
		//data_model.geodesics_stopping.outlier_factor = 8.0;
		////data_model.ruled_width = data_model.target_max_dimension*0.5;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.1;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;
		////these are the settings for the orthogonally going patches over the top
		////data_model.ruled_vertex_indices = { 737,776 }; //manually found
		//data_model.geodesics_directions = { 90 };
		//data_model.ruled_width = data_model.target_max_dimension*1.1;
		//data_model.geodesics_stopping.winding_factor = 1.5;


		////Bumpy side, recreate for error map
		//data_model.ruled_vertex_indices = { 377, 876, 749, 1022 };
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		//data_model.geodesics_stopping.absolute_threshold = 0.7;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.0; // 0.25;
		//data_model.optimization_settings.coverage_threshold = 0.8;
		//data_model.optimization_settings.weight_bending_energy = 2.5;
		//data_model.optimization_settings.weight_fitting_energy = 1.5;
		

		//LILIUM vary granularity
		//// -- 4-5 patches --
		////data_model.label_selection_smoothness = 3000;
		////data_model.ruled_vertex_indices = { 677, 584, 1963, 2214 }; //...selected indices : 677, 584, 1963, 2214
		//data_model.label_selection_smoothness = 2000;
		////...selected indices : 2069, 2427, 1220, 2871, 2158 -- GREAT!
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		////data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.ruled_width = data_model.target_max_dimension * 0.75;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.1;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;

		////lower
		////...selected indices : 1834, 507, 1383, 601, 572, 1235, 244, 1997, 288, 271
		////data_model.label_selection_smoothness = 500;
		//data_model.ruled_vertex_indices = { 1834, 507, 1383, 601, 572, 1235, 244, 1997, 288, 271 };
		//data_model.label_selection_smoothness = 5;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 3.0;
		////data_model.optimization_settings.weight_bending_energy = 3.0;
		////data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.ruled_width = data_model.target_max_dimension * 0.75;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.1; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;


		//data_model.label_selection_smoothness = 2000;
		//data_model.use_geodesics_flat_detetction = false;
		//data_model.geodesics_stopping.outlier_factor = 5.0;
		//data_model.geodesics_stopping.window_size = 10.0;
		//data_model.geodesics_stopping.winding_factor = 1.5;
		//data_model.ruled_width = data_model.target_max_dimension * 0.75;
		//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * 0.02;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.5; //used absolute 0.3
		//data_model.optimization_settings.outlier_normal_threshold = 1.0;
		////smoothness: 2000 ...selected indices : 2196, 454, 33, 2848
		////data_model.ruled_vertex_indices = { 2196, 454, 33, 2848 };
		////smoothness: 500 ...selected indices : 2651, 562, 2402, 1643, 2900, 191
		////data_model.ruled_vertex_indices = { 2651, 562, 2402, 1643, 2900, 191 };
		//////smoothness: 300 ...selected indices : 3035, 2838, 1785, 483, 1147, 2900, 2402, 191
		////smoothness: 250 ...selected indices : 1535, 3360, 1425, 80, 1208, 1400, 165, 1963, 663, 2691, 1508, 371, 2769
		//data_model.ruled_vertex_indices = { 1535, 3360, 1425, 80, 1208, 1400, 165, 1963, 663, 2691, 1508, 371, 2769 };
		////smoothness: 100 gives 16 geodesics, too much!
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;


		//data_model.optimization_settings.weight_bending_energy = 10.0;
		//smoothness: 750
		//...selected indices : 1828, 1383, 594, 1112, 111, 1977, 1539, 2958
		//...selected indices : 3157, 1383, 479, 111, 1977, 2896, 1112, 256, 2958
		//vary size



		////over constraining DOGs (use cover target constraints only)
		//data_model.ruled_vertex_indices = { 1055 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//data_model.optimization_settings.outlier_threshold = 10;
		//data_model.optimization_settings.outlier_normal_threshold = 1.57;
		//data_model.geodesics_stopping.winding_factor = 0.5;










		//data_model.ruled_vertex_indices = { 3073, 3139, 3119, 3252, 43, 232, 281, 1079, 958, 1896, 1866, 755, 1266, 1741, 1765, 1609, 1412, 2499, 2405, 2689, 2880, 195, 2544 }; //ear vertices from graph cut: 3167, 2967
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;



		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 0.33; //0.2; //~ edgelength * 0.33
		//data_model.optimization_settings.coverage_threshold = data_model.target.average_edge_length * 1.5; // try to avoid adding corners


		////dog model
		////1049, 1145, 1735, 1272, 1461, 1918, 1138, 1344, 18, 727, 463, 173, 574, 1397, 1510, 752, 1469 // from graph cut
		//	//+ 75 1350 133 1688 54 505 130 681
		//	//- 1510 1461 1344
		////data_model.ruled_vertex_indices = { 1049, 1145, 1735, 1272, 1918, 1138, 18, 727, 463, 173, 574, 1397, 752, 1469, 75, 1350, 133, 1688, 54, 505, 130, 681 }; 
		////data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//data_model.ruled_width = 15;
		//data_model.label_selection_smoothness = 200;

		//data_model.geodesics_stopping.winding_factor = 3.0;
		//data_model.geodesics_stopping.window_size = 8.0;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.0;
		//data_model.optimization_settings.coverage_threshold = 0.6;



		//////bunny edge length: 0.86...
		////data_model.ruled_vertex_indices = { 3073, 3139, 3119, 3252, 1840 }; 
		//data_model.ruled_vertex_indices = { 3073, 3139, 3119, 3252, 43, 232, 281, 1079, 958, 1896, 1866, 755, 1266, 1741, 1765, 1609, 1412, 2499, 2405, 2689, 2880, 195, 2544 }; //ear vertices from graph cut: 3167, 2967
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * 1.0;
		//data_model.optimization_settings.coverage_threshold = 0.8;
		//data_model.geodesics_directions = { 0 };



		////bumpy side 2477  (METHOD EXAMLE!) edge length: 0.516255
		//data_model.ruled_vertex_indices = { 876, 1022, 749, 377 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 1000;
		//data_model.optimization_settings.outlier_threshold = 0.25; //~ edgelength * 0.5
		////--> changed relax factor from 3 to 1!!
		//data_model.optimization_settings.coverage_threshold = 0.25; //~ edgelength * 1.0
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_stopping.window_size = 10;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;
		//data_model.geodesics_directions = { 0 };

		//data_model.is_pausing_on_adding_patch = true;
		//data_model.is_pausing_on_changing_constraint_strategy = true;


		////bumpy side 2477  (granularity 9 patches + 1 hole) edge length: 0.516255
		//data_model.ruled_vertex_indices = { 1932, 640, 135, 2273, 1628, 537, 1703, 1793, 626 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;
		//data_model.optimization_settings.outlier_threshold = 0.13; //~ edgelength * 0.25
		////data_model.optimization_settings.coverage_threshold = 0.13; //~ edgelength * 1.0
		//data_model.optimization_settings.coverage_threshold = 0.25; //~ edgelength * 1.0
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_stopping.window_size = 10;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;
		//data_model.geodesics_directions = { 0 };

		//data_model.is_pausing_on_adding_patch = true;
		//data_model.is_pausing_on_changing_constraint_strategy = true;



		// --> in loading above!!
		//// lilium  edge length: 0.41339
		//data_model.ruled_vertex_indices = { 670, 1496, 3038, 2739, 745 };
		////data_model.ruled_vertex_indices = { 1779, 562, 2453, 2365, 2739, 2858, 1379, 2158, 1562  };
		////data_model.ruled_vertex_indices = { 1779, 16, 514, 446, 1882, 3282, 2684, 2403, 802, 1409, 2369, 2672, 2739, 1914, 3037, 2514, 3091, 3145, 199, 614, 1442, 728, 710, 718, 3006, 2223, 1562, 1775, 320, 3025 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;
		//data_model.optimization_settings.outlier_threshold = 0.2; //~ edgelength * 0.5
		//data_model.optimization_settings.coverage_threshold = 0.41; //~ edgelength * 1.0
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_directions = { 0 };
		////data_model.geodesics_stopping.window_size = 10;
		////data_model.geodesics_stopping.bounds_min_absolute = 0.05;

		////data_model.is_pausing_on_adding_patch = true;
		////data_model.is_pausing_on_changing_constraint_strategy = true;







		////bunny head, fixing unfolding
		//data_model.ruled_vertex_indices = { 500, 601, 251 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//data_model.is_pausing_on_adding_patch = true;
		//data_model.is_pausing_on_changing_constraint_strategy = true;

		////bunny, fixing unfolding
		//data_model.ruled_vertex_indices = { 2904, 3042, 2641 }; //116, 322, 807, 982, 1424, 1487, 2028, 2532, 1952, 2301, 2988, 2904, 3042
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;
		//data_model.is_pausing_on_adding_patch = true;
		//data_model.is_pausing_on_changing_constraint_strategy = true;



		////TODO make video! (dropbox/2020-01-18_19-33-16)
		//// bumpy-side-irregular  edge length: 0.599774
		//data_model.ruled_vertex_indices = { 803, 222, 964, 1290, 2303 };
		//data_model.selected_ruled_vertex_indices = data_model.ruled_vertex_indices;

		//data_model.ruled_width = 10;
		//data_model.label_selection_smoothness = 500;
		//data_model.optimization_settings.outlier_threshold = 0.2; //~ edgelength * 0.33
		////data_model.optimization_settings.coverage_threshold = 0.6; //~ edgelength * 1.0
		//data_model.optimization_settings.coverage_threshold = 0.9; // try to avoid adding corners
		//data_model.optimization_settings.weight_bending_energy = 4.0;

		//data_model.geodesics_directions = { 0 };
		//data_model.geodesics_stopping.window_size = 10;
		//data_model.geodesics_stopping.bounds_min_absolute = 0.05;

		////data_model.is_pausing_on_adding_patch = true;
		////data_model.is_pausing_on_changing_constraint_strategy = true;




		//data_model.ruled_vertex_indices = { 3017, 2739, 1801, 3177, 2007, 2464, 31, 2813, 2567, 2293, 2632, 179, 3247, 2025 };


		//bumpy
		//data_model.ruled_vertex_indices = { 3926, 34, 3576, 3604, 1702, 4296, 178, 1796, 1815, 3017 };


		////double cone
		//data_model.ruled_vertex_indices = { 895, 335, 947, 922 };
		//data_model.ruled_width = 12.0;
			   

		////bunny
		//data_model.ruled_vertex_indices.push_back(2902);
		//data_model.label_selection_smoothness = 500;
		//data_model.ruled_width = 6.0;






		if (is_starting_optimization)
			optimization.initialize();
	}
	else //start empty, load via UI
	{
		data_model.models_folder = "";

		GlobalSettings::max_average_edge = 1e6; //prevent upsampling
		view.view_model.target_view_settings.standard_fill_color = Eigen::RowVector3d(179.0 / 255.0, 176.0 / 255.0, 174.0 / 255.0); 
	}

	viewer_register_callbacks();
	data_model.viewer.launch();
}
