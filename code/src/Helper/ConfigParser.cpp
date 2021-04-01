#include "ConfigParser.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <iostream>
#include <fstream>

#include "DataModel.h"
#include "Utils.h"
#include "Logger.h"

using namespace std;


std::string config_model_path;

void string_to_vector(const string& input, vector<int>& out_list, const char deliminator = ' ')
{
	std::stringstream iss(input);

	int number;
	out_list.clear();
	while (iss >> number)
		out_list.push_back(number);
}

std::string parse_string(std::string content)
{
	std::string variable = content;
	//write_log(3) << content << ": " << variable << endl;

	return variable;
}

bool parse_bool(std::string content)
{
	int variable = atoi(content.c_str());
	//write_log(3) << content << ": " << boolalpha << variable << endl;

	return variable;
}

int parse_int(std::string content)
{
	int variable = atoi(content.c_str());
	//write_log(3) << content << ": " << variable << endl;

	return variable;
}

float parse_float(std::string content)
{
	float variable = atof(content.c_str());
	//write_log(3) << content << ": " << variable << endl;

	return variable;
}

std::vector<int> parse_int_vector(std::string content, const char deliminator = ' ')
{
	std::vector<int> variable;
	string_to_vector(content, variable, deliminator);
	//write_log(3) << content << ": " << list_to_string(variable) << endl;

	return variable;
}

void StartConfig::read_startup_config(const std::string& config_filepath, DataModel& data_model)
{
	// Open the File
	std::ifstream in(config_filepath.c_str());

	// Check if object is valid
	if (!in)
	{
		std::cerr << "Cannot open the STARTUP config file : " << config_filepath << std::endl;
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
		else if (qualifier == "config_model")
		{
			config_model_path = parse_string(content);
			break;
		}

		//if (qualifier == "stop_reading")
		//{
		//	break;
		//}
		//else if (qualifier == "start_empty")
		//{
		//	is_loading_scene = false;
		//	is_loading_target = false;
		//
		//	write_log(3) << "starting empty " << endl;
		//	write_log(3) << "-- stopping reading config here." << endl;
		//	break;
		//}
		//else if (qualifier == "load_scene")
		//{
		//	is_loading_scene = true;
		//	scene_file = parse_string(content);
		//	write_log(3) << "-- stopping reading config here." << endl;
		//	break;
		//}
		//else if (qualifier == "target")
		//{
		//	is_loading_target = true;
		//	data_model.target_filename = parse_string(content);
		//}
		//else if (qualifier == "models_folder")
		//{
		//	data_model.models_folder = parse_string(content);
		//}
		//else if (qualifier == "target_max_dimension")
		//{
		//	data_model.target_max_dimension = parse_int(content);
		//}
		//else if (qualifier == "max_average_edge")
		//{
		//	GlobalSettings::max_average_edge = parse_int(content);
		//}

	}

	//Close The File
	in.close();

	
	
	if (!config_model_path.empty())
		read_model_config(get_folder_path(config_filepath) + config_model_path, data_model);
	
	write_log(2) << endl << "DONE reading config settings. " << endl << endl;
}

void StartConfig::read_model_config(const std::string& config_filepath, DataModel& data_model)
{
	// Open the File
	std::ifstream in(config_filepath.c_str());

	// Check if object is valid
	if (!in)
	{
		std::cerr << "Cannot open the MODEL config file : " << config_filepath << std::endl;
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

		if (qualifier == "stop_reading")
		{
			break;
		}
		else if (qualifier == "start_empty")
		{
			is_loading_scene = false;
			is_loading_target = false;

			write_log(3) << "starting empty " << endl;
			write_log(3) << "-- stopping reading config here." << endl;
			break;
		}
		else if (qualifier == "load_scene")
		{
			is_loading_scene = true;
			scene_file = parse_string(content);

			write_log(3) << "loading scene " << endl;
			write_log(3) << "-- stopping reading config here." << endl;
			break;
		}
		else if (qualifier == "target")
		{
			is_loading_target = true;
			data_model.target_filename = parse_string(content);
		}
		else if (qualifier == "initialize_optimization")
			is_starting_optimization = parse_bool(content);
		else if (qualifier == "models_folder")
			data_model.models_folder = parse_string(content);
		else if (qualifier == "target_max_dimension")
			data_model.target_max_dimension = parse_int(content);
		else if (qualifier == "max_average_edge")
			GlobalSettings::max_average_edge = parse_int(content);

		else if (qualifier == "label_selection_smoothness")
			data_model.label_selection_smoothness = parse_int(content);
		else if (qualifier == "ruled_vertex_indices")
			data_model.ruled_vertex_indices = parse_int_vector(content, ',');
		else if (qualifier == "number_random_points")
			data_model.number_random_points = parse_int(content);
		else if (qualifier == "number_random_points_factor")
			number_random_points_factor = parse_float(content);
			//data_model.number_random_points = data_model.target.V.rows() * parse_float(content);
		else if (qualifier == "ruled_width")
			data_model.ruled_width = parse_int(content);
		else if (qualifier == "ruled_width_factor")
			ruled_width_factor = parse_float(content);
			//data_model.ruled_width = data_model.target_dimensions.maxCoeff() * parse_int(content);
		
		else if (qualifier == "geodesics_directions")
			data_model.geodesics_directions = parse_int_vector(content, ',');
		else if (qualifier == "use_geodesics_flat_detetction")
			data_model.use_geodesics_flat_detetction = parse_int(content);
		else if (qualifier == "geodesics_stopping.outlier_factor")
			data_model.geodesics_stopping.outlier_factor = parse_float(content);
		else if (qualifier == "geodesics_stopping.window_size")
			data_model.geodesics_stopping.window_size = parse_float(content);
		else if (qualifier == "geodesics_stopping.winding_factor")
			data_model.geodesics_stopping.winding_factor = parse_float(content);
		
		else if (qualifier == "optimization_settings.coverage_threshold")
			data_model.optimization_settings.coverage_threshold = parse_float(content);
		else if (qualifier == "optimization_settings.coverage_threshold_factor")
			coverage_threshold_factor = parse_float(content);
			//data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * parse_float(content);
		else if (qualifier == "optimization_settings.outlier_threshold")
			data_model.optimization_settings.outlier_threshold = parse_float(content);
		else if (qualifier == "optimization_settings.outlier_threshold_factor")
			outlier_threshold_factor = parse_float(content);
			//data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * parse_float(content);
		else if (qualifier == "optimization_settings.outlier_normal_threshold")
			data_model.optimization_settings.outlier_normal_threshold = parse_float(content);
		else if (qualifier == "optimization_settings.weight_bending_energy")
			data_model.optimization_settings.weight_bending_energy = parse_float(content);
		else if (qualifier == "optimization_settings.convergence_developablity_threshold")
			data_model.optimization_settings.convergence_developablity_threshold = parse_float(content);
	}

	//Close The File
	in.close();
}

void StartConfig::set_mesh_dependent_settings(DataModel& data_model)
{
	data_model.update_flat_threshold();

	if(!is_equal(number_random_points_factor, -1))
		data_model.number_random_points = data_model.target.V.rows() * number_random_points_factor;
	if(!is_equal(ruled_width_factor, -1))
		data_model.ruled_width = data_model.target_dimensions.maxCoeff() * ruled_width_factor;
	if(!is_equal(coverage_threshold_factor, -1))
		data_model.optimization_settings.coverage_threshold = data_model.target_max_dimension * coverage_threshold_factor;
	if(!is_equal(outlier_threshold_factor, -1))
		data_model.optimization_settings.outlier_threshold = data_model.target.average_edge_length * outlier_threshold_factor;
}

/*
void read_config(DataModel& data_model)
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
*/