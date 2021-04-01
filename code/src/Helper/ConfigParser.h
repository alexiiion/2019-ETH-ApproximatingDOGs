#pragma once

#include <string>
class DataModel;

struct StartConfig
{
public:
	bool is_loading_scene = false;
	bool is_loading_target = false;
	bool is_starting_optimization = false;
	std::string scene_file;

	double number_random_points_factor = -1;
	double ruled_width_factor = -1;
	double coverage_threshold_factor = -1;
	double outlier_threshold_factor = -1;

	void read_startup_config(const std::string& config_filepath, DataModel& data_model);
	void set_mesh_dependent_settings(DataModel& data_model);

private:
	void read_model_config(const std::string& config_filepath, DataModel& data_model);

};
