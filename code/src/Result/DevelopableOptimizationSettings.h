#pragma once


struct DevelopableOptimizationSettings
{
	float label_selection_smoothness = 20000;

	//energy weights
	float weight_stiching_objective = 1.0;
	float weight_developable_objective = 0.001;

	float weight_proximity_objective = 0.001;
	float weight_boundary_objective = 0.001;

	float weight_smooth_objective = 0.5;

	float weight_proximity_inner_objective = 1.0;
	float weight_boundary_seam_objective = 5.0;

	//thresholds
	float angle_defect_threshold = 0.01;

	//debug mode
	bool use_debug_mode = false;
};