#pragma once

struct OptimizationSettings
{

	//energy weights
	float weight_isometry_energy = 1.0;
	float weight_regularizing_energy = 0.0;
	float weight_bending_energy = 2.5;
	float weight_fitting_energy = 1.5;

	//double fitting_objective;
	//double developable_energy;

	//thresholds
	float convergence_fitting_threshold = 0.01; //0.002;
	float convergence_developablity_threshold = 0.0001;

	float outlier_threshold = 0.2;
	float outlier_normal_threshold = 1.0;
	float coverage_threshold = 0.1;
	//float geodesics_coverage_threshold = 0.5;

	bool should_converge = false;

	//old energies
	//float laplacian_multiplier = 1e-6;
	//float weight_laplacian_energy = 0;
	//float weight_bilaplacian_energy = 0;
	//float weight_smoothness_energy = 1;
};
