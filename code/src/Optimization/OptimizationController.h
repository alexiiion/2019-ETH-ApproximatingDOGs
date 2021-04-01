#pragma once

#include "DataModel.h"
#include "GeodesicsController.h"

#include "RuledGeodesicsController.h"
//#include "FeatureGeodesicsController.h"
//#include "FeatureGeodesics.h"

#include "Result/DevelopableResultController.h"

#include "CompositeObjective.h"
#include "QuadraticConstraintsSumObjective.h"
#include "BarycentricPositionalConstraints.h"
#include "LengthRegularizerConstraints.h"
#include "SimplifiedBendingObjective.h"
#include "IsometryObjective.h"
#include "DOGConstraints.h"
//#include "NewtonKKT.h"
#include "EqSQP.h"




class OptimizationController
{
public:
	OptimizationController(DataModel& data_model) : data_model(data_model), settings(data_model.optimization_settings)
	{
		geodesics_controller = new RuledGeodesicsController(data_model);
		result_controller = new DevelopableResultController(data_model);
		//result_controller = new DevelopableResultController(data_model.developable_model, data_model.target, data_model.patches, data_model.is_optimizing, data_model.viewer);
		//result_controller = new DevelopableResultController(data_model.target, data_model.patches, data_model.is_optimizing, data_model.viewer);
		//result_controller = new DevelopableResultController(data_model);
	};
	~OptimizationController() 
	{ 
		delete geodesics_controller; 

		delete solver;
		delete dog_constraints;

		delete composite_objective;
		delete fitting_constraints;
		delete fitting_objective;
		delete bending_energy;
		delete iso_energy;
	};
	

	RuledGeodesicsController* geodesics_controller;
	//GeodesicsController* geodesics_controller;
	DevelopableResultController* result_controller;

	EqSQP* solver;
	DogConstraints* dog_constraints;

	CompositeObjective* composite_objective;
	BarycentricPositionalConstraints* fitting_constraints;
	QuadraticConstraintsSumObjective* fitting_objective;
	LengthRegularizerConstraints* length_constraints;
	QuadraticConstraintsSumObjective* length_objective;
	SimplifiedBendingObjective* bending_energy;
	IsometryObjective* iso_energy;


	void initialize();
	void run_optimiztaion();

	void optimize(Patch& patch);
	void update_position_constraints();
	void update_patch_rendering(bool do_update_all = false);
	void try_add_patch();

	void reset_mesh_optimiztaion();


private:
	DataModel& data_model;
	OptimizationSettings& settings;
	OptimizationSettings initial_settings;


	bool is_initialized = false;
	bool do_update_objectives = true;
	int current_patch;

	double patch_init_time = 0;

	bool should_adapt_weights = false;
	int valid_DOG_iterations;
	int max_valid_DOG_iterations;
	int max_iterations;

	bool is_first = true;
	igl::Timer timer;
	double init_time;

	
	//double initial_fitting_weight;
	//double initial_bending_weight;
	//double initial_isometry_weight;

	//solver settings
	double solver_target_developability_error = 1e-3; //infeasability_epsilon
	double solver_max_developability_error = 0.01; //UNUSED infeasability_filter
	int solver_max_iterations = 15;
	double solver_merit_p = 10;


	void run_DOG_optimiztaion();
	void run_mesh_optimiztaion();

	void add_patch();

	void update_objectives(Patch& patch);
	void update_target_positions(Patch& patch);
	void update_coverage();

	void update_fast_patch_rendering();
	//void update_precise_patch_rendering();

	//void merge_patches();
	void process_optimization_result(Patch& patch);
	void process_result_adaptive(Patch& patch);
	void process_result_static(Patch& patch);

	//void update_render_mesh(Patch& patch, const std::vector<Eigen::MatrixXd>& intersection_curves);

	//bool find_intersection_curve(const WrapperMesh& wrapper, const WrapperMesh& other_wrapper, Eigen::MatrixXd& out_intersection_curve);
	//void find_min_error_point(const Patch& patch, int& out_vertex_index);
};

