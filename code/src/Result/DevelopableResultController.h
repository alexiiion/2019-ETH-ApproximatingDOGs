#pragma once

#include <vector>
#include <Eigen/Core>

#include <igl/Timer.h>

#include "MeshModel.h"
#include "GaussianCurvature.h"
#include "DevelopableModel.h"
#include "DataModel.h"
#include "DevelopableOptimizationSettings.h"

#include "CompositeObjective.h"
#include "QuadraticConstraintsSumObjective.h"
#include "AngleDeficitConstraints.h"
#include "LaplacianObjective.h"
#include "PointPairConstraints.h"
#include "PositionalConstraints.h"

#include "Newton.h"


enum OptimizationState
{
	Projection, 
	Developability,
	Finished
};

class DevelopableResultController
{
public:
	//DevelopableResultController(DataModel& data_model) : data_model(data_model)
	//{ /* empty on purpose */ };
	//DevelopableResultController(Mesh& target, std::vector<Patch*>& patches, bool& is_optimizing) : patches(patches), is_optimizing(is_optimizing), viewer()
	//{
	//	developable_model.non_developable = target;
	//};

	DevelopableResultController(DataModel& data_model)
		: developable_model(data_model.developable_model), target(data_model.target), patches(data_model.patches), result_coverage(data_model.result_coverage), is_optimizing(data_model.is_optimizing), viewer(data_model.viewer)
	{
		developable_model.non_developable = target;
		original_settings = developable_model.optimization_settings;
		dimension = 3;
	};
	//TODO remove!!
	//DevelopableResultController(Developable& developable_model, Mesh& target, std::vector<Patch*>& patches, bool& is_optimizing, igl::opengl::glfw::Viewer& viewer) 
	//	: developable_model(developable_model), target(target), patches(patches), is_optimizing(is_optimizing), viewer(viewer)
	//{
	//	developable_model.non_developable = target;
	//	original_settings = developable_model.optimization_settings;
	//	dimension = 3;
	//};

	~DevelopableResultController()
	{
		reset();
	};

	//should be just one method, add detection if target_developable has changed (subdivide, new assignment, etc), otherwise just run shape optimization
	void reset();
	void initialize(bool use_debug_mode = false);
	void run_optimiztaion();
	void finish_optimiztaion();


private:
	//DataModel& data_model;
	Mesh& target;
	std::vector<Patch*>& patches;
	bool& is_optimizing;
	bool is_initialized = false;

	//TODO remove!!
	igl::opengl::glfw::Viewer& viewer; 
	std::vector<int> developable_parts_view_indices;
	int concatenated_view_index = -1;

	int dimension;

	bool is_first = true;
	igl::Timer timer;
	double init_time;


	Developable& developable_model;
	//GaussCurvature* gaussian_curvature;
	Coverage*& result_coverage;


	DevelopableOptimizationSettings original_settings;
	OptimizationState optimization_state;
	int iteration = 0;
	int max_iterations = 30;
	int inter_iteration = 0;
	int max_inter_iterations = 20;


	Eigen::SparseMatrix<double> L_initial;
	Eigen::SparseMatrix<double> M_initial;
	Eigen::VectorXd inner_indices_mass; //masses for inner indices only (for AngleDeficitObjective)
	Eigen::MatrixXd all_vertices_projected;


	bool are_objectives_initialized = false;
	CompositeObjective* composite_objective;

	PointPairConstraints* stiching_constraints;
	QuadraticConstraintsSumObjective* stiching_objective;

	AngleDeficitConstraints* developable_constraints;
	QuadraticConstraintsSumObjective* developable_objective;

	PositionalConstraints* proximity_constraints;
	QuadraticConstraintsSumObjective* proximity_objective;

	PositionalConstraints* boundary_constraints;
	QuadraticConstraintsSumObjective* boundary_objective;


	PositionalConstraints* inner_proximity_constraints;
	QuadraticConstraintsSumObjective* inner_proximity_objective;

	PositionalConstraints* boundary_seam_constraints;
	QuadraticConstraintsSumObjective* boundary_seam_objective;

	LaplacianObjective* smooth_objective;

	//Newton solver;



	void update_concatenated_mesh(); //update after each optimization step

	void initialize_optimization();
	void reset_objectives();
	void update_objectives(const Eigen::VectorXd& x0);
	void process_optimization_result(const Eigen::VectorXd& x0, const Eigen::VectorXd& x);

	void initialize_weights();
	void update_projection_objectives(const Eigen::VectorXd& x0);
	void update_developability_objectives(const Eigen::VectorXd & x0);

	std::vector<int> get_optimized_face_labeling(bool use_debug_mode);
	void find_seam_indices(const std::vector<int>& face_patch_assignment);
	Eigen::MatrixXi get_seam_pairs(std::vector<std::vector<int>>& global_patch_correspondance, std::vector<int>& seam_patch_lookup);

	void create_massmatrix();
	void create_laplacian();

	void find_inner_indices();
	void find_boundary_indices();

	void create_meshes(const std::vector<int>& face_patch_assignment);
	void create_developables_parts_meshes(std::vector<std::vector<int>>& label_vertex_assignment, std::vector<std::vector<int>>& label_face_assignment);

	void project_developable_parts();

	void build_correspondance(const std::vector<int>& face_patch_assignment);
	std::vector<int> get_assigned_labels(std::vector<std::vector<int>>& label_vertex_assignment);
	std::vector<int> get_seam_patch_lookup(std::vector<std::vector<int>>& global_patch_correspondance);
	std::vector<int> get_cumulative_offsets();
	int get_closest_index_on(Eigen::RowVectorXd vertex, const int& patch_index);

};