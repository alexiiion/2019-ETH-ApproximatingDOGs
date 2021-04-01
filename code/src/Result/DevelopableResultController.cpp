#include "DevelopableResultController.h"

#include <algorithm>

#include <igl/point_mesh_squared_distance.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/upsample.h>

#include <igl/facet_components.h>
#include <igl/remove_duplicate_vertices.h>

#include <igl/slim.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/MappingEnergyType.h>
#include <igl/flipped_triangles.h>


#include "PatchModel.h"
#include "GraphCut.h"
#include "MeshController.h"

#include "Utils.h"
#include "ViewUtils.h"
#include "Colors.h"
#include "Logger.h"

using namespace std;

void uniform_laplacian(const Eigen::MatrixXi F, Eigen::SparseMatrix<double>& L);
void stitch_mesh(Developable& developable_model);
void flatten_patches(Developable& developable_model);

void add_cwise(int offset, std::vector<std::vector<int>>& list);
void add_cwise(int offset, std::vector<int>& list);

void flatten_indices(const Eigen::MatrixXi& matrix, const int element_offset, const int dimension, Eigen::MatrixXi& out_flat);
void flatten_indices(const Eigen::VectorXi& input, const int element_offset, const int dimension, Eigen::VectorXi& out_flat);

void flatten_indices(const std::vector<int>& input, const int element_offset, const int dimension, std::vector<int>& out_flat);
void flatten_indices(const std::vector<std::vector<int>>& input, const int element_offset, const int dimension, std::vector<std::vector<int>>& out_flat);

//bool has_only_finite(const std::vector<int>& list)
//{
//	int non_finite = std::count_if(std::begin(list), std::end(list), 
//					 [](const auto& value) { return !std::isfinite(value); });
//
//	return non_finite == 0;
//}


void DevelopableResultController::run_optimiztaion()
{
	if (!is_initialized)
	{
		//reset();
		initialize(developable_model.optimization_settings.use_debug_mode);
		initialize_weights();
		is_optimizing = false;
		//DataModel::log_level_optimization = 4; 
	}

	if (!is_initialized)
		return;
	if (!is_optimizing)
		return;

	if (is_first)
	{
		init_time = timer.getElapsedTime();
		is_first = false;
		return;
	}

	if (optimization_state == OptimizationState::Finished)
	{
		double t = timer.getElapsedTime();
		write_log(2) << endl << endl << "Finished optimizing DEVELOPABLE -- ELAPSED TIME: " << t - init_time << endl << endl;

		is_optimizing = false;
		return;
	}

	Eigen::VectorXd x0;
	mat2_to_vec(developable_model.concatenated.V, x0);

	update_objectives(x0);

	// solve
	Eigen::VectorXd x;
	Newton solver;

	solver.solve(x0, *composite_objective, x);
	process_optimization_result(x0, x);

	// update mesh
	vec_to_mat2(x, developable_model.concatenated.V);
	update_concatenated_mesh();
	result_coverage->compute(developable_model.concatenated);

	//TODO remove, make output global and set in objective
	auto angle_defect = developable_constraints->Vals(x);
	developable_model.output_angle_defect = angle_defect.norm();
	developable_model.output_stiching = stiching_constraints->Vals(x).norm();

	//update K
	developable_model.gaussian_curvature.update(developable_model.inner_indices_mask);

	write_log(4) << "angle defect:  norm = " << developable_model.output_angle_defect << ", sum = " << angle_defect.sum();
	write_log(4) << " --> K:: avg = " << developable_model.gaussian_curvature.average() << " (" << developable_model.gaussian_curvature.get().minCoeff() << " - " << developable_model.gaussian_curvature.get().maxCoeff() << ") " << endl;
}

void DevelopableResultController::reset()
{
	if (!is_initialized)
		return;


	developable_model.reset();
	developable_model.non_developable = target;

	reset_objectives();

	is_initialized = false;

	iteration = 0;
	inter_iteration = 0;

	double smooth = developable_model.optimization_settings.label_selection_smoothness;
	developable_model.optimization_settings = original_settings;
	developable_model.optimization_settings.label_selection_smoothness = smooth;



	//TODO move to view!
	for (int i = 0; i < developable_parts_view_indices.size(); i++)
		viewer.data_list[developable_parts_view_indices[i]].clear();

	developable_parts_view_indices.clear();

	if (concatenated_view_index >= 0 && concatenated_view_index < viewer.data_list.size())
	{
		viewer.data_list[concatenated_view_index].clear();
		concatenated_view_index = -1;
	}
}

void DevelopableResultController::reset_objectives()
{
	if (!are_objectives_initialized)
		return;

	delete composite_objective;
	delete stiching_constraints;
	delete stiching_objective;
	delete developable_constraints;
	delete developable_objective;

	delete proximity_constraints;
	delete proximity_objective;
	delete boundary_constraints;
	delete boundary_objective;

	delete inner_proximity_constraints;
	delete inner_proximity_objective;
	delete boundary_seam_constraints;
	delete boundary_seam_objective;

	delete smooth_objective;

	are_objectives_initialized = false;
}

void DevelopableResultController::finish_optimiztaion()
{
	optimization_state = OptimizationState::Finished;

	stitch_mesh(developable_model);
	flatten_patches(developable_model);
}


void DevelopableResultController::initialize_optimization()
{
	Eigen::VectorXd x0; //TODO check: this is the only values that used to change
	mat2_to_vec(developable_model.concatenated.V, x0);


	//flat boundary indices
	Eigen::VectorXi boundary_indices_flat;
	flatten_indices(developable_model.boundary_indices, developable_model.number_concatenated_vertices, dimension, boundary_indices_flat);

	//flat boundary vertices
	Eigen::VectorXd boundary_vertices_flat;
	Eigen::MatrixXd bv;
	index_to_value(developable_model.boundary_indices, developable_model.concatenated.V, bv);
	mat2_to_vec(bv, boundary_vertices_flat);


	//flat seams
	Eigen::MatrixXi seam_pairs_flat;
	flatten_indices(developable_model.seam_pairs, developable_model.number_concatenated_vertices, dimension, seam_pairs_flat);


	////flat inner vertices
	//std::vector<int> inner_indices_flat;
	//std::vector<std::vector<int>> inner_indices_adjacency_flat;
	//flatten_indices(developable_model.inner_indices, developable_model.number_concatenated_vertices, dimension, inner_indices_flat);
	//flatten_indices(developable_model.inner_indices_adjacency, developable_model.number_concatenated_vertices, dimension, inner_indices_adjacency_flat);



	composite_objective = new CompositeObjective();

	////int f = std::count_if(std::begin(developable_model.inner_indices), std::end(developable_model.inner_indices), [](const auto& value) { return std::isnan(value); });
	////write_log(0) << "finite? inner_indices: " << f << endl;
	//
	////write_log(0) << "finite? inner_indices: " << has_only_finite(developable_model.inner_indices) << endl;
	//write_log(0) << "finite? inner_indices_mass: " << (inner_indices_mass.hasNaN() == false) << endl;
	//write_log(0) << "inner_indices_mass: " << linebreak << inner_indices_mass << linebreak << endl;
	//
	//write_log(0) << "inner_indices: " << linebreak << list_to_string(developable_model.inner_indices) << linebreak << endl;
	//
	//write_log(0) << "inner_indices_adjacency: " << endl;
	//for(int i = 0; i < developable_model.inner_indices_adjacency.size(); i++)
	//	write_log(0) << "  [" << list_to_string(developable_model.inner_indices_adjacency[i]) << "]" << endl;


	developable_constraints = new AngleDeficitConstraints(developable_model.inner_indices, developable_model.inner_indices_adjacency, inner_indices_mass);
	developable_objective = new QuadraticConstraintsSumObjective(*developable_constraints, x0);
	composite_objective->add_objective(developable_objective, developable_model.optimization_settings.weight_developable_objective);
	//cout << "developable_objective.grad(x0).norm() " << developable_objective.grad(x0).norm() << endl;

	stiching_constraints = new PointPairConstraints(seam_pairs_flat);
	stiching_objective = new QuadraticConstraintsSumObjective(*stiching_constraints, x0);
	composite_objective->add_objective(stiching_objective, developable_model.optimization_settings.weight_stiching_objective);

	/*
	Eigen::VectorXd all_vertices_flat;
	mat2_to_vec(developable_model.concatenated.V, all_vertices_flat);
	int concatenated_rows = all_vertices_flat.rows();
	Eigen::VectorXi all_indices(concatenated_rows);
	all_indices.setLinSpaced(concatenated_rows, 0, concatenated_rows - 1);

	proximity_constraints = new PositionalConstraints(all_indices, all_vertices_flat);
	proximity_objective = new QuadraticConstraintsSumObjective(*proximity_constraints, x0);
	composite_objective->add_objective(proximity_objective, developable_model.optimization_settings.weight_proximity_objective);
	*/

	Eigen::VectorXd all_vertices_flat;
	mat2_to_vec(all_vertices_projected, all_vertices_flat);

	int concatenated_rows = all_vertices_flat.rows();
	Eigen::VectorXi all_indices(concatenated_rows);
	all_indices.setLinSpaced(concatenated_rows, 0, concatenated_rows - 1);

	proximity_constraints = new PositionalConstraints(all_indices, all_vertices_flat);
	proximity_objective = new QuadraticConstraintsSumObjective(*proximity_constraints, x0);
	composite_objective->add_objective(proximity_objective, developable_model.optimization_settings.weight_proximity_objective);

	boundary_constraints = new PositionalConstraints(boundary_indices_flat, boundary_vertices_flat);
	boundary_objective = new QuadraticConstraintsSumObjective(*boundary_constraints, x0);
	composite_objective->add_objective(boundary_objective, developable_model.optimization_settings.weight_boundary_objective);

	smooth_objective = new LaplacianObjective(L_initial*-1);
	composite_objective->add_objective(smooth_objective, developable_model.optimization_settings.weight_smooth_objective);
}

void DevelopableResultController::process_optimization_result(const Eigen::VectorXd& x0, const Eigen::VectorXd& x)
{
	if (optimization_state == OptimizationState::Finished)
		return;

	if (optimization_state == OptimizationState::Developability)
	{
		if ((developable_model.gaussian_curvature.max() <= 0.005 && developable_model.output_stiching < 1e-3) || iteration >= max_iterations)
		{
			finish_optimiztaion();
			return;
		}
	}


	if (x.hasNaN())
	{
		write_log(1) << endl << "ERROR at optimization: result has NAN!" << endl << endl;
		is_optimizing = false;
		return;
	}

	double stiching_progress = stiching_constraints->Vals(x).norm();
	double proximity_progress = proximity_constraints->Vals(x).norm();

	write_log(4) << "stitching norm = " << stiching_progress << ", proximity norm = " << proximity_progress << endl;

	const double epsilon = 1e-4;
	double progress = (x - x0).squaredNorm();
	
	//const double epsilon = 3;
	//double progress = proximity_progress;

	inter_iteration++;

	if (progress < epsilon || inter_iteration >= max_inter_iterations)
	{
		iteration++;
		inter_iteration = 0;

		if (optimization_state == OptimizationState::Projection)
		{
			project_developable_parts();

			optimization_state = OptimizationState::Developability;

			developable_model.optimization_settings.weight_developable_objective = 5.0;
			developable_model.optimization_settings.weight_stiching_objective = 1.0;
			developable_model.optimization_settings.weight_proximity_objective = 1.0;
			developable_model.optimization_settings.weight_boundary_objective = 1.0;
			//developable_model.optimization_settings.weight_smooth_objective = 1.0;

			developable_model.optimization_settings.weight_proximity_inner_objective = 0.0;
			developable_model.optimization_settings.weight_boundary_seam_objective = 0.0;

			iteration = 1;
		}
	}
	
	write_log(0) << "after solve: progress = " << progress << ", inter_iteration = " << inter_iteration << ", iteration = " << iteration << endl;

	if (optimization_state == OptimizationState::Projection)
		update_projection_objectives(x0);
	else if (optimization_state == OptimizationState::Developability)
		update_developability_objectives(x0);
}


void DevelopableResultController::initialize_weights()
{
	developable_model.optimization_settings.weight_stiching_objective = 0.001;
	developable_model.optimization_settings.weight_developable_objective = 0.001;

	//developable_model.optimization_settings.weight_proximity_objective = 1.0;
	//developable_model.optimization_settings.weight_boundary_objective = 1.0;

	developable_model.optimization_settings.weight_smooth_objective = 0.5;

	developable_model.optimization_settings.weight_proximity_inner_objective = 1.0;
	developable_model.optimization_settings.weight_boundary_seam_objective = 5.0;
}

void DevelopableResultController::update_projection_objectives(const Eigen::VectorXd& x0)
{
	return;

	////just keep initial values because it runs only for 1 iteration
	//developable_model.optimization_settings.weight_stiching_objective = 5.0;
	//developable_model.optimization_settings.weight_proximity_objective = 2.0;
	//developable_model.optimization_settings.weight_boundary_objective = 1.0;
	//developable_model.optimization_settings.weight_smooth_objective = 1.0;
	//
	//developable_model.optimization_settings.weight_proximity_inner_objective = 0.0;
	//developable_model.optimization_settings.weight_boundary_seam_objective = 0.0;

	double mu = 1.005;
	developable_model.optimization_settings.weight_stiching_objective *= mu;
	developable_model.optimization_settings.weight_proximity_objective *= mu;
	//developable_model.optimization_settings.weight_boundary_objective = 1.0;

	project_developable_parts();

	return;


	//developable_model.optimization_settings.weight_stiching_objective = 0.001;
	//developable_model.optimization_settings.weight_developable_objective = 0.001;
	//
	//developable_model.optimization_settings.weight_proximity_objective = 0.001;
	//developable_model.optimization_settings.weight_boundary_objective = 0.001;
	//
	//developable_model.optimization_settings.weight_smooth_objective = 0.5;
	//
	//developable_model.optimization_settings.weight_proximity_inner_objective = 1.0;
	//developable_model.optimization_settings.weight_boundary_seam_objective = 5.0;
}

void DevelopableResultController::update_developability_objectives(const Eigen::VectorXd& x0)
{
	if (iteration < 2)
		return;

	float next_w_stitch = pow(iteration, 2);
 	developable_model.optimization_settings.weight_stiching_objective = max(next_w_stitch, developable_model.optimization_settings.weight_stiching_objective);
	developable_model.optimization_settings.weight_developable_objective = pow(iteration, 3);

	//developable_model.optimization_settings.weight_proximity_objective = 0.001;
	//developable_model.optimization_settings.weight_boundary_objective = 0.001;
	//
	//developable_model.optimization_settings.weight_smooth_objective = 0.5;
	//
	//developable_model.optimization_settings.weight_proximity_inner_objective = 1.0;
	//developable_model.optimization_settings.weight_boundary_seam_objective = 5.0;
}


void DevelopableResultController::update_objectives(const Eigen::VectorXd& x0)
{
	reset_objectives();


	//if (optimization_state == OptimizationState::Developability)
	//	project_developable_parts();


	//all vertices
	Eigen::VectorXd all_vertices_flat;
	mat2_to_vec(all_vertices_projected, all_vertices_flat);

	//all indices
	int concatenated_rows = all_vertices_flat.rows();
	Eigen::VectorXi all_indices(concatenated_rows);
	all_indices.setLinSpaced(concatenated_rows, 0, concatenated_rows - 1);


	//inner vertices
	Eigen::MatrixXd inner_vertices_projected(developable_model.inner_indices.size(), dimension);
	for (int i = 0; i < developable_model.inner_indices.size(); i++)
		inner_vertices_projected.row(i) = all_vertices_projected.row(developable_model.inner_indices[i]);

	Eigen::VectorXd inner_vertices_flat;
	mat2_to_vec(inner_vertices_projected, inner_vertices_flat);

	//inner indices
	Eigen::VectorXi inner_indices = Eigen::Map<Eigen::VectorXi>(developable_model.inner_indices.data(), developable_model.inner_indices.size());
	Eigen::VectorXi inner_indices_flat;
	flatten_indices(inner_indices, developable_model.number_concatenated_vertices, dimension, inner_indices_flat);


	//flat boundary indices
	Eigen::VectorXi boundary_indices_flat;
	flatten_indices(developable_model.boundary_indices, developable_model.number_concatenated_vertices, dimension, boundary_indices_flat);

	//flat boundary vertices (original, undeformed ones but indices refer to concatenated)
	Eigen::VectorXd boundary_vertices_flat;
	mat2_to_vec(developable_model.boundary_vertices_original, boundary_vertices_flat);

	////flat boundary vertices (from deformed concatenated mesh)
	//Eigen::VectorXd boundary_vertices_flat;
	//Eigen::MatrixXd boundary_vertices;
	//index_to_value(developable_model.boundary_indices, all_vertices_projected, boundary_vertices);
	////index_to_value(developable_model.boundary_indices, developable_model.concatenated.V, bv);
	//mat2_to_vec(boundary_vertices, boundary_vertices_flat);


	//flat seams
	Eigen::MatrixXi seam_pairs_flat;
	flatten_indices(developable_model.seam_pairs, developable_model.number_concatenated_vertices, dimension, seam_pairs_flat);


	//boundary + seams
	std::vector<int> bi(developable_model.boundary_indices.data(), developable_model.boundary_indices.data() + developable_model.boundary_indices.size());
	int number_seams = developable_model.seam_pairs.rows();
	Eigen::VectorXi s1 = developable_model.seam_pairs.block(0, 0, number_seams, 1);
	Eigen::VectorXi s2 = developable_model.seam_pairs.block(0, 1, number_seams, 1);
	std::vector<int> si1(s1.data(), s1.data() + s1.size());
	std::vector<int> si2(s2.data(), s2.data() + s2.size());

	bi.insert(bi.end(), si1.begin(), si1.end());
	bi.insert(bi.end(), si2.begin(), si2.end());

	sort(bi.begin(), bi.end());
	bi = remove_duplicates(bi);

	//boundary + seam indices
	Eigen::VectorXi boundary_seam_indices = Eigen::Map<Eigen::VectorXi>(bi.data(), bi.size());
	Eigen::VectorXi boundary_seam_indices_flat;
	flatten_indices(boundary_seam_indices, developable_model.number_concatenated_vertices, dimension, boundary_seam_indices_flat);

	//boundary + seam vertices
	Eigen::MatrixXd boundary_seam_vertices;
	index_to_value(bi, all_vertices_projected, boundary_seam_vertices);
	
	Eigen::VectorXd boundary_seam_vertices_flat;
	mat2_to_vec(boundary_seam_vertices, boundary_seam_vertices_flat);





	composite_objective = new CompositeObjective();


	developable_constraints = new AngleDeficitConstraints(developable_model.inner_indices, developable_model.inner_indices_adjacency, inner_indices_mass);
	developable_objective = new QuadraticConstraintsSumObjective(*developable_constraints, x0);
	composite_objective->add_objective(developable_objective, developable_model.optimization_settings.weight_developable_objective);
	//cout << "developable_objective.grad(x0).norm() " << developable_objective.grad(x0).norm() << endl;

	stiching_constraints = new PointPairConstraints(seam_pairs_flat);
	stiching_objective = new QuadraticConstraintsSumObjective(*stiching_constraints, x0);
	composite_objective->add_objective(stiching_objective, developable_model.optimization_settings.weight_stiching_objective);

	proximity_constraints = new PositionalConstraints(all_indices, all_vertices_flat);
	proximity_objective = new QuadraticConstraintsSumObjective(*proximity_constraints, x0);
	composite_objective->add_objective(proximity_objective, developable_model.optimization_settings.weight_proximity_objective);

	boundary_constraints = new PositionalConstraints(boundary_indices_flat, boundary_vertices_flat);
	boundary_objective = new QuadraticConstraintsSumObjective(*boundary_constraints, x0);
	composite_objective->add_objective(boundary_objective, developable_model.optimization_settings.weight_boundary_objective);


	inner_proximity_constraints = new PositionalConstraints(inner_indices_flat, inner_vertices_flat);
	inner_proximity_objective = new QuadraticConstraintsSumObjective(*inner_proximity_constraints, x0);
	composite_objective->add_objective(inner_proximity_objective, developable_model.optimization_settings.weight_proximity_inner_objective);

	if (boundary_seam_indices_flat.hasNaN())
		write_log(0) << "boundary_seam_indices_flat.hasNaN()" << endl;
	if (boundary_seam_vertices_flat.hasNaN())
		write_log(0) << "boundary_seam_vertices_flat.hasNaN()" << endl;

	boundary_seam_constraints = new PositionalConstraints(boundary_seam_indices_flat, boundary_seam_vertices_flat);
	boundary_seam_objective = new QuadraticConstraintsSumObjective(*boundary_seam_constraints, x0);
	composite_objective->add_objective(boundary_seam_objective, developable_model.optimization_settings.weight_boundary_seam_objective);


	smooth_objective = new LaplacianObjective(L_initial*-1);
	composite_objective->add_objective(smooth_objective, developable_model.optimization_settings.weight_smooth_objective);

	are_objectives_initialized = true;



	//TODO delete!
	concatenated_view_index = viewer.selected_data_index;
	viewer.data().add_points(boundary_seam_vertices, Eigen::RowVector3d(0.001, 0.001, 0.001));
	//viewer.data().add_edges(all_vertices_projected, developable_model.concatenated.V, Eigen::RowVector3d(0.001, 0.001, 0.001));
}


void DevelopableResultController::initialize(bool use_debug_mode)
{
	if (is_initialized)
		return;

	//igl::upsample(developable_model.non_developable.V, developable_model.non_developable.F, developable_model.non_developable.V, developable_model.non_developable.F, 1);
	//igl::loop(data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.V, data_model.target_developable.F, 2);


	dimension = target.V.cols();
	std::vector<int> face_patch_assignment = get_optimized_face_labeling(use_debug_mode);
	write_log(0) << "face_patch_assignment: " << linebreak << list_to_string(face_patch_assignment, "", true) << endl;

	create_meshes(face_patch_assignment);
	developable_model.number_concatenated_vertices = developable_model.cumulative_offset_vertices.back() + developable_model.developable_parts.back().V.rows();
	
	build_correspondance(face_patch_assignment);
	find_seam_indices(face_patch_assignment);

	////project developable part to their assigned patch
	//if (!use_debug_mode)
		project_developable_parts(); 

	developable_model.concatenated = meshhelper::concatenate_meshes(developable_model.developable_parts);
	igl::per_vertex_normals(developable_model.concatenated.V, developable_model.concatenated.F, developable_model.concatenated.normals_vertices);

	////ONLY DEBUG (contains very small number eg < 1e-14)
	//Eigen::VectorXd concatenated_F_areas = meshhelper::compute_face_areas(developable_model.concatenated.V, developable_model.concatenated.F);
	//write_log(0) << "concatenated_F_areas: (sum = " << concatenated_F_areas.sum() << ")" << linebreak << concatenated_F_areas << endl;
	//if (concatenated_F_areas.hasNaN())
	//	write_log(0) << "concatenated_F_areas.hasNaN!" << endl;
	////end ONLY DEBUG


	find_inner_indices();
	find_boundary_indices();

	create_laplacian();
	create_massmatrix();

	developable_model.gaussian_curvature.update(developable_model.inner_indices_mask);
	//gaussian_curvature = new GaussCurvature(developable_model.concatenated);


	if (use_debug_mode)
	{
		is_optimizing = false;
		DataModel::log_level_optimization = 4; 
	}


	is_initialized = true;
	is_first = true;



	//TODO should be done in view
	for (int i = 0; i < developable_model.assigned_labels.size(); i++)
	{
		Mesh* mesh = &developable_model.developable_parts[i];
		//meshhelper::add_mesh(data_model.viewer, data_model.developable_parts[i].V, data_model.developable_parts[i].F, false, false);
		meshhelper::add_mesh(viewer, mesh->V, mesh->F, false, false);
		developable_parts_view_indices.push_back(viewer.selected_data_index);
		
		viewer.data().compute_normals();

		Eigen::Vector4f line_color(Colors::GRAY_MID(0), Colors::GRAY_MID(1), Colors::GRAY_MID(2), 1.0f);
		viewer.data().line_color = line_color;
		viewer.data().line_width = 0.5;
		viewer.data().point_size = 5.0;

		viewer.data().set_colors(get_color(i, developable_model.number_patches));

		//int color_step = 80.0 / developable_model.assigned_labels.size();
		//color_step = max(color_step, 10);
		//double gray_value = (255.0 - color_step*i) / 255.0;
		//Eigen::RowVector3d rgb(gray_value, gray_value, gray_value);
		//viewer.data().set_colors(rgb);
	}



	//DEBUG stuff
	write_log(4) << "assigned " << developable_model.assigned_labels.size() << " labels: " << endl; log_list(4, developable_model.assigned_labels, "", false);


	//for quick debug view
	viewer.append_mesh();
	viewer.data().point_size = 10;
	concatenated_view_index = viewer.selected_data_index;

	//viewer.data().add_points(all_vertices_projected, Eigen::RowVector3d(0.001, 0.001, 0.001));
	//viewer.data().add_edges(all_vertices_projected, developable_model.concatenated.V, Eigen::RowVector3d(0.001, 0.001, 0.001));

	//Eigen::MatrixXd inner_vertices_points;
	//index_to_value(developable_model.inner_indices, developable_model.concatenated.V, inner_vertices_points);
	//viewer.data().add_points(inner_vertices_points, Eigen::RowVector3d(0.001, 0.001, 0.001));

	//Eigen::MatrixXd boundary_vertices_points;
	//index_to_value(boundary_indices, concatenated.V, boundary_vertices_points);
	//data_model.viewer.data().add_points(boundary_vertices_points, Eigen::RowVector3d(0.001, 0.001, 1));

	//Eigen::MatrixXd boundary_vertices_points_from_flat;
	//vec_to_mat2(boundary_vertices_flat, boundary_vertices_points_from_flat);fla
	//data_model.viewer.data().add_points(boundary_vertices_points_from_flat, Eigen::RowVector3d(0.001, 1, 1));

	optimization_state = OptimizationState::Projection;
	is_initialized = true;
}

void DevelopableResultController::create_laplacian()
{
	igl::cotmatrix(developable_model.concatenated.V, developable_model.concatenated.F, L_initial);

	//gets NaN likely due to bad coverage and tiny (almost zero area) triangles
	//bool nan = L_initial.toDense().hasNaN();
	double sum = L_initial.sum();
	if (!std::isnan(sum))
		return;

	write_log(2) << linebreak << "  WARN  cotan Laplacian sum = " << sum << " --> use uniform Laplacian instead." << endl;
	
	//use uniform laplacian instead
	uniform_laplacian(developable_model.concatenated.F, L_initial);

	sum = L_initial.sum();
	if (std::isnan(sum))
	{
		write_log(1) << linebreak << "ERROR  L:: sum = " << sum << linebreak << endl;
		exit(ERROR);
	}

	write_log(2) << "  INFO  uniform Laplacian is valid. Sum = " << sum << linebreak << endl;
}

void uniform_laplacian(const Eigen::MatrixXi F, Eigen::SparseMatrix<double>& L)
{
	Eigen::SparseMatrix<double> A;
	igl::adjacency_matrix(F, A);
	
	// sum each row
	Eigen::SparseVector<double> Asum;
	igl::sum(A, 1, Asum);
	
	// Convert row sums into diagonal of sparse matrix
	Eigen::SparseMatrix<double> Adiag;
	igl::diag(Asum, Adiag);
	
	L = A - Adiag;
}

void DevelopableResultController::create_massmatrix()
{
	igl::massmatrix(developable_model.concatenated.V, developable_model.concatenated.F, igl::MASSMATRIX_TYPE_VORONOI, M_initial);
	
	//bool nan = M_initial.toDense().hasNaN();
	double sum = M_initial.sum();
	if (std::isnan(sum))
		write_log(2) << linebreak << "  ERROR  Mass matrix sum = " << sum << " --> looks like bad coverage! (can re-run labelling (graph cut), but this is not a fix!)" << endl;


	Eigen::VectorXd M_diagonal = M_initial.diagonal();
	//write_log(0) << "M_diagonal:: has NaN? " << M_diagonal.hasNaN() << endl;

	int count = 0;
	inner_indices_mass.resize(developable_model.inner_indices.size());
	for (int inner_vi : developable_model.inner_indices)
	{
		inner_indices_mass(count) = M_diagonal(inner_vi);
		count++;
	}
	//write_log(0) << "inner_indices_mass:: has NaN? " << inner_indices_mass.hasNaN() << endl;


	//write_log(0) << "M_initial: " << endl << M_initial << endl << endl;
	//write_log(0) << "M_initial.diagonal() size = " << MV.rows() <<  " : " << endl << MV << endl << endl;
	//write_log(0) << "inner_vertices_mass size = " << inner_indices_mass.rows() <<  " : " << endl << inner_indices_mass << endl << endl;
}

void DevelopableResultController::update_concatenated_mesh()
{
	if (developable_model.concatenated.V.rows() < 1)
		return; 

	for (int i = 0; i < developable_model.cumulative_offset_vertices.size(); i++)
	{
		int start = developable_model.cumulative_offset_vertices[i];
		int end = i >= developable_model.cumulative_offset_vertices.size() - 1 ? developable_model.concatenated.V.rows() : developable_model.cumulative_offset_vertices[i + 1];
		int size = end - start;

		write_log(6) << "[" << i << "]: start = " << start << ", end = " << end << ", size = " << size << endl;

		//Mesh* mesh = &data_model.developable_parts[i];
		//mesh->V = concatenated.V.block(start, 0, size, 3);
		developable_model.developable_parts[i].V = developable_model.concatenated.V.block(start, 0, size, 3);

		//TODO in view
		int view_index = developable_parts_view_indices[i];
		//data_model.viewer.data_list[view_index].set_vertices(mesh->V);
		viewer.data_list[view_index].set_vertices(developable_model.developable_parts[i].V);
		viewer.data_list[view_index].compute_normals();
	}
}

std::vector<int> DevelopableResultController::get_optimized_face_labeling(bool use_debug_mode)
{
	//if (is_initialized)
	//	return;
	//if (data_model.patches.size() < 1)
	//	return;


	Eigen::MatrixXd F_midpoints = meshhelper::compute_face_midpoints(developable_model.non_developable.V, developable_model.non_developable.F);
	Eigen::VectorXd F_areas = meshhelper::compute_face_areas(developable_model.non_developable.V, developable_model.non_developable.F);

	//normalize area
	double max_area = F_areas.maxCoeff();
	F_areas /= max_area;


	//DEBUG for testing angle deficit constraints
	std::vector<Mesh> meshes;
	if (use_debug_mode)
		meshes.push_back(developable_model.non_developable);
	else
		meshes = meshhelper::get_wrapper_meshes(patches);

	for (Mesh& mesh : meshes)
		igl::upsample(mesh.V, mesh.F, 2);

	developable_model.number_patches = meshes.size();

	meshhelper::project_to_meshes(F_midpoints, F_areas, meshes, developable_model.patch_face_distances, developable_model.patch_closest_points);
	Eigen::MatrixXi F_edges = meshhelper::face_edges(developable_model.non_developable.F);

	std::vector<double> resulting_distances;
	std::vector<int> face_patch_assignment = graphcut::compute_graph_cut_face_labeling(F_midpoints, F_edges, developable_model.patch_face_distances, developable_model.optimization_settings.label_selection_smoothness, resulting_distances);
	developable_model.face_patch_assignment = face_patch_assignment;

	return face_patch_assignment;
}

void DevelopableResultController::create_meshes(const std::vector<int>& face_patch_assignment)
{
	vector<vector<int>> label_face_assignment(developable_model.number_patches);
	vector<vector<int>> label_vertex_assignment(developable_model.number_patches);

	//get vertices and faces per label
	for (int fi = 0; fi < face_patch_assignment.size(); fi++)
	{
		int label = face_patch_assignment[fi];
		label_face_assignment[label].push_back(fi);

		Eigen::RowVectorXi v = developable_model.non_developable.F.row(fi);
		vector<int> vertices(v.data(), v.data() + v.size());
		label_vertex_assignment[label].insert(label_vertex_assignment[label].end(), vertices.begin(), vertices.end());
	}

	developable_model.assigned_labels = get_assigned_labels(label_vertex_assignment);
	create_developables_parts_meshes(label_vertex_assignment, label_face_assignment);
	developable_model.cumulative_offset_vertices = get_cumulative_offsets();
}

void DevelopableResultController::build_correspondance(const std::vector<int>& face_patch_assignment)
{
	const int loglevel = 5;

	developable_model.global_patch_correspondance.clear();
	developable_model.global_patch_correspondance.resize(developable_model.non_developable.V.rows());

	for (int fi = 0; fi < face_patch_assignment.size(); fi++)
	{
		int label = face_patch_assignment[fi];
		Eigen::RowVectorXi v = developable_model.non_developable.F.row(fi);
		vector<int> vertices(v.data(), v.data() + v.size());

		for (int vi : vertices)
			developable_model.global_patch_correspondance[vi].push_back(label);
	}

	write_log(loglevel+1) << endl << "RAW: " << endl;
	for (int i = 0; i < developable_model.global_patch_correspondance.size(); i++)
		log_list(loglevel, developable_model.global_patch_correspondance[i], "global_patch_correspondance[" + to_string(i) + "]", false);
}

void DevelopableResultController::find_seam_indices(const std::vector<int>& face_patch_assignment)
{
	if (developable_model.number_patches < 2)
		return;
	if (developable_model.assigned_labels.size() < 2)
		return;

	const int loglevel = 5;
	
	vector<int> seam_patch_lookup = get_seam_patch_lookup(developable_model.global_patch_correspondance);
	developable_model.seam_pairs = get_seam_pairs(developable_model.global_patch_correspondance, seam_patch_lookup);

	//flatten_indices(developable_model.seam_pairs, developable_model.number_concatenated_vertices, developable_model.non_developable.V.cols(), seam_pairs_flat);
	//write_log(loglevel) << "final: seam_pairs_flat: " << endl << seam_pairs_flat << endl << endl;
}

void DevelopableResultController::find_inner_indices()
{
	developable_model.inner_indices.clear();
	developable_model.inner_indices_adjacency.clear();

	for (int i = 0; i < developable_model.developable_parts.size(); i++)
	{
		Mesh* mesh = &developable_model.developable_parts[i];

		std::vector<int> indices(mesh->V.rows());
		std::iota(indices.begin(), indices.end(), 0);

		vector<vector<int>> adjacency_VV;
		igl::adjacency_list(mesh->F, adjacency_VV, true);

		std::vector<int> boundary = meshhelper::boundary_loop_flat(mesh->F);
		std::sort(boundary.begin(), boundary.end());


		for (int j = 0; j < boundary.size(); j++)
		{
			int bi = boundary[j];
			indices.erase(indices.begin() + bi - j);
			adjacency_VV.erase(adjacency_VV.begin() + bi - j);
		}

		int offset = developable_model.cumulative_offset_vertices[i];
		add_cwise(offset, indices);
		add_cwise(offset, adjacency_VV);

		developable_model.inner_indices.insert(developable_model.inner_indices.end(), indices.begin(), indices.end());
		developable_model.inner_indices_adjacency.insert(developable_model.inner_indices_adjacency.end(), adjacency_VV.begin(), adjacency_VV.end());
	}

	//flatten_indices(inner_indices, number_concatenated_vertices, dimension, inner_indices_flat);
	//flatten_indices(inner_indices_adjacency, number_concatenated_vertices, dimension, inner_indices_adjacency_flat);

	developable_model.inner_indices_mask.resize(developable_model.number_concatenated_vertices);
	developable_model.inner_indices_mask.setZero();
	for (int i = 0; i < developable_model.inner_indices.size(); i++)
		developable_model.inner_indices_mask(developable_model.inner_indices[i]) = 1;
}

void DevelopableResultController::find_boundary_indices()
{
	const int loglevel = 6;

	vector<int> global_flat_boundary = meshhelper::boundary_loop_flat(developable_model.non_developable.F);
	std::vector<int> boundary_concatenated;

	for (int bi : global_flat_boundary)
	{
		vector<int> patch_indices = developable_model.global_patch_correspondance[bi];
		write_log(loglevel) << bi << ": patch_indices: "; log_list(loglevel, patch_indices, "", false);

		Eigen::RowVectorXd vertex = developable_model.non_developable.V.row(bi);

		for (int pi : patch_indices)
		{
			int index = get_closest_index_on(vertex, pi);
			boundary_concatenated.push_back(index);
		}
	}
	write_log(loglevel) << "#V=" << developable_model.number_concatenated_vertices << ", boundary: " << list_to_string(boundary_concatenated, "", true) << endl;


	//store indices
	developable_model.boundary_indices = Eigen::Map<Eigen::VectorXi>(boundary_concatenated.data(), boundary_concatenated.size());

	//store vertices
	index_to_value(developable_model.boundary_indices, developable_model.concatenated.V, developable_model.boundary_vertices_original);
}

std::vector<int> DevelopableResultController::get_seam_patch_lookup(std::vector<std::vector<int>>& global_patch_correspondance)
{
	const int loglevel = 5;

	std::vector<int> seam_patch_lookup;
	for (int i = 0; i < global_patch_correspondance.size(); i++)
	{
		global_patch_correspondance[i] = remove_duplicates(global_patch_correspondance[i]);
		log_list(loglevel, global_patch_correspondance[i], "global_patch_correspondance[" + to_string(i) + "]", false);

		if (global_patch_correspondance[i].size() > 1)
		{
			seam_patch_lookup.push_back(i);
			write_log(loglevel) << "   --> store seam pair" << endl;
		}
	}
	log_list(loglevel, seam_patch_lookup, "seam_patch_lookup: ", false);

	return seam_patch_lookup;
}

std::vector<int> DevelopableResultController::get_assigned_labels(std::vector<std::vector<int>>& label_vertex_assignment)
{
	std::vector<int> assigned_labels;

	for (int i = 0; i < developable_model.number_patches; i++)
	{
		if(label_vertex_assignment[i].size() <= 0)
			continue;

		assigned_labels.push_back(i);
	}

	return assigned_labels;
}

void DevelopableResultController::create_developables_parts_meshes(std::vector<std::vector<int>>& label_vertex_assignment, std::vector<std::vector<int>>& label_face_assignment)
{
	const int loglevel = 5;

	developable_model.developable_parts.clear();

	for (int i = 0; i < developable_model.number_patches; i++)
	{
		vector<int> vertices = label_vertex_assignment[i];
		if (vertices.size() <= 0)
			continue;

		label_vertex_assignment[i] = remove_duplicates(vertices);

		Mesh mesh = meshhelper::sub_mesh(label_vertex_assignment[i], label_face_assignment[i], developable_model.non_developable.V, developable_model.non_developable.F);
		developable_model.developable_parts.push_back(mesh);
	}
}

Eigen::MatrixXi DevelopableResultController::get_seam_pairs(std::vector<std::vector<int>>& global_patch_correspondance, std::vector<int>& seam_patch_lookup)
{
	const int loglevel = 5;

	Eigen::MatrixXi seam_pairs(developable_model.non_developable.V.rows(), 2);
	seam_pairs.setConstant(-1);
	int number_seam_pairs = 0;

	//vector<vector<int>> global_local_correspondance(data_model.target_developable.V.rows()); //#V x max 2  

	for (int i = 0; i < seam_patch_lookup.size(); i++)
	{
		int vi = seam_patch_lookup[i];
		vector<int> patch_indices = global_patch_correspondance[vi];

		write_log(loglevel) << vi << ": patch_indices: "; log_list(loglevel, patch_indices, "", false);

		Eigen::RowVectorXd vertex = developable_model.non_developable.V.row(vi);

		int end = patch_indices.size();
		if(patch_indices.size() > 2)
			end = patch_indices.size() + 1;

		for (int j = 1; j < end; j++)
		{
			int pi_b = patch_indices[j - 1];
			int concatenated_index_b = get_closest_index_on(vertex, pi_b);
			//int mesh_index_b = index_of(pi_b, developable_model.assigned_labels);
			//int local_index_b = get_closest_vertex(developable_model.developable_parts[mesh_index_b].V, vertex);
			//int concatenated_index_b = local_index_b + developable_model.cumulative_offset_vertices[mesh_index_b];
			////int concatenated_index_b = local_index_b + cumulative_offset_vertices[pi_b];
			////global_local_correspondance[vi].push_back(local_index_b);
			//write_log(loglevel) << "   BACK:  mesh_index = " << mesh_index_b << ", local_index = " << local_index_b << ", concatenated_index = " << concatenated_index_b << endl;

			int pi_f = patch_indices[j%patch_indices.size()];
			int concatenated_index_f = get_closest_index_on(vertex, pi_f);
			//int mesh_index_f = index_of(pi_f, developable_model.assigned_labels);
			//int local_index_f = get_closest_vertex(developable_model.developable_parts[mesh_index_f].V, vertex);
			//int concatenated_index_f = local_index_f + developable_model.cumulative_offset_vertices[mesh_index_f];
			////int concatenated_index_f = local_index_f + cumulative_offset_vertices[pi_f];
			////global_local_correspondance[vi].push_back(local_index_f);
			//write_log(loglevel) << "   FRONT: mesh_index = " << mesh_index_f << ", local_index = " << local_index_f << ", concatenated_index = " << concatenated_index_f << endl;

			if(concatenated_index_b < 0 || concatenated_index_f < 0)
				write_log(loglevel) << "      Error" << endl;

			seam_pairs(number_seam_pairs, 0) = concatenated_index_b;
			seam_pairs(number_seam_pairs, 1) = concatenated_index_f;
			number_seam_pairs++;
		}
	}

	write_log(loglevel) << "final: concatenated_seam_pairs: " << endl << seam_pairs << endl << endl;
	seam_pairs.conservativeResize(number_seam_pairs, Eigen::NoChange);
	write_log(loglevel) << "resized: concatenated_seam_pairs: " << endl << seam_pairs << endl << endl;

	return seam_pairs;
}

std::vector<int> DevelopableResultController::get_cumulative_offsets()
{
	std::vector<int> cumulative_offset_vertices;
	cumulative_offset_vertices.push_back(0);

	for (int i = 1; i < developable_model.developable_parts.size(); i++)
	{
		int previous_offset = cumulative_offset_vertices[i - 1];
		int current_offset = developable_model.developable_parts[i - 1].V.rows();
		cumulative_offset_vertices.push_back(previous_offset + current_offset);
	}

	return cumulative_offset_vertices;
}

void DevelopableResultController::project_developable_parts()
{
	//for (int i = 0; i < developable_parts_view_indices.size(); i++)
	//	viewer.data_list[developable_parts_view_indices[i]].clear();

	//developable_parts_view_indices.clear();


	all_vertices_projected.resize(0, dimension);
	//all_vertices_projected.resize(developable_model.number_concatenated_vertices, dimension);

	for (int i = 0; i < developable_model.assigned_labels.size(); i++)
	{
		int label = developable_model.assigned_labels[i];
		Mesh* mesh = &developable_model.developable_parts[i];

		//TODO store points from previous computation (for graph cut)
		WrapperMesh& wrapper = patches[label]->wrapper;
		Eigen::VectorXi face_indices;
		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd closest_points;
		igl::point_mesh_squared_distance(mesh->V, wrapper.V, wrapper.F, distances_squared, face_indices, closest_points);

		append_matrix(closest_points, all_vertices_projected);
		//all_vertices_projected << closest_points;
		//mesh->V = closest_points;

	}
	//write_log(0) << "all_vertices_projected: " << linebreak << all_vertices_projected << linebreak << endl;


	/*
	//project developable part to it's assigned DOG
	for (int i = 0; i < developable_model.assigned_labels.size(); i++)
	{
		int label = developable_model.assigned_labels[i];
		Mesh* mesh = &developable_model.developable_parts[i];

		//TODO store points from previous computation (for graph cut)
		WrapperMesh& wrapper = patches[label]->wrapper;
		Eigen::VectorXi face_indices;
		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd closest_points;
		igl::point_mesh_squared_distance(mesh->V, wrapper.V, wrapper.F, distances_squared, face_indices, closest_points);

		mesh->V = closest_points;
		
		////TODO check!!! --> doesn't work because somehow indices dont correspond anymore. investigate why??
		//mesh->V = patch_closest_points[i];



		////TODO should be done in view
		////meshhelper::add_mesh(data_model.viewer, data_model.developable_parts[i].V, data_model.developable_parts[i].F, false, false);
		//meshhelper::add_mesh(data_model.viewer, mesh->V, mesh->F, false, false);
		//data_model.developable_parts_view_indices.push_back(data_model.viewer.selected_data_index);

		////data_model.viewer.data().show_lines = false;
		////data_model.viewer.data().show_vertid = true;
		//data_model.viewer.data().set_colors(get_color(i, number_patches));
	}
	*/
}

void add_cwise(int offset, std::vector<std::vector<int>>& list)
{
	if (offset == 0)
		return;

	for (int i = 0; i < list.size(); i++)
		add_cwise(offset, list[i]);
}

void add_cwise(int offset, std::vector<int>& list)
{
	if (offset == 0)
		return;

	for (int i = 0; i < list.size(); i++)
		list[i] += offset;
}

//move to Utils
void flatten_indices(const Eigen::MatrixXi& matrix, const int element_offset, const int dimension, Eigen::MatrixXi& out_flat)
{
	const int size = matrix.rows();
	const int number_vertices = element_offset;

	const int num_cols = matrix.cols();
	out_flat.resize(size * dimension, num_cols);

	Eigen::RowVectorXi ones(num_cols);
	ones.setOnes();

	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < dimension; j++) 
		{
			Eigen::RowVectorXi vertex = matrix.row(i);
			out_flat.row(i + size * j) = vertex + j * element_offset*ones;

			write_log(6) << "flatten_indices[" << i << "][" << j << "](" << i + size * j << ") = " << out_flat.row(i + size * j) << endl;

			if (i + size * j >= size * dimension)
				write_log(1) << "ERROR at flatten_indices" << endl;
		}
	}
}

void flatten_indices(const Eigen::VectorXi& input, const int element_offset, const int dimension, Eigen::VectorXi& out_flat)
{
	const int size = input.rows();
	const int number_vertices = element_offset;

	out_flat.resize(size * dimension);
	out_flat.setZero();

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			out_flat(i + j*size) = input(i) + j * number_vertices;
			write_log(6) << "flatten_indices[" << i << "][" << j << "](" << i + size * j << ") = " << out_flat.row(i + size * j) << endl;

			if (i + size * j >= size * dimension)
				write_log(1) << "ERROR at flatten_indices" << endl;
		}
	}
}

void flatten_indices(const std::vector<int>& input, const int element_offset, const int dimension, std::vector<int>& out_flat)
{
	const int size = input.size();
	const int number_vertices = element_offset;

	out_flat.clear();
	out_flat.resize(size * dimension);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			out_flat[i + j*size] = input[i] + j * number_vertices;
			write_log(6) << "flatten_indices[" << i << "][" << j << "](" << i + size * j << ") = " << out_flat[i + size * j] << endl;

			if (i + size * j >= size * dimension)
				write_log(1) << "ERROR at flatten_indices" << endl;
		}
	}
}

void flatten_indices(const std::vector<std::vector<int>>& input, const int element_offset, const int dimension, std::vector<std::vector<int>>& out_flat)
{
	const int size = input.size();
	const int number_vertices = element_offset;

	out_flat.clear();
	out_flat.resize(size * dimension);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			std::vector<int> offset_list = input[i];
			add_cwise(j * number_vertices, offset_list);

			out_flat[i + j * size] = offset_list;
		}
	}
}

int DevelopableResultController::get_closest_index_on(Eigen::RowVectorXd vertex, const int& patch_index)
{
	int mesh_index = index_of(patch_index, developable_model.assigned_labels);
	int local_index = get_closest_vertex(developable_model.developable_parts[mesh_index].V, vertex);
	int concatenated_index = local_index + developable_model.cumulative_offset_vertices[mesh_index];

	return concatenated_index;
}

/*
void DevelopableResultController::boundary_stuff()
{
	const int loglevel = 5;

	//Eigen::VectorXi concatenated_boundary_indices;
	//Eigen::MatrixXd concatenated_boundary_vertices;
	get_boundary_on_concatenated_developable(global_patch_correspondance, concatenated_boundary_indices, concatenated_boundary_vertices);

	mat2_to_vec(concatenated_boundary_vertices, flattened_boundary_vertices);

	flatten_indices(concatenated_boundary_indices, number_concatenated_vertices, data_model.target_developable.V.cols(), flattend_boundary_indices);
	write_log(loglevel) << "final: concatenated_boundary_indices: " << endl << concatenated_boundary_indices << endl << endl;
	write_log(loglevel) << "final: flattend_boundary_indices: " << endl << flattend_boundary_indices << endl << endl;



	//setup inner vertrices for angle deficit constraints

	inner_indices.clear();
	inner_indices_adjacency.clear();

	for (int i = 0; i < data_model.developable_parts.size(); i++)
	{
		Mesh* mesh = &data_model.developable_parts[i];

		std::vector<int> indices(mesh->V.rows());
		std::iota(indices.begin(), indices.end(), 0);

		vector<vector<int>> adjacency_VV;
		igl::adjacency_list(mesh->F, adjacency_VV, true);

		std::vector<int> boundary = meshhelper::boundary_loop_flat(mesh->F);
		std::sort(boundary.begin(), boundary.end());


		//write_log(loglevel) << endl << "BEFORE erasing boundary: " << endl;
		//log_list(loglevel, boundary, "boundary: ", false);
		//log_list(loglevel, indices, "indices: ", false);
		//write_log(loglevel) << "adjacency: " << endl;
		//for (auto a : adjacency_VV)
		//	log_list(loglevel, a, "", false);


		for (int j = 0; j < boundary.size(); j++)
		{
			int bi = boundary[j];
			indices.erase(indices.begin() + bi - j);
			adjacency_VV.erase(adjacency_VV.begin() + bi - j);

			write_log(loglevel) << j << ": delete bi = " << bi << ", which is at position " << bi - j << endl;
		}

		//write_log(loglevel) << endl << "AFTER erasing boundary: " << endl;
		//log_list(loglevel, boundary, "boundary: ", false);
		//log_list(loglevel, indices, "indices: ", false);
		//write_log(loglevel) << "adjacency: " << endl;
		//for (auto a : adjacency_VV)
		//	log_list(loglevel, a, "", false);



		int offset = cumulative_offset_vertices[i];
		//add_cwise(offset, boundary);
		add_cwise(offset, indices);
		add_cwise(offset, adjacency_VV);

		write_log(loglevel) << endl << "AFTER offset: " << endl;

		//log_list(loglevel, boundary, "boundary: ", false);
		log_list(loglevel, indices, "indices: ", false);
		write_log(loglevel) << "adjacency: " << endl;
		for (auto a : adjacency_VV)
			log_list(loglevel, a, "", false);


		inner_indices.insert(inner_indices.end(), indices.begin(), indices.end());
		inner_indices_adjacency.insert(inner_indices_adjacency.end(), adjacency_VV.begin(), adjacency_VV.end());



		//write_log(0) << endl << "CHECK: " << endl;
		//std::vector<int> all_indices(mesh->V.rows());
		//std::iota(all_indices.begin(), all_indices.end(), 0);
		//log_list(0, all_indices, "all_indices: ", false);

		//std::vector<int> local_inner_vertices(mesh->V.rows());
		//auto iterator = std::set_difference(all_indices.begin(), all_indices.end(), boundary.begin(), boundary.end(), local_inner_vertices.begin());
		//local_inner_vertices.resize(iterator - local_inner_vertices.begin());
		//log_list(0, local_inner_vertices, "local_inner_vertices: ", false);

	}

	flatten_indices(inner_indices, number_concatenated_vertices, data_model.target_developable.V.cols(), flattend_inner_indices);
	flatten_indices(inner_indices_adjacency, number_concatenated_vertices, data_model.target_developable.V.cols(), flattend_inner_indices_adjacency);
}
*/

/*
void DevelopableResultController::get_boundary_on_concatenated_developable(const std::vector<std::vector<int>>& global_patch_correspondance, Eigen::VectorXi& out_concatenated_boundary_indices, Eigen::MatrixXd& out_concatenated_boundary_vertices)
{
	const int loglevel = 5;

	vector<int> global_flat_boundary = meshhelper::boundary_loop_flat(data_model.target_developable.F);
	vector<int> concatenated_boundary_indices;
	out_concatenated_boundary_vertices.resize(global_flat_boundary.size() * 3, 3);
	//Eigen::MatrixXd concatenated_boundary_vertices(global_flat_boundary.size() * 3, 3);

	inner_vertices_mask.resize(number_concatenated_vertices);
	inner_vertices_mask.setOnes();

	for (int bi : global_flat_boundary)
	{
		vector<int> patch_indices = global_patch_correspondance[bi];
		write_log(loglevel) << bi << ": patch_indices: "; log_list(loglevel, patch_indices, "", false);

		Eigen::RowVectorXd vertex = data_model.target_developable.V.row(bi);

		for (int pi : patch_indices)
		{
			int mesh_index_b = index_of(pi, data_model.assigned_labels);
			int local_index_b = get_closest_vertex(data_model.developable_parts[mesh_index_b].V, vertex);
			int concatenated_index_b = local_index_b + cumulative_offset_vertices[pi];

			concatenated_boundary_indices.push_back(concatenated_index_b);
			//concatenated_boundary_vertices.row(concatenated_boundary_indices.size() - 1) = vertex;
			out_concatenated_boundary_vertices.row(concatenated_boundary_indices.size() - 1) = data_model.developable_parts[mesh_index_b].V.row(local_index_b);

			inner_vertices_mask(concatenated_index_b) = 0;
		}
	}

	out_concatenated_boundary_vertices.conservativeResize(concatenated_boundary_indices.size(), Eigen::NoChange);
	out_concatenated_boundary_indices = Eigen::Map<Eigen::VectorXi>(concatenated_boundary_indices.data(), concatenated_boundary_indices.size());

	//write_log(0) << "number_concatenated_V = " << number_concatenated_V << ", inner_vertices_mask: " << endl << inner_vertices_mask << endl;
}
*/

/*
void DevelopableResultController::split_vertices(std::vector<int>& boundary_indices, std::vector<int>& inner_indices, std::vector<std::vector<int>>& inner_indices_adjacency)
{
	boundary_indices.clear();
	inner_indices.clear();
	inner_indices_adjacency.clear();

	for (int i = 0; i < data_model.developable_parts.size(); i++)
	{
		Mesh* mesh = &data_model.developable_parts[i];

		std::vector<int> indices(mesh->V.rows());
		std::iota(indices.begin(), indices.end(), 0);

		vector<vector<int>> adjacency_VV;
		igl::adjacency_list(mesh->F, adjacency_VV, true);

		std::vector<int> boundary = meshhelper::boundary_loop_flat(mesh->F);
		std::sort(boundary.begin(), boundary.end());


		for (int j = 0; j < boundary.size(); j++)
		{
			int bi = boundary[j];
			indices.erase(indices.begin() + bi - j);
			adjacency_VV.erase(adjacency_VV.begin() + bi - j);
		}

		int offset = cumulative_offset_vertices[i];
		add_cwise(offset, boundary);
		add_cwise(offset, indices);
		add_cwise(offset, adjacency_VV);

		boundary_indices.insert(boundary_indices.end(), boundary.begin(), boundary.end());
		inner_indices.insert(inner_indices.end(), indices.begin(), indices.end());
		inner_indices_adjacency.insert(inner_indices_adjacency.end(), adjacency_VV.begin(), adjacency_VV.end());
	}
}
*/

void stitch_mesh(Developable& developable_model) 
{
	Eigen::MatrixXd V = developable_model.concatenated.V;
	Eigen::MatrixXi F = developable_model.concatenated.F;
	
	Eigen::MatrixXi C;
	igl::facet_components(F, C);

	// snap vertices to each other
	Eigen::MatrixXi seams = developable_model.seam_pairs;
	for (auto i = 0; i < seams.rows(); i++) 
		V.row(seams(i, 0)) = V.row(seams(i, 1));

	Eigen::MatrixXd newV;
	Eigen::MatrixXi newVI, newVJ, newF;
	igl::remove_duplicate_vertices(V, F, 1e-7, newV, newVI, newVJ, newF);

	// count connected components of the stitched mesh
	igl::facet_components(newF, C);
	write_log(4) << "stiched result mesh: " << newV.rows() << " vertices and " << newF.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << endl;

	developable_model.result.V = newV;
	developable_model.result.F = newF;
}

void flatten_single_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& uv)
{
	Eigen::VectorXi bnd;
	Eigen::MatrixXd bnd_uv;
	igl::boundary_loop(F, bnd);
	igl::map_vertices_to_circle(V, bnd, bnd_uv);

	igl::harmonic(V, F, bnd, bnd_uv, 1, uv);
	if (igl::flipped_triangles(uv, F).size() != 0) {
		igl::harmonic(F, bnd, bnd_uv, 1, uv); // use uniform laplacian
	}

	igl::SLIMData sData;
	sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
	
	Eigen::VectorXi b; 
	Eigen::MatrixXd bc;
	igl::slim_precompute(V, F, uv, sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);
	igl::slim_solve(sData, 20); // 20 iters
	uv = sData.V_o;
}

void flatten_patches(Developable& developable_model)
{
	Eigen::MatrixXd V = developable_model.concatenated.V;
	Eigen::MatrixXi F = developable_model.concatenated.F;

	developable_model.developable_parts_2D.clear();
	developable_model.developable_parts_2D.resize(developable_model.developable_parts.size());

	for(int i = 0; i < developable_model.developable_parts.size(); i++)
	{
		Eigen::MatrixXd subV = developable_model.developable_parts[i].V;
		Eigen::MatrixXi subF = developable_model.developable_parts[i].F;

		Eigen::MatrixXd subUV;
		flatten_single_mesh(subV, subF, subUV);

		// To save it as an OBJ we need to have 3 columns, i.e. x,y,z coordinates, so we will just save the 'z' as 0
		Eigen::MatrixXd subUV3d(subUV.rows(), 3);
		subUV3d.setZero(); 
		subUV3d.col(0) = subUV.col(0); 
		subUV3d.col(1) = subUV.col(1);

		Mesh mesh2d;
		mesh2d.V = subUV3d;
		mesh2d.F = subF;
		developable_model.developable_parts_2D[i] = mesh2d;
	}
}