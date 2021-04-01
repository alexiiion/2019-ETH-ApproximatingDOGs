#include "OptimizationController.h"

#include <igl/point_mesh_squared_distance.h>
//#include <igl/triangle_triangle_adjacency.h>
//#include <igl/adjacency_list.h>
#include <igl/Timer.h>

#include "Logger.h"
#include "Utils.h"
#include "MeshController.h"

#include "PatchModel.h"
#include "ConstraintsStrategy.h"
#include "Coverage.h"


//#include "Initialization/GeodesicSurfacesController.h"



using namespace std;


void OptimizationController::initialize()
{
	if (is_initialized)
		return;


	//TODO delete, only a test
	//GeodesicSurfacesController* surface_initialization = new GeodesicSurfacesController(data_model.target);
	//GeodesicSurfacesController* surface_initialization = new GeodesicSurfacesController(data_model.target, data_model.geodesics, data_model.ruled_surfaces);
	//int n = data_model.ruled_surfaces.settings.number_ruled_surfaces;
	//surface_initialization->build_surfaces(n);



	data_model.is_optimizing_DOG = true;
	data_model.is_optimizing_mesh = false;

	should_adapt_weights = false;
	do_update_objectives = true;

	valid_DOG_iterations = 0;
	max_valid_DOG_iterations = 10;
	max_iterations = 50;

	initial_settings = settings;

	bool success = geodesics_controller->initialize();
	if (!success)
		return;

	is_initialized = true;
	data_model.optimization_settings.should_converge = false;
	//add_patch();
}

void OptimizationController::run_optimiztaion()
{
	if (!data_model.is_optimizing)
		return;

	if (data_model.is_optimizing_DOG)
		run_DOG_optimiztaion();
	else if (data_model.is_optimizing_mesh)
		run_mesh_optimiztaion();
}

void OptimizationController::run_DOG_optimiztaion()
{
	if (!is_initialized)
		return;

	if (is_first)
	{
		init_time = timer.getElapsedTime(); 
		is_first = false;
	}

	DataModel::is_local_minimum = false;

	if(data_model.patches.size() < 1)
		add_patch();

	Patch* patch = data_model.patches[current_patch];
	patch->update();
	
	if (patch->constraints.paired_points_from.rows() < 1)
	{
		//finish last constraints (eg reset altered optimization weights)
		bool can_update = patch->get_current_strategy()->can_update();

		//add new patch
		add_patch();
		patch = data_model.patches[current_patch];
		patch->update();
	}

	//update_coverage();
	optimize(*patch);

	update_patch_rendering();
	//update_precise_patch_rendering();
	if (patch->has_strategy_changed && data_model.is_pausing_on_changing_constraint_strategy)
		data_model.is_optimizing = false;

	process_optimization_result(*patch);
}

void OptimizationController::run_mesh_optimiztaion()
{
	//if (result_controller)
	//	delete result_controller;

	//if (!result_controller)
	//	result_controller = new DevelopableResultController(data_model.target, data_model.patches, data_model.is_optimizing, data_model.viewer);
	
	result_controller->run_optimiztaion();
}

void OptimizationController::reset_mesh_optimiztaion()
{
	if(result_controller)
		result_controller->reset();
}

void OptimizationController::try_add_patch()
{
	if (!is_initialized)
		return;

	//if (data_model.target.is_covered())
	//if (data_model.patch_coverage->is_covered())
	//{
	//	//merge_patches();
	//	data_model.is_optimizing = false;

	//	double t = timer.getElapsedTime();
	//	write_log(0) << endl << endl << "Finished approximating DOGS -- ELAPSED TIME: " << t - init_time << endl << endl;

	//	return;
	//}

	add_patch();
}

void OptimizationController::add_patch()
{
	if (!is_initialized)
		return;

	//update_precise_patch_rendering();

	//meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

	//check_coverage();

	Eigen::MatrixXd target_curve;
	Eigen::MatrixXd wrapper_curve;
	int vertex;
	bool success = geodesics_controller->get_next_constraints(target_curve, wrapper_curve, vertex);


	//bool init_by_random_patches; if (init_by_random_patches && !success) do graph cut
	assert(success);
	//if (!success)
	//{
	//	check_coverage_holes();
	//	//TODO if found holes, add geodesics to data_model
	//	success = geodesics_controller->get_next_constraints(target_curve, wrapper_curve);
	//}

	//no more constraints, therefore we are done with the DOGs
	if (!success)
	{
		double t = timer.getElapsedTime();
		write_log(1) << endl << endl << "Finished approximating -- ELAPSED TIME: " << t - init_time << endl << endl;

		data_model.is_optimizing = false;
		data_model.is_optimizing_DOG = false;
		data_model.is_optimizing_mesh = true;
		is_first = true;


		//meshhelper::update_target_developable_optimized(data_model);

		return;
	}


	if(data_model.patches.size() > 0)
		write_log(1) << endl << endl << "Finished approximating PATCH after: " << timer.getElapsedTime() - patch_init_time << endl << endl;
	patch_init_time = timer.getElapsedTime();

	//Patch *patch = new Patch(data_model, target_curve, wrapper_curve);
	Patch *patch = new Patch(data_model, target_curve, vertex);
	data_model.patches.push_back(patch);
	current_patch = data_model.patches.size() - 1;

	do_update_objectives = true;

	//meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

	if (data_model.is_pausing_on_adding_patch)
		data_model.is_optimizing = false;

	settings.weight_fitting_energy = initial_settings.weight_fitting_energy;
}

void OptimizationController::update_position_constraints()
{
	//TODO do this in view, so I get which patch should be updated
	
	if (data_model.patches.size() < 1)
		return;

	data_model.patches[current_patch]->constraints.update_positions(); 
}

void OptimizationController::optimize(Patch& patch)
{
	if (!data_model.is_optimizing)
		return;

	write_log(3) << endl << "----- optimizing for DOG constraints -----" << endl << endl;
	write_log(3) << endl << "-- solving constraints" << endl << endl;


	Eigen::VectorXd x0;
	mat2_to_vec(patch.wrapper.V, x0);

	do_update_objectives = do_update_objectives || patch.has_strategy_changed || patch.current_strategy > 0;
	if(do_update_objectives)
		update_objectives(patch);
	else
		update_target_positions(patch);

	Eigen::VectorXd x;
	double obj = solver->solve_constrained(x0, *composite_objective, *dog_constraints, x);
	write_log(4) << "     obj = " << obj << endl;

	if (x.hasNaN())
	{
		write_log(1) << endl << "ERROR at optimization: result has NAN!" << endl << endl;
		data_model.is_optimizing = false;
		return;
	}

	vec_to_mat2(x, patch.wrapper.V);

	//double change_after_solving = (x - x0).squaredNorm();
	//write_log(4) << endl << "change_after_solving: " << change_after_solving << endl << endl;
	

  	patch.constraints.update_positions(); //update fitting energy to check for convergence

	data_model.output_fitting_energy = data_model.output_fitting_energy / patch.constraints.paired_points_from.rows(); //normalized fitting error
	//do_update_objectives = DataModel::output_DOG_objective < data_model.optimization_settings.convergence_deve lopablity_threshold;

	patch.DOG_error = data_model.output_DOG_objective;
	patch.fitting_error = data_model.output_fitting_energy;
	patch.is_local_minimum = data_model.is_local_minimum;
}

void OptimizationController::update_objectives(Patch& patch)
{
	write_log(4) << endl << "-------------" << endl << "UPDATE objectives" << endl << "-------------" << endl << endl;

	delete dog_constraints;

	delete composite_objective;
	delete fitting_constraints;
	delete fitting_objective;
	delete length_constraints;
	delete length_objective;
	delete bending_energy;
	delete iso_energy;

	composite_objective = new CompositeObjective();

	Eigen::MatrixXi wrapper_barycentric_indices_flat;
	Eigen::MatrixXd wrapper_barycentric_weights_flat;
	Eigen::VectorXd target_vertices_flat;
	patch.constraints.get_flattend_constraints(wrapper_barycentric_indices_flat, wrapper_barycentric_weights_flat, target_vertices_flat);

	fitting_constraints = new BarycentricPositionalConstraints(wrapper_barycentric_indices_flat, wrapper_barycentric_weights_flat, target_vertices_flat);
	write_log(0) << "wrapper_barycentric_indices_flat.maxCoeff() = " << wrapper_barycentric_indices_flat.maxCoeff() << endl;
	fitting_objective = new QuadraticConstraintsSumObjective(*fitting_constraints, patch.wrapper.initial_x0);
	composite_objective->add_objective(fitting_objective, settings.weight_fitting_energy);

	length_constraints = new LengthRegularizerConstraints(patch.wrapper.quad_topology);
	length_objective = new QuadraticConstraintsSumObjective(*length_constraints, patch.wrapper.initial_x0); 
	composite_objective->add_objective(length_objective, settings.weight_regularizing_energy);

	bending_energy = new SimplifiedBendingObjective(patch.wrapper.quad_topology, patch.wrapper.initial_x0);
	composite_objective->add_objective(bending_energy, settings.weight_bending_energy);

	iso_energy = new IsometryObjective(patch.wrapper.quad_topology, patch.wrapper.initial_x0);
	composite_objective->add_objective(iso_energy, settings.weight_isometry_energy);

	dog_constraints = new DogConstraints(patch.wrapper.quad_topology);

	delete solver;
	solver = new EqSQP(solver_target_developability_error, solver_max_developability_error, solver_max_iterations, solver_merit_p);

	do_update_objectives = false;
}

void OptimizationController::update_target_positions(Patch& patch)
{
	Eigen::MatrixXi wrapper_barycentric_indices_flat;
	Eigen::MatrixXd wrapper_barycentric_weights_flat;
	Eigen::VectorXd target_vertices_flat;
	patch.constraints.get_flattend_constraints(wrapper_barycentric_indices_flat, wrapper_barycentric_weights_flat, target_vertices_flat);

	fitting_constraints->update_coords(target_vertices_flat);
}

void OptimizationController::process_optimization_result(Patch& patch)
{
	if (should_adapt_weights)
		process_result_adaptive(patch);
	else
		process_result_static(patch);
}

void OptimizationController::process_result_static(Patch & patch)
{
	if (patch.is_done_approximating() || patch.constraints.constrained_barycentric_indices.rows() < 1)
	{
		do_update_objectives = false;
		try_add_patch();
	}

	//if (patch.has_converged())
	//{
	//	write_log(3) << endl << "----- has CONVERGED -----" << endl << endl;
	//	data_model.iterations_since_converge = 0;
	//
	//	if (patch.is_last_strategy())
	//		try_add_patch();
	//}
	//else if (DataModel::is_local_minimum)
	//{
	//	write_log(3) << endl << "----- is in local MINIMUM? moving on... -----" << endl << endl;
	//	data_model.iterations_since_converge = 0;
	//}
	else
	{
		data_model.iterations_since_converge++;
	}
}

void OptimizationController::process_result_adaptive(Patch & patch)
{
	//TODO further ideas for adapting weights: 
	//(1) use smoothness metric for a smoothing step, i.e. high bending energy when not smooth enough
	//(2) at last step of constraints strategy, relax isometry energy to allow for stretch as long as we keep good fit

	if (data_model.output_DOG_objective < settings.convergence_developablity_threshold)
		valid_DOG_iterations++;
	else
		valid_DOG_iterations = 0;

	if (patch.has_converged()) //is DOG && is fitted
	{
		write_log(3) << endl << "----- has CONVERGED -----" << endl << endl;
		data_model.iterations_since_converge = 0;

		if (settings.weight_fitting_energy < initial_settings.weight_fitting_energy)
			settings.weight_fitting_energy = initial_settings.weight_fitting_energy;

		if (patch.is_done_approximating())
		{
			write_log(3) << "----- try adding new PATCH -----" << endl << endl;
			try_add_patch();
		}
	}
	else if (patch.is_last_strategy())
	{
		write_log(3) << endl << "----- has NOT converged, is last strategy -----" << endl;
		if (DataModel::is_local_minimum)
		{
			write_log(3) << endl << "----- is local MINIMUM " << endl;

			if (data_model.output_DOG_objective < settings.convergence_developablity_threshold)  //is DOG
			{
				write_log(3) << "&& is DOG --> increasing fitting weight -----" << endl << endl;
				settings.weight_fitting_energy *= 1.1;
			}
			else
			{
				write_log(3) << "&& is NOT DOG --> decreasing fitting weight -----" << endl << endl;
				settings.weight_fitting_energy *= 0.85;
			}

			//data_model.iterations_since_converge = 0;
			data_model.iterations_since_converge++;
		}
		else if (data_model.iterations_since_converge > max_iterations) // over max iterations
		{
			write_log(3) << "----- is over max iterations " << endl;

			if (data_model.output_DOG_objective < settings.convergence_developablity_threshold) //is DOG
			{
				write_log(3) << "&& is DOG --> try adding new PATCH -----" << endl << endl;

				try_add_patch();
				data_model.iterations_since_converge = 0;
			}
			else // is not DOG
			{
				write_log(3) << "&& is no DOG --> decreasing fitting weight -----" << endl << endl;
				settings.weight_fitting_energy *= 0.7;
				data_model.iterations_since_converge = 0;
				//data_model.iterations_since_converge++;
			}
		}
		else //not over max iterations
		{
			if (data_model.output_DOG_objective < settings.convergence_developablity_threshold) //is DOG
			{
				if (valid_DOG_iterations >= max_valid_DOG_iterations) //is longtime DOG
				{
					write_log(3) << "&& is LONGTIME DOG --> increasing fitting weight -----" << endl << endl;
					settings.weight_fitting_energy *= 1.1;
				}
			}
			//else // is not DOG
			//{
			//	write_log(3) << "----- is NOT DOG  --> decreasing fitting weight -----" << endl;
			//	data_model.weight_fitting_energy *= 0.9;
			//}
			
			data_model.iterations_since_converge++;
		}
	}
	else if (data_model.output_DOG_objective < settings.convergence_developablity_threshold) //is DOG
	{
		write_log(3) << endl << "----- has NOT converged, is not last strategy, but IS DOG -----" << endl;

		if (valid_DOG_iterations >= max_valid_DOG_iterations) //is longtime DOG
		{
			write_log(3) << "&& is LONGTIME DOG --> increasing fitting weight -----" << endl << endl;
			settings.weight_fitting_energy *= 1.1;
		}

		data_model.iterations_since_converge++;
	}
	else if(DataModel::is_local_minimum)
	{
		write_log(3) << endl << "----- is local MINIMUM " << endl;

		if (data_model.output_DOG_objective < settings.convergence_developablity_threshold)  //is DOG
		{
			write_log(3) << "&& is DOG --> increasing fitting weight -----" << endl << endl;
			settings.weight_fitting_energy *= 1.1;
		}
		else
		{
			write_log(3) << "&& is NOT DOG --> decreasing fitting weight -----" << endl << endl;
			settings.weight_fitting_energy *= 0.75;
		}

		data_model.iterations_since_converge++;
	}
	else
	{
		data_model.iterations_since_converge++;
	}

	/*
	if (data_model.output_DOG_objective < data_model.convergence_developablity_threshold)
		valid_DOG_iterations++;
	else
		valid_DOG_iterations = 0;

	if (patch.has_converged())
	{
		write_log(3) << endl << "----- has CONVERGED -----" << endl << endl;
		data_model.iterations_since_converge = 0;

		//data_model.weight_fitting_energy = initial_fitting_weight;
		if (data_model.weight_fitting_energy < initial_fitting_weight)
			data_model.weight_fitting_energy = clip(data_model.weight_fitting_energy * 1.5, 0, initial_fitting_weight);

		//data_model.weight_bending_energy *= 1.05;
		//data_model.weight_isometry_energy *= 1.03;

		//if (patch.is_last_strategy())
		if (patch.is_done_approximating())
		{
			write_log(3) << "----- try adding new PATCH -----" << endl << endl;
			try_add_patch();
		}
	}
	else if (DataModel::is_local_minimum)
	{
		write_log(3) << endl << "----- local MINIMUM --> changing weights. -----" << endl << endl;

		if (data_model.output_DOG_objective > data_model.convergence_developablity_threshold)
			data_model.weight_fitting_energy *= 0.75;
		else 
			data_model.weight_fitting_energy *= 1.1;

		data_model.iterations_since_converge++;
	}
	else if (valid_DOG_iterations >= max_valid_DOG_iterations)
	{
		if (patch.is_last_strategy())
		{
			write_log(3) << endl << "----- not converged, but developable --> try adding new PATCH -----" << endl << endl;
			//data_model.weight_fitting_energy *= 1.25;
			try_add_patch();
			data_model.iterations_since_converge = 0;
		}
		else
		{
			data_model.weight_fitting_energy *= 1.1;
		}
	}
	else if (data_model.iterations_since_converge > max_iterations && patch.is_last_strategy())
	{
		if (data_model.output_DOG_objective > data_model.convergence_developablity_threshold)
		{
			data_model.weight_fitting_energy *= 0.95;
			data_model.iterations_since_converge = 0;
		}
		else
		{
			//add_patch();
			try_add_patch();
			data_model.iterations_since_converge = 0;
		}
	}
	else
	{
		data_model.iterations_since_converge++;
	}
	*/
}

void OptimizationController::update_coverage()
{
	data_model.patch_coverage->compute(meshhelper::get_wrapper_meshes(data_model.patches));
}

void OptimizationController::update_patch_rendering(bool do_update_all)
{
	if (do_update_all || data_model.is_rendering_patch_coverage)
	{
		meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

		update_coverage();
		//update_fast_patch_rendering();
	}
	//else if (data_model.patches[current_patch]->wrapper.render_F.rows() != data_model.patches[current_patch]->wrapper.F.rows())
	//{
	//	for (int i = 0; i < data_model.patches.size(); i++)
	//	{
	//		WrapperMesh& wrapper = data_model.patches[i]->wrapper;
	//		wrapper.render_F = wrapper.F;
	//		wrapper.render_F_quad = wrapper.F_quad;
	//	}
	//}
}

/*
void OptimizationController::update_fast_patch_rendering()
{
	vector<int> patch_F3_lengths;
	vector<vector<int>> visible_faces;

	for (int i = 0; i < data_model.patches.size(); i++)
	{
		WrapperMesh& wrapper = data_model.patches[i]->wrapper;

		int number_F3 = wrapper.F.rows();
		visible_faces.push_back(vector<int>{});

		if (i == 0)
			patch_F3_lengths.push_back(number_F3);
		else
			patch_F3_lengths.push_back(number_F3 + patch_F3_lengths[patch_F3_lengths.size() - 1]);
	}

	for (int i = 0; i < data_model.target.covering_wrapper_faces_distances.rows(); i++)
	{
		double distance = data_model.target.covering_wrapper_faces_distances(i);

		if (distance > settings.coverage_threshold)
			continue;
		  
		int face_concatenated = i;
		int patch_index = 0;
		int faces_start = 0;

		//coverage faces are stored for concatenated patches, need to find the face index per patch
		for (patch_index; patch_index < patch_F3_lengths.size(); patch_index++)
		{
			if (face_concatenated >= faces_start && face_concatenated < patch_F3_lengths[patch_index])
				break;

			faces_start = patch_F3_lengths[patch_index];
		}

		int face = face_concatenated - faces_start;
		visible_faces[patch_index].push_back(face);
	}

	//set render faces
	for (int i = 0; i < data_model.patches.size(); i++)
	{
		WrapperMesh& wrapper = data_model.patches[i]->wrapper;

		vector<int> visible_quads;
		triangle_to_quad_faces(visible_faces[i], visible_quads);
		index_to_value(visible_quads, wrapper.F_quad, wrapper.render_F_quad);

		//convert from quads back to triangles to avoid half quads from showing
		vector<int> visible_triangles;
		quad_to_triangle_faces(visible_quads, visible_triangles);
		index_to_value(visible_triangles, wrapper.F, wrapper.render_F);
		//index_to_value(visible_faces[i], wrapper.F, wrapper.render_F);
	}
}
*/

/*
void OptimizationController::update_precise_patch_rendering()
{
	for (int i = 0; i < data_model.patches.size(); i++)
	{
		WrapperMesh& wrapper = data_model.patches[i]->wrapper;
		std::vector<Eigen::MatrixXd> intersection_curves;

		for (int j = 0; j < data_model.patches.size(); j++)
		{
			if (i == j)
				continue;

			WrapperMesh& other_wrapper = data_model.patches[j]->wrapper;

			Eigen::MatrixXd intersection_curve;
			bool does_intersect = find_intersection_curve(wrapper, other_wrapper, intersection_curve);

			if (does_intersect)
				intersection_curves.push_back(intersection_curve);

			write_log(0) << "wrapper: " << i << ", other: " << j << "   does intersect: " << boolalpha << does_intersect << endl;
		}

		update_render_mesh(*data_model.patches[i], intersection_curves);
	}
}

void OptimizationController::update_render_mesh(Patch& patch, const std::vector<Eigen::MatrixXd>& intersection_curves)
{
	if (intersection_curves.size() < 1)
		return;

	Eigen::MatrixXd closest_points_unused;
	Eigen::VectorXd sqrD_unused;

	Eigen::VectorXi triangle_face_indices;
	vector<int> all_triangles;

	for (int i = 0; i < intersection_curves.size(); i++)
	{
		igl::point_mesh_squared_distance(intersection_curves[i], patch.wrapper.V, patch.wrapper.F, sqrD_unused, triangle_face_indices, closest_points_unused);
		
		//append_matrix(intersection_curves[i], data_model.debug_points);

		vector<int> triangles(triangle_face_indices.data(), triangle_face_indices.data() + triangle_face_indices.rows());
		all_triangles.insert(all_triangles.end(), triangles.begin(), triangles.end());
	}


	vector<int> intersection_quads;
	triangle_to_quad_faces(all_triangles, intersection_quads);

	if (LOG_LEVEL >= 5)
	{
		write_log(0) << endl << "intersection_quads:" << endl;
		for (int i = 0; i < intersection_quads.size(); i++)
			write_log(0) << "[" << intersection_quads[i] << "]: " << patch.wrapper.F_quad.row(intersection_quads[i]) << endl;
		write_log(0) << endl << endl;
	}

	//vector<vector<int>> QQA;
	//quad_quad_adjacency(F, quad_topology, QQA);

	//int min_error_vertex_index;
	//find_min_error_point(patch, min_error_vertex_index);
	//
	//int min_error_face = patch.wrapper.quad_topology.A[min_error_vertex_index][0];
	int min_error_face = (patch.wrapper.quad_height-2) / 2 * (patch.wrapper.quad_width-1);
	//write_log(0) << endl << "min_error_face = " << min_error_face << endl << endl;
	//debug_show_label(patch.wrapper.V.row(patch.wrapper.F_quad(min_error_face, 0)), to_string(min_error_face), data_model);

	vector<int> component;

	//std::vector<std::vector<int>> QQA;
	//quad_quad_adjacency(patch.wrapper.F, patch.wrapper.quad_topology, QQA);
	//dfs(QQA, min_error_face, intersection_quads, component);

	dfs(patch.wrapper.quad_topology.QQA, min_error_face, intersection_quads, component);
	log_list(6, component, "component quads:", false);

	vector<int> visible_quads = component;

	//vector<int> visible_quads;
	//visible_quads.insert(visible_quads.begin(), component.begin(), component.end());
	//visible_quads.insert(visible_quads.end(), intersection_quads.begin(), intersection_quads.end());
	
	//vector<int> visible_quads(component.size() + intersection_quads.size());
	//std::sort(component.begin(), component.end());
	//std::sort(intersection_quads.begin(), intersection_quads.end());
	//auto iterator = std::set_union(component.begin(), component.end(), intersection_quads.begin(), intersection_quads.end(), visible_quads.begin());
	log_list(5, visible_quads, "visible quads:", false);


	vector<int> render_triangles;
	quad_to_triangle_faces(visible_quads, render_triangles);

	index_to_value(visible_quads, patch.wrapper.F_quad, patch.wrapper.render_F_quad);
	index_to_value(render_triangles, patch.wrapper.F, patch.wrapper.render_F);
}

bool OptimizationController::find_intersection_curve(const WrapperMesh& wrapper, const WrapperMesh& other_wrapper, Eigen::MatrixXd& out_intersection_curve)
{
	std::vector<WrapperMesh> wrappers = { wrapper, other_wrapper };

	int intersect_direction = 0;
	double max_distance = 1e3;
	double max_distance_change = 1e3;
	double threshold = 1e-3;

	int number_iterations = 0;
	int number_iterations_no_change = 0;
	const int max_number_iterations = 50;
	double intersection_threshold = 1e-3;
	double mean_distance = 1e3;

	Eigen::MatrixXd C;
	Eigen::MatrixXd P;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::VectorXi I;

	const int loglevel = 5;
	write_log(loglevel) << endl << "merge_patches: " << endl;

	while (max_distance > threshold)
	{
		if (number_iterations > max_number_iterations && mean_distance > intersection_threshold)
			return false;

		const int i_1 = intersect_direction;
		const int i_2 = intersect_direction == 0 ? 1 : 0;

		P = C.rows() == 0 ? wrappers[i_1].V : C;
		V = wrappers[i_2].V;
		F = wrappers[i_2].F;

		Eigen::VectorXd sqrD;
		igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

		double current_max_distance = sqrt(sqrD.maxCoeff());
		max_distance_change = abs(max_distance - current_max_distance);
		max_distance = current_max_distance;
		intersect_direction = i_2;

		if (max_distance_change > 1e-3)
			number_iterations_no_change = 0;
		else
			number_iterations_no_change++;

		if (number_iterations_no_change > 5)
			break; //meaning that there is an intersection??
			//return false;

		mean_distance = sqrD.mean();
		number_iterations++;

		write_log(loglevel) << number_iterations << ": max_distance = " << max_distance << " sqrD.mean = " << mean_distance << " (change=" << max_distance_change << " iterations=" << number_iterations_no_change << ")" << endl;
	}

	//out_intersection_curve.resize(P.rows() + C.rows(), 3);
	//out_intersection_curve << P, C;

	out_intersection_curve = C;
	//remove_straight_segments(out_intersection_curve);

	//data_model.debug_points = out_intersection_curve;
	//debug_show_curve(out_intersection_curve, data_model);
	//data_model.debug_faces = I;

	return true;
}

void OptimizationController::find_min_error_point(const Patch& patch, int& out_vertex_index)
{
	//find face on wrapper with smallest position error
	double min_distance = 1e3;
	int min_vertex_index;

	for (int i = 0; i < patch.constraints.paired_points_from.rows(); i++)
	{
		auto point_to = patch.constraints.paired_points_to.row(i);
		auto distance = (patch.constraints.paired_points_from.row(i) - point_to).squaredNorm();
		write_log(0) << i << ": distance = " << distance << endl;

		if (distance > min_distance)
			continue;

		min_distance = distance;
		(patch.constraints.pair_to_V->rowwise() - point_to).rowwise().squaredNorm().minCoeff(&min_vertex_index);

		write_log(0) << "---> is minimum! vertex index = " << min_vertex_index << endl;
	}

	out_vertex_index = min_vertex_index;
	write_log(0) << endl << "min_distance = " << min_distance << ", vertex index = " << out_vertex_index << endl;
}
*/