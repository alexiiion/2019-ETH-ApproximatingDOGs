#include "GrowTargetConstraints.h"

#include "PairingModes.h"
#include "Logger.h"
#include "Utils.h"

using namespace std;

void GrowTargetConstraints::configure()
{
	write_log(3) << endl << "    -- configuring grow target constraints " << endl;

	do_update = true;
	is_first_update = true;

	constraints.should_limit_face_constraints = true;
	constraints.should_remove_outliers = true;
	constraints.should_remove_normal_outliers = true;
	constraints.is_keeping_pairs_constant = false; //false
	constraints.use_previous_constraints = true; //true
	constraints.should_foster_stretch = false;

	previous_weight_fitting_energy = constraints.optimization_settings.weight_fitting_energy;
	previous_should_converge = constraints.optimization_settings.should_converge;

	//constraints.optimization_settings.weight_fitting_energy *= 2.0; //this was not used
	//constraints.optimization_settings.should_converge = true;


	//check since we are using normal outlier removal here, check if normals are oriented correctly
	int mid_patch = patch->wrapper.quad_width * patch->wrapper.quad_height / 2; //this works b/c the dog is constructed in an ordered way
	int mid_target = patch->target_curve.rows() % 2 == 0 ? patch->target_curve.rows() / 2 - 1 : patch->target_curve.rows() / 2;
	Eigen::RowVectorXd nt = constraints.target.normals_vertices.row(mid_target);
	Eigen::RowVectorXd np = patch->wrapper.normals_vertices.row(mid_patch);
	
	write_log(0) << "        mid patch index : " << mid_patch << ", mid_target: " << mid_target << ", normal target: " << nt << ", normal wrapper: " << np << endl;

	if (nt.dot(np) < 0)
		patch->wrapper.normals_vertices *= -1;


	//construct constraints
	constraints.pairing_target = PairingTarget::to_SURFACE;
	constraints.pairing_direction = PairingDirection::TARGET_to_DOG; //original setting was this one
	//constraints.pairing_direction = PairingDirection::DOG_to_TARGET;

	value_to_index(patch->target_curve, constraints.target.V, last_added_vertices);
	constrainted_vertices.insert(constrainted_vertices.begin(), last_added_vertices.begin(), last_added_vertices.end());
	constraints.add(patch->target_curve); //TODO is added twice at this point, could refactor if becomes a problem, but this curve is rather short

	//HACK!
	current_grow_iteration = 0;
	max_number_grow_iterations = patch->wrapper.quad_height * 1; //5

	DataModel::current_constraints_step = 0.0;
}

void GrowTargetConstraints::update()
{
	write_log(3) << endl << "    -- update grow target constraints " << endl;
	write_log(3)         << "       update:: current_grow_iteration = " << current_grow_iteration << endl << endl;

	//TODO check change in coverage stop when it is not changing enough anymore 
	if (current_grow_iteration >= max_number_grow_iterations-1)
		do_update = false;

 	//constraints.optimization_settings.weight_bending_energy *= 1.15; 

	//get neighbors of last added vertices
	vector<int> neighbors;
	for (int i = 0; i < last_added_vertices.size(); i++)
	{
		int vertex_index = last_added_vertices[i];
		vector<int> vertex_neighbors = constraints.target.adjacency_VV[vertex_index];

		neighbors.insert(neighbors.begin(), vertex_neighbors.begin(), vertex_neighbors.end());
	}

	//get the difference set between newly found neighbors and already constrained vertices
	std::vector<int> constrainted_vertices_copy(constrainted_vertices.size());
	std::copy(constrainted_vertices.begin(), constrainted_vertices.end(), constrainted_vertices_copy.begin());

	sort(constrainted_vertices_copy.begin(), constrainted_vertices_copy.end());
	sort(neighbors.begin(), neighbors.end());
	
	neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

	vector<int> difference(neighbors.size());
	auto iterator = std::set_difference(neighbors.begin(), neighbors.end(), constrainted_vertices_copy.begin(), constrainted_vertices_copy.end(), difference.begin());
	difference.resize(iterator - difference.begin());

	//if could not find any more unique vertices to constrain, then can't update anymore
	if (difference.size() <= 0)
	{
		do_update = false;
		return;
	}

	//set the unique, non overlapping vertices as constraints
	last_added_vertices = difference;
	constrainted_vertices.insert(constrainted_vertices.begin(), last_added_vertices.begin(), last_added_vertices.end());

	Eigen::MatrixXd from_points;
	index_to_value(constrainted_vertices, constraints.target.V, from_points);
	//index_to_value(last_added_vertices, constraints.target->V, from_points);
	constraints.set_current_constraints(from_points);
	if (constraints.paired_points_from.rows() < 1)
		do_update = false;

	current_grow_iteration++;
	DataModel::current_constraints_step++;
}

bool GrowTargetConstraints::can_update()
{
	write_log(3) << endl << "    -- can update grow target constraints? " << boolalpha << do_update << endl;

	if (!do_update)
	{
		constraints.optimization_settings.weight_fitting_energy = previous_weight_fitting_energy;
		constraints.optimization_settings.should_converge = previous_should_converge;
	}

	return do_update;
}

