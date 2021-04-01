#include "PositionConstraintsModel.h"

#include <igl/point_mesh_squared_distance.h>

#include "PairingModes.h"
#include "CoordinateConverter.h"
#include "Utils.h"
#include "Logger.h"


using namespace std;

void flatten_barycentrix_indices(const Eigen::MatrixXi& matrix, const int size, const int element_offset, Eigen::MatrixXi& out_flat);
void flatten_barycentrix_matrix(const Eigen::MatrixXd& matrix, const int size, const int element_offset, Eigen::MatrixXd& out_flat);
void flatten_matrix(const Eigen::MatrixXd& matrix, const int size, Eigen::VectorXd& out_flat);

PositionConstraints::PositionConstraints(Mesh& target, WrapperMesh& wrapper, OptimizationSettings& optimization_settings)
	: target(target), wrapper(wrapper), optimization_settings(optimization_settings)
{
	pairing_target = PairingTarget::to_SURFACE;
	set_pairing_direction(PairingDirection::TARGET_to_DOG);
}

//void PositionConstraints::initialize(Mesh* target, WrapperMesh* wrapper, Eigen::MatrixXd& from, Eigen::MatrixXd& to)
//{
//	this->target = target;
//	this->wrapper = wrapper;
//
//	pairing_target = PairingTarget::to_SURFACE;
//	set_pairing_direction(PairingDirection::TARGET_to_DOG);
//}

//void PositionConstraints::initialize(Mesh* target, WrapperMesh* wrapper, OptimizationSettings* optimization_settings)
//{
//	this->target = target;
//	this->wrapper = wrapper;
//	this->optimization_settings = optimization_settings;
//
//	pairing_target = PairingTarget::to_SURFACE;
//	set_pairing_direction(PairingDirection::TARGET_to_DOG);
//}

void PositionConstraints::add(Eigen::MatrixXd& from)
{
	Eigen::MatrixXd to;
	find_closest_points(from, to);

	add(from, to);
}

void PositionConstraints::add(Eigen::MatrixXd& from, Eigen::MatrixXd& to)
{
	add_block_settings(from.rows());

	if (!use_previous_constraints)
		for(int i = 0; i < blocks.size() - 1; i++)
			blocks[i].is_in_use = false;
	
	current_block++;
	number_all_constraints += from.rows();

	Eigen::MatrixXi bary_indices;
	Eigen::MatrixXd bary_weigths;
	Eigen::MatrixXd& wrapper_points = PairingDirection::TARGET_to_DOG ? to : from;
	//CoordinateConverter::barycentric_coords_from_points(wrapper_points, wrapper->V, wrapper->F, bary_indices, bary_weigths);
	CoordinateConverter::barycentric_coords_from_points(wrapper_points, wrapper.V, wrapper.F, bary_indices, bary_weigths);

	blocks[current_block].points_from = from;
	blocks[current_block].points_to = to;
	blocks[current_block].barycentric_indices = bary_indices;
	blocks[current_block].barycentric_weights = bary_weigths;

	update_positions();
}

void PositionConstraints::update_positions()
{
	for (BlockConstraintsConfig& block : blocks)
		update_block(block);

	construct_optimized_constraints();
}

void PositionConstraints::update_block(BlockConstraintsConfig& block)
{
	if (!block.is_in_use && !use_all_constraints)
		return;


	set_pairing_direction(block.pairing_direction);

	if (block.has_constant_pairing) //refresh to-points from barycentric coords (all in block!)
	{
		Eigen::MatrixXd &points_on_wrapper = block.pairing_direction == PairingDirection::TARGET_to_DOG ? block.points_to : block.points_from;
		//CoordinateConverter::points_from_barycentric_coords(block.barycentric_indices, block.barycentric_weights, wrapper->V, points_on_wrapper); //TODO check referencing!
		CoordinateConverter::points_from_barycentric_coords(block.barycentric_indices, block.barycentric_weights, wrapper.V, points_on_wrapper); //TODO check referencing!
	}
	else //store points from barycentric in all_points_to
	{
		Eigen::MatrixXd &surface = pairing_target == PairingTarget::to_VERTEX ? block.points_to : block.points_from;
		Eigen::VectorXi closest_faces;
		find_closest_points(surface, block.points_to, closest_faces);

		Eigen::MatrixXd &points_on_wrapper = block.pairing_direction == PairingDirection::TARGET_to_DOG ? block.points_to : block.points_from;
		//CoordinateConverter::barycentric_coords_from_points(points_on_wrapper, closest_faces, wrapper->V, wrapper->F, block.barycentric_indices, block.barycentric_weights);
		CoordinateConverter::barycentric_coords_from_points(points_on_wrapper, closest_faces, wrapper.V, wrapper.F, block.barycentric_indices, block.barycentric_weights);
	}

	block.filtered_constraints.clear();
	block.filtered_constraints.resize(block.points_from.rows(), true);

	if (block.is_removing_outliers)
		label_outliers(block);

	if (block.is_limiting_face_constraints)
		label_exceeding_face_constraints(block);

	if (block.is_removing_normal_outliers)
		label_normal_outliers(block);
}

void PositionConstraints::construct_optimized_constraints()
{
	paired_points_from.resize(number_all_constraints, 3);
	paired_points_to.resize(number_all_constraints, 3);

	constrained_barycentric_indices.resize(number_all_constraints, 3);
	constrained_barycentric_weights.resize(number_all_constraints, 3);

	int added_constraints = 0;

	for (BlockConstraintsConfig& block : blocks)
	{
		if (!block.is_in_use && !use_all_constraints)
			continue;

		for (int j = 0; j < block.filtered_constraints.size(); j++)
		{
			if (!block.filtered_constraints[j])
				continue;

			paired_points_from.row(added_constraints) = block.points_from.row(j);
			paired_points_to.row(added_constraints) = block.points_to.row(j);

			constrained_barycentric_indices.row(added_constraints) = block.barycentric_indices.row(j);
			constrained_barycentric_weights.row(added_constraints) = block.barycentric_weights.row(j);

			added_constraints++;
		}
	}

	if (added_constraints == 0)
	{
		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> number of constraints are 0!" << endl << endl;
		//abort();
	}

	paired_points_from.conservativeResize(added_constraints, Eigen::NoChange);
	paired_points_to.conservativeResize(added_constraints, Eigen::NoChange);

	constrained_barycentric_indices.conservativeResize(added_constraints, Eigen::NoChange);
	constrained_barycentric_weights.conservativeResize(added_constraints, Eigen::NoChange);
}

void PositionConstraints::set_current_constraints(Eigen::MatrixXd& from)
{
	if (current_block < 0)
		return;

	number_all_constraints += from.rows() - blocks[current_block].points_from.rows();
	blocks[current_block].points_from.resize(number_all_constraints, Eigen::NoChange);
	blocks[current_block].points_from = from;
	update_positions();
}

void PositionConstraints::get_flattend_constraints(Eigen::MatrixXi& out_wrapper_barycentric_indices_flat, Eigen::MatrixXd& out_wrapper_barycentric_weights_flat, Eigen::VectorXd& out_target_points_flat)
{
	//const int number_constrained_points = paired_points_from.rows();
	const int number_constrained_points = paired_points_to.rows();
	//flatten_barycentrix_indices(constrained_barycentric_indices, number_constrained_points, wrapper->V.rows(), out_wrapper_barycentric_indices_flat);
	flatten_barycentrix_indices(constrained_barycentric_indices, number_constrained_points, wrapper.V.rows(), out_wrapper_barycentric_indices_flat);
	flatten_barycentrix_matrix(constrained_barycentric_weights, number_constrained_points, 0, out_wrapper_barycentric_weights_flat);

	Eigen::MatrixXd &points_on_target = blocks[current_block].pairing_direction == PairingDirection::TARGET_to_DOG ? paired_points_from : paired_points_to;
	flatten_matrix(points_on_target, number_constrained_points, out_target_points_flat);

	write_log(5) << "barycentric indides flat: " << endl << out_wrapper_barycentric_indices_flat << endl;
	write_log(5) << "barycentric weights flat: " << endl << out_wrapper_barycentric_weights_flat << endl;
	write_log(5) << "points on target flat: " << endl << out_target_points_flat << endl << endl;
}

void PositionConstraints::remove_all_constraints()
{
	blocks.clear();
	current_block = -1;
	number_all_constraints = 0;

	paired_points_from.resize(number_all_constraints, Eigen::NoChange);
	paired_points_to.resize(number_all_constraints, Eigen::NoChange);

	constrained_barycentric_indices.resize(number_all_constraints, Eigen::NoChange);
	constrained_barycentric_weights.resize(number_all_constraints, Eigen::NoChange);
}

void PositionConstraints::label_outliers(BlockConstraintsConfig& block)
{
	const int number_constraints = block.points_from.rows();
	const double max_distance_sq = (optimization_settings.outlier_threshold) * (optimization_settings.outlier_threshold);

	for (int i = 0; i < number_constraints; i++)
	{
		const auto from = block.points_from.row(i);
		const auto to = block.points_to.row(i);

		if ((from - to).squaredNorm() > max_distance_sq)
			block.filtered_constraints[i] = false;
	}

	int active_constraints = std::count(block.filtered_constraints.begin(), block.filtered_constraints.end(), true);
	write_log(4) << "<" << __FUNCTION__ << "> active constraints: " << active_constraints << endl;
}

void PositionConstraints::label_normal_outliers(BlockConstraintsConfig& block)
{
	const int number_constraints = block.points_from.rows();
	//const double max_deviation = 0.7; //TODO add GlobalSettings::crease_threshold
	const double max_deviation = optimization_settings.outlier_normal_threshold;
	/*
	Mesh from = block.pairing_direction == PairingDirection::TARGET_to_DOG ? target  : wrapper;
	Mesh to   = block.pairing_direction == PairingDirection::TARGET_to_DOG ? wrapper : target;
	*/

	Eigen::MatrixXd from_V	 = block.pairing_direction == PairingDirection::TARGET_to_DOG ? target.V : wrapper.V;
	Eigen::MatrixXi from_F	 = block.pairing_direction == PairingDirection::TARGET_to_DOG ? target.F : wrapper.F;
	Eigen::MatrixXd from_NV  = block.pairing_direction == PairingDirection::TARGET_to_DOG ? target.normals_vertices : wrapper.normals_vertices;
	std::vector<std::vector<int>> from_adjacency_VF = block.pairing_direction == PairingDirection::TARGET_to_DOG ? target.adjacency_VF : wrapper.adjacency_VF;

	Eigen::MatrixXd to_V	 = block.pairing_direction == PairingDirection::TARGET_to_DOG ? wrapper.V : target.V ;
	Eigen::MatrixXi to_F	 = block.pairing_direction == PairingDirection::TARGET_to_DOG ? wrapper.F : target.F;
	Eigen::MatrixXd to_NV    = block.pairing_direction == PairingDirection::TARGET_to_DOG ? wrapper.normals_vertices : target.normals_vertices;
	std::vector<std::vector<int>> to_adjacency_VF = block.pairing_direction == PairingDirection::TARGET_to_DOG ? wrapper.adjacency_VF : target.adjacency_VF;

	for (int i = 0; i < number_constraints; i++)
	{
		//if the constraint was already filtered, then skip normal computation
		if (!block.filtered_constraints[i]) 
			continue; 

		//otherwise compute normal deviation
		const auto from = block.points_from.row(i);
		const auto to = block.points_to.row(i);

		const auto from_normal = get_surface_normal(from, from_V, from_F, from_NV, from_adjacency_VF);
		const auto to_normal   = get_surface_normal(to,   to_V,   to_F,   to_NV,   to_adjacency_VF);

		double normal_angle = angle(from_normal, to_normal);
		if (normal_angle > max_deviation)
			block.filtered_constraints[i] = false;

		//write_log(0) << i << ": filtered? " << !(normal_angle > max_deviation) << " (normal_angle = " << normal_angle << ")" << endl;
	}

	int active_constraints = std::count(block.filtered_constraints.begin(), block.filtered_constraints.end(), true);
	write_log(4) << "<" << __FUNCTION__ << "> active constraints: " << active_constraints << endl;
}

void PositionConstraints::label_exceeding_face_constraints(BlockConstraintsConfig& block)
{
	/*
	const int number_constraints = block.points_from.rows();
	vector<int> faces_with_constraints;

	for (int i = 0; i < number_constraints; i++)
	{
		if (block.filtered_constraints[i] == false)
			continue;

		const auto face = block.barycentric_indices.row(i);

		Eigen::MatrixXd::Index face_index;
		double d = (wrapper.F.rowwise() - face).rowwise().squaredNorm().minCoeff(&face_index);
		//double d = (wrapper->F.rowwise() - face).rowwise().squaredNorm().minCoeff(&face_index);

		auto iterator = std::find(faces_with_constraints.begin(), faces_with_constraints.end(), face_index);
		bool is_found = iterator != faces_with_constraints.end();

		if (is_found)
			block.filtered_constraints[i] = false; //face already has constraint
		else 
			faces_with_constraints.push_back(face_index);


	//	write_log(0) << "[" << i << "] has constraint? " << boolalpha << is_found << endl;
	//	write_log(0) << "   found face index: " << face_index << " (distance = " << d << ")" << endl;
	}

	//log_list(0, faces_with_constraints, "faces_with_constraints:", true, true);

	//write_log(0) << "block.filtered_constraints: " << endl;
	//for(int i = 0; i < block.filtered_constraints.size(); i++)
	//	write_log(0) << "[" << i << "]: " << boolalpha << block.filtered_constraints[i] << endl;


	block.filtered_constraints[number_constraints] = true; //add last constraint to pull boundary
	*/

	const int number_constraints = block.points_from.rows();
	vector<vector<int>> constraint_faces(wrapper.F.rows());
	vector<vector<double>> constraint_distances(wrapper.F.rows());

	//assign all constraints to face
	for (int i = 0; i < number_constraints; i++)
	{
		if (block.filtered_constraints[i] == false)
			continue;

		const auto face = block.barycentric_indices.row(i);

		Eigen::MatrixXd::Index face_index;
		double d = (wrapper.F.rowwise() - face).rowwise().squaredNorm().minCoeff(&face_index);

		const auto from = block.points_from.row(i);
		const auto to = block.points_to.row(i);

		constraint_faces[face_index].push_back(i);
		constraint_distances[face_index].push_back((from - to).squaredNorm());
	}

	for (int i = 0; i < constraint_faces.size(); i++)
	{
		//write_log(0) << "face [" << i << "] constraints: " << list_to_string(constraint_faces[i]) << endl;
		//write_log(0) << "face [" << i << "] distances  : " << list_to_string(constraint_distances[i]) << endl;

		if (constraint_faces[i].size() <= 1)
			continue;

		//on this face, find constraint with largest distance
		double distance = *std::min_element(constraint_distances[i].begin(), constraint_distances[i].end());
		if(should_foster_stretch)
			distance = *std::max_element(constraint_distances[i].begin(), constraint_distances[i].end());

		int index = index_of(distance, constraint_distances[i]);

		//remove all other constraints that are on this face from consideration
		for(int j = 0; j < constraint_faces[i].size(); j++)
			if(j != index)
				block.filtered_constraints[constraint_faces[i][j]] = false;
	}

	int active_constraints = std::count(block.filtered_constraints.begin(), block.filtered_constraints.end(), true);
	write_log(4) << "<" << __FUNCTION__ << "> active constraints: " << active_constraints << endl;
}

void PositionConstraints::add_block_settings(int block_size)
{
	BlockConstraintsConfig config;
	config.index = blocks.size();
	//config.start_index = all_points_from.rows();
	//config.block_size = block_size;

	config.is_in_use = true;
	config.has_constant_pairing = is_keeping_pairs_constant;
	config.is_removing_outliers = should_remove_outliers;
	config.is_removing_normal_outliers = should_remove_normal_outliers;
	config.is_limiting_face_constraints = should_limit_face_constraints;

	config.pairing_direction = pairing_direction;
	config.pairing_target = pairing_target;

	blocks.push_back(config);
}

void PositionConstraints::find_closest_points(const Eigen::MatrixXd& from, Eigen::MatrixXd& out_to)
{
	Eigen::VectorXi closest_faces;
	Eigen::VectorXd sq_d;
	Eigen::MatrixXd to;
	igl::point_mesh_squared_distance(from, *pair_to_V, *pair_to_F, sq_d, closest_faces, out_to);
}

void PositionConstraints::find_closest_points(const Eigen::MatrixXd& from, Eigen::MatrixXd& out_to, Eigen::VectorXi& out_closest_faces)
{
	Eigen::VectorXd sq_d;
	Eigen::MatrixXd to;
	igl::point_mesh_squared_distance(from, *pair_to_V, *pair_to_F, sq_d, out_closest_faces, out_to);
}

void PositionConstraints::set_pairing_direction(int pairing_direction)
{
	if (pairing_direction == PairingDirection::DOG_to_TARGET)
	{
		pair_from_V = &wrapper.V;
		pair_from_F = &wrapper.F;
		pair_to_V = &target.V;
		pair_to_F = &target.F;
		//pair_from_V = &wrapper->V;
		//pair_from_F = &wrapper->F;
		//pair_to_V = &target->V;
		//pair_to_F = &target->F;
	}
	else if (pairing_direction == PairingDirection::TARGET_to_DOG)
	{
		pair_from_V = &target.V;
		pair_from_F = &target.F;
		pair_to_V = &wrapper.V;
		pair_to_F = &wrapper.F;
		//pair_from_V = &target->V;
		//pair_from_F = &target->F;
		//pair_to_V = &wrapper->V;
		//pair_to_F = &wrapper->F;
	}

	this->pairing_direction = pairing_direction;
}

void PositionConstraints::invert_pairing_direction()
{
	write_log(3) << endl << "- invert pairing direction" << endl;

	int new_direction = pairing_direction == PairingDirection::TARGET_to_DOG ? PairingDirection::DOG_to_TARGET : PairingDirection::TARGET_to_DOG;

	string log = pairing_direction == PairingDirection::TARGET_to_DOG ? "TARGET_to_DOG --> DOG_to_TARGET" : "DOG_to_TARGET --> TARGET_to_DOG";
	write_log(4) << endl << log << endl << endl;


	set_pairing_direction(new_direction);
	blocks[current_block].pairing_direction = new_direction;

	auto temp_to = blocks[current_block].points_to;
	blocks[current_block].points_to = blocks[current_block].points_from;
	blocks[current_block].points_from = temp_to;
}

void flatten_barycentrix_indices(const Eigen::MatrixXi& matrix, const int size, const int element_offset, Eigen::MatrixXi& out_flat)
{
	const int num_cols = matrix.cols();
	out_flat.resize(size * num_cols, num_cols);

	Eigen::RowVectorXi ones(num_cols);
	ones.setOnes();

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < num_cols; j++) {
			Eigen::RowVectorXi vertex = matrix.row(i);
			out_flat.row(i + size * j) = vertex + j * element_offset*ones;

			write_log(6) << "flatten_barycentrix_indices["<< i << "][" << j << "](" << i + size * j << ") = " << out_flat.row(i + size * j) << endl;

			if (i + size * j >= size * num_cols)
				write_log(1) << "ERROR at flatten_barycentrix_indices" << endl;
		}
	}
}

void flatten_barycentrix_matrix(const Eigen::MatrixXd& matrix, const int size, const int element_offset, Eigen::MatrixXd& out_flat)
{
	const int num_cols = matrix.cols();
	out_flat.resize(size * num_cols, num_cols);

	Eigen::RowVectorXd ones(num_cols);
	ones.setOnes();

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < num_cols; j++) {
			Eigen::RowVectorXd vertex = matrix.row(i);
			out_flat.row(i + size * j) = vertex + j * element_offset*ones;

			write_log(6) << "flatten_barycentrix_matrix(" << i + size * j << ") = " << out_flat.row(i + size * j) << endl;
		}
	}
}

void flatten_matrix(const Eigen::MatrixXd& matrix, const int size, Eigen::VectorXd& out_flat)
{
	const int num_cols = matrix.cols();
	out_flat.resize(size * num_cols);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < num_cols; j++) {
			auto vertex = matrix(i, j);
			out_flat(i + j * size) = vertex;

			write_log(6) << "flatten_matrix(" << i + j * size << ") = " << out_flat(i + j * size) << endl;
		}
	}
}
