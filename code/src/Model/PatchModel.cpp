#include "PatchModel.h"

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/procrustes.h>
#include <igl/point_mesh_squared_distance.h>

#include "CurveInterpolationConstraints.h"
#include "GrowTargetConstraints.h"
#include "CoverTargetConstraints.h"

#include "CoordinateConverter.h"
#include "MeshController.h"
#include "Logger.h"
#include "Utils.h"
#include "CurveHelper.h"

using namespace std;

const Eigen::RowVector3d sample_direction = Eigen::RowVector3d::UnitX();


Patch::Patch(DataModel& data_model, Eigen::MatrixXd& target_curve, Eigen::MatrixXd& wrapper_curve, int& vertex)
	: data_model(data_model), target_curve(target_curve), wrapper_curve(wrapper_curve), geodesic_vertex(vertex), constraints(data_model.target, wrapper, data_model.optimization_settings)
{
	add_wrapper();

	constraints_strategies.push_back(new CurveInterpolationConstraints(this, constraints));
	constraints_strategies.push_back(new GrowTargetConstraints(this, constraints));
	constraints_strategies.push_back(new CoverTargetConstraints(this, constraints));

	current_strategy = 0;
	constraints_strategies[current_strategy]->configure();
	renew_position_constraints();
	has_strategy_changed = true;

	DOG_error = 1e3;
	fitting_error = 1e3;
	is_local_minimum = true;
}
Patch::Patch(DataModel& data_model, Eigen::MatrixXd& target_curve, int& vertex)
	: data_model(data_model), target_curve(target_curve), geodesic_vertex(vertex), constraints(data_model.target, wrapper, data_model.optimization_settings)
{
	add_wrapper();


	////for overconstraining DOG scenario
	//constraints.should_limit_face_constraints = true;
	//constraints.should_remove_outliers = true;
	//constraints.should_remove_normal_outliers = false;
	//constraints.is_keeping_pairs_constant = false;
	//constraints.use_previous_constraints = false;
	//constraints_strategies.push_back(new CoverTargetConstraints(this, constraints));


	constraints_strategies.push_back(new CurveInterpolationConstraints(this, constraints));
	constraints_strategies.push_back(new GrowTargetConstraints(this, constraints));
	constraints_strategies.push_back(new CoverTargetConstraints(this, constraints));

	current_strategy = 0;
	constraints_strategies[current_strategy]->configure();
	renew_position_constraints();
	has_strategy_changed = true;

	DOG_error = 1e3;
	fitting_error = 1e3;
	is_local_minimum = true;
}

Patch::~Patch()
{
	for (auto strategy : constraints_strategies)
		delete strategy;

	constraints_strategies.clear();
}

bool Patch::is_done_approximating()
{
	if (constraints.paired_points_from.rows() < 1 || constraints.paired_points_from.rows() < 1)
		return true;

	bool can_update = constraints_strategies[current_strategy]->can_update();
	return is_last_strategy() && can_update == false;
}

bool Patch::has_converged()
{
	return (fitting_error < data_model.optimization_settings.convergence_fitting_threshold  
		    && DOG_error < data_model.optimization_settings.convergence_developablity_threshold);
}

bool Patch::is_last_strategy()
{
	return current_strategy >= constraints_strategies.size() - 1;
}

void Patch::update()
{
	write_log(4) << "  Patch:: update()" << endl;

	//TODO optimization settings are reset in strategy before strategy is advanced --> FIX!!
	has_strategy_changed = false;

	igl::per_vertex_normals(wrapper.V, wrapper.F, wrapper.normals_vertices);
	//igl::per_face_normals(wrapper.V, wrapper.F, wrapper.normals_faces);

	if (!data_model.optimization_settings.should_converge || (data_model.optimization_settings.should_converge && has_converged()))
		renew_position_constraints();
	else
		constraints.update_positions();
}

void Patch::renew_position_constraints()
{
	write_log(4) << "  Patch:: renew_position_constraints()" << endl;

	//if (is_last_strategy())
	//	return;

	if (constraints_strategies.size() > 1
		&& current_strategy < constraints_strategies.size() - 1
		&& constraints_strategies[current_strategy]->can_update() == false)
	{
		//if(data_model.is_pausing_on_adding_patch)
		//	data_model.is_optimizing = false;
		
		advance_constraints_strategy();
	}

	//if (is_last_strategy())
	//	has_strategy_changed = true;

	constraints_strategies[current_strategy]->update();
	//constraints.update_positions();
}

void Patch::advance_constraints_strategy()
{
	write_log(4) << "  Patch:: advance_constraints_strategy()" << endl;

	current_strategy++;
	constraints_strategies[current_strategy]->configure();

	DOG_error = 1e3;
	fitting_error = 1e3;
	is_local_minimum = true;

	has_strategy_changed = true;
}

void Patch::add_wrapper()
{
	create_wrapper();
	align_wrapper();


	meshhelper::add_wrapper(data_model, wrapper);

	igl::per_vertex_normals(wrapper.V, wrapper.F, wrapper.normals_vertices);
	igl::per_face_normals(wrapper.V, wrapper.F, wrapper.normals_faces);
	meshhelper::compute_mesh_properties(wrapper, false);


	/*
	//add wrapper of adequate size to cover all constraints
	const auto min_point = wrapper_curve.colwise().minCoeff();
	const auto max_point = wrapper_curve.colwise().maxCoeff();
	const auto dimensions = max_point - min_point;


	//TODO if curve is very curved (and short??) then make double size patch and scale down??

	const int width = ceil(dimensions(0)) + 1;
	int height = data_model.ruled_width;
	if (height >= width)
		height = 0.8 * width;

	height = height % 2 == 0 ? height + 1 : height;
	//const int height = 5;
	//int height = width * data_model.wrapper_width_scale;
	//height = height < 3 ? 3 : height;



	write_log(4) << "min_point: " << min_point << ", max_point: " << max_point << ", dimensions: " << dimensions << endl;
	write_log(4) << "dimensions(0): " << dimensions(0) << ", quad_width: " << width << endl;
	write_log(4) << "dimensions(1): " << dimensions(1) << ", quad_height: " << height << endl;
	write_log(4) << "wrapper_curve.rows(): " << wrapper_curve.rows() << endl;
	write_log(5) << "wrapper_curve: " << endl << wrapper_curve << endl;

	//meshhelper::add_wrapper(width, height, data_model, wrapper);
	meshhelper::create_wrapper(width, height, data_model, wrapper);
	wrapper.V = wrapper.V.rowwise() + (Eigen::RowVector3d(0, (-height) / 2, 0)); //center wrapper only in y-direction

	align_wrapper();
	*/
}

void Patch::create_wrapper()
{
	/*
	// get wrapper dimesions and resolution

	//add wrapper of adequate size to cover all constraints
	const auto min_point = wrapper_curve.colwise().minCoeff();
	const auto max_point = wrapper_curve.colwise().maxCoeff();
	const auto dimensions = max_point - min_point;

	//if curve is very short, make double-resolution DOG
	const int height_default = data_model.ruled_width;
	const int width_original = ceil(dimensions(0)) + 1;
	//const double resolution_multiplier = width_original < data_model.ruled_width ? 2.0 : 1.0;
	const double resolution_multiplier = 1.0;

	const int width = width_original * resolution_multiplier;
	int height = height_default; // *resolution_multiplier;
	height = height % 2 == 0 ? height + 1 : height;
	
	write_log(4) << "min_point: " << min_point << ", max_point: " << max_point << ", dimensions: " << dimensions << endl;
	write_log(4) << "resolution_multiplier: " << resolution_multiplier << ", quad_width: " << width << ", quad_height: " << height << endl;
	write_log(4) << "wrapper_curve.rows(): " << wrapper_curve.rows() << endl;
	write_log(5) << "wrapper_curve: " << endl << wrapper_curve << endl;
	*/

	// get wrapper dimesions and resolution
	double target_length = curve::compute_length(target_curve);

	//if curve is very short, make double-resolution DOG
	//const int height_default = data_model.ruled_width;
	int height_default = data_model.ruled_width;

	int ruled_index = index_of(geodesic_vertex, data_model.ruled_vertex_indices);
	if (ruled_index > -1)
	{
		height_default = std::ceil(data_model.ruled_developables[ruled_index]->get_width());
		height_default = std::max((float)height_default, data_model.ruled_width / 2);
	}

	const int longer = 2;
	const int width_original = ceil(target_length) + 1+longer;
	//const double resolution_multiplier = width_original < data_model.ruled_width ? 2.0 : 1.0;
	const double resolution_multiplier = 1.0;

	const int width = width_original * resolution_multiplier;
	int height = height_default; // *resolution_multiplier;
	height = height % 2 == 0 ? height + 1 : height;


	resample_constraint_curves(1.0/resolution_multiplier);

	write_log(4) << "target_length: " << target_length << ", resolution_multiplier: " << resolution_multiplier << ", quad_width: " << width << ", quad_height: " << height << endl;
	write_log(4) << "wrapper_curve.rows(): " << wrapper_curve.rows() << endl;
	write_log(5) << "wrapper_curve: " << endl << wrapper_curve << endl;


	//create the wapper mesh
	wrapper.quad_width = width;
	wrapper.quad_height = height;

	get_planar_square_mesh(wrapper.V, wrapper.F, height, width);
	wrapper.F_quad = F_to_Fsqr(wrapper.F);
	// get mesh connectivity data structure (edges, stars, etc..)
	quad_topology(wrapper.V, wrapper.F_quad, wrapper.quad_topology);

	if (!wrapper.F.rows() || !wrapper.F_quad.rows())
		write_log(1) << "error at creating wrapper mesh!" << endl;


	//transform wrapper
	wrapper.V = wrapper.V.rowwise() + (Eigen::RowVector3d(-longer/2, (-height) / 2, 0)); //center wrapper only in y-direction
	wrapper.V /= resolution_multiplier;

	mat2_to_vec(wrapper.V, wrapper.initial_x0);
}


void Patch::align_wrapper()
{
	//align wrapper in 3D (by using *face* normal)

	/*
	int mid = wrapper_curve.rows() % 2 == 0 ? wrapper_curve.rows() / 2 - 1 : wrapper_curve.rows() / 2;
	Eigen::RowVectorXd normal_wrapper;
	Eigen::RowVectorXd normal_target;
	get_face_normal(wrapper_curve.row(mid), wrapper.V, wrapper.F, normal_wrapper);
	get_face_normal(target_curve.row(mid), data_model.target.V, data_model.target.F, normal_target);

	double normal_scale = 1.0;
	int number_added_points = 1;

	Eigen::MatrixXd to_transform(wrapper_curve.rows() + number_added_points, wrapper_curve.cols());
	to_transform << (normal_wrapper * normal_scale), wrapper_curve;
	//write_log(0) << "    add patch:: wrapper_curve = " << endl << wrapper_curve << endl << endl;
	write_log(0) << "    add patch:: normal_wrapper = " << endl << normal_wrapper << endl << endl;
	//write_log(0) << "    add patch:: to_transform = " << endl << to_transform << endl << endl << endl;

	Eigen::MatrixXd transform_target(target_curve.rows() + number_added_points, target_curve.cols());
	transform_target << (normal_target * normal_scale), target_curve;
	//write_log(0) << "    add patch:: target_curve = " << endl << target_curve << endl << endl;
	write_log(0) << "    add patch:: normal_target = " << endl << normal_target << endl << endl;
	//write_log(0) << "    add patch:: transform_target = " << endl << transform_target << endl << endl << endl;
	*/

	int mid = wrapper_curve.rows() % 2 == 0 ? wrapper_curve.rows() / 2 - 1 : wrapper_curve.rows() / 2;
	Eigen::MatrixXd frame_wrapper = get_surface_curve_frame_at(mid, wrapper_curve, wrapper.V, wrapper.F); //TNB
	Eigen::MatrixXd frame_target = get_surface_curve_frame_at(mid, target_curve, data_model.target.V, data_model.target.F);

	write_log(0) << "frame_wrapper [T N B]: " << linebreak << frame_wrapper << endl;
	frame_wrapper = frame_wrapper.colwise() + wrapper_curve.row(mid).transpose();
	write_log(0) << "frame_wrapper [T N B]: " << linebreak << frame_wrapper << linebreak << endl;

	write_log(0) << "frame_target [T N B]: "  << linebreak << frame_target  << endl;
	frame_target = frame_target.colwise() + target_curve.row(mid).transpose();
	write_log(0) << "frame_target [T N B]: "  << linebreak << frame_target  << linebreak << endl;



	double normal_scale = 1.0;
	int number_added_points = 1;

	Eigen::MatrixXd to_transform(wrapper_curve.rows() + number_added_points, wrapper_curve.cols());
	to_transform << (frame_wrapper.col(1).transpose() * normal_scale), wrapper_curve;

	Eigen::MatrixXd transform_target(target_curve.rows() + number_added_points, target_curve.cols());
	transform_target << (frame_target.col(1).transpose() * normal_scale), target_curve;

	//int number_added_points = 2;

	//Eigen::MatrixXd to_transform(wrapper_curve.rows() + number_added_points, wrapper_curve.cols());
	//to_transform << (frame_wrapper.col(1).transpose() * normal_scale), (frame_wrapper.col(2).transpose() * normal_scale), wrapper_curve;
	//write_log(0) << "    add patch:: normal_wrapper = " << endl << frame_wrapper.col(1).transpose() << endl << endl;

	//Eigen::MatrixXd transform_target(target_curve.rows() + number_added_points, target_curve.cols());
	//transform_target << (frame_target.col(1).transpose() * normal_scale), (frame_target.col(2).transpose() * normal_scale), target_curve;
	//write_log(0) << "    add patch:: normal_target = " << endl << frame_target.col(1).transpose() << endl << endl;


	Eigen::RowVector3d t;
	Eigen::Matrix3d R;
	//igl::procrustes(wrapper_curve, target_curve, R, t);
	igl::procrustes(to_transform, transform_target, false, false, R, t);

	wrapper.V = (wrapper.V * R).rowwise() + t;
	wrapper_curve = (wrapper_curve * R).rowwise() + t;
}

void Patch::resample_constraint_curves(double segment_length)
{
	Eigen::MatrixXd original_target_curve = target_curve;
	double origianl_target_length = curve::compute_length(original_target_curve);

	double resampled_target_length;
	curve::resample_uniformly(original_target_curve, segment_length, target_curve, resampled_target_length);

	double target_off = abs(origianl_target_length - resampled_target_length);
	if (target_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;


	double resampled_wrapper_length;
	curve::resample_matching_curve(target_curve, 0, sample_direction, wrapper_curve, resampled_wrapper_length);

	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
	if (wrapper_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;
}
