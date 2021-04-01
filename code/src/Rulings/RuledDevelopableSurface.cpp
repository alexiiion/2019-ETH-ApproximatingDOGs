#include "RuledDevelopableSurface.h"

#include <igl/Timer.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

#include "AnalyticRulings.h"
#include "CompoundRulings.h"
#include "DiscreteRulings.h"

#include "GlobalSettings.h"
#include "CurveHelper.h"
#include "Logger.h"


using namespace std;

int sign(double value)
{
	if (value > 0) return 1;
	if (value < 0) return -1;
	return 0;
}

RuledDevelopableSurface::RuledDevelopableSurface(Mesh& target, Eigen::MatrixXd& geodesic, const int at_vertex, const RulingsType& rulings_type) : geodesic(geodesic), vertex(at_vertex), rulings_type(rulings_type)
{
	if (rulings_type == RulingsType::Analytic)
		rulings_controller = new AnalyticRulings(target, geodesic);
	else if (rulings_type == RulingsType::Discrete)
		rulings_controller = new DiscreteRulings(target, geodesic);
	else if(rulings_type == RulingsType::Compound)
		rulings_controller = new CompoundRulings(target, geodesic);
}

RuledDevelopableSurface::~RuledDevelopableSurface()
{
	delete rulings_controller;
}

Mesh RuledDevelopableSurface::create(const double width)
{
	this->width = width;
	GlobalSettings::min_ruling_length = width * 0.5;

	rulings_controller->compute();
	developable = create_mesh(rulings_controller->ruling_directions, rulings_controller->ruling_lengths, width);

	return developable;
}

//Mesh RuledDevelopableSurface::create(const double width)
//{
//	//igl::Timer timer;
//	//double t = timer.getElapsedTime();
//
//	rulings_controller->compute_curve();
//	//write_log(4) << "    compute curve in " << timer.getElapsedTime() - t << endl;
//	//t = timer.getElapsedTime();
//
//	rulings = rulings_controller->compute_ruling_directions();
//	//write_log(4) << "    compute ruling directions in " << timer.getElapsedTime() - t << endl;
//	//t = timer.getElapsedTime();
//
//	int max_k_index;
//	rulings_controller->curvature.maxCoeff(&max_k_index);
//
//	direction_normalizer = normalize_rulings_side(rulings, max_k_index);
//	ruling_lengths = compute_ruling_lengths(rulings);
//	//write_log(4) << "    compute ruling lengths in " << timer.getElapsedTime() - t << endl;
//	//t = timer.getElapsedTime();
//
//
//	rulings_controller->compute();
//	developable = create_mesh(rulings_controller->ruling_directions, rulings_controller->ruling_lengths, width);
//	//write_log(4) << "    create developable mesh in " << timer.getElapsedTime() - t << endl;
//	//t = timer.getElapsedTime();
//
//	return developable;
//}

Mesh RuledDevelopableSurface::create_mesh(const Eigen::MatrixXd& rulings, const Eigen::VectorXd& ruling_lengths, const double width)
{
	//compute surface: V & F 
	const int number_samples = rulings.rows();
	int number_vertices = number_samples * 2;

	const double max_width = width * 0.5;
	const double min_width = width * 0.125; //GlobalSettings::min_ruling_length; //width / 8.0; 	//2.0;
	int added_vertices = 0;

	write_log(5) << "RuledDevelopableSurface::create_mesh() ruling_lengths: " << linebreak << ruling_lengths << endl;
	write_log(5) << "RuledDevelopableSurface::create_mesh() rulings: " << linebreak << rulings << endl;

	Eigen::MatrixXd V(number_vertices, 3);

	for (int i = 0; i < number_samples; i++)
	{
		double ruling_length = ruling_lengths(i);
		if (abs(ruling_length) < min_width)
			continue;

		//int sign_multiplier = sign(ruling_length);
		double clamped_length = abs(ruling_length) > max_width ? max_width * sign(ruling_length) : ruling_length;

		double positive_length;
		double negative_length;

		//add ruling lengths for both sides of the surface
		if (clamped_length < 0) //intersection is negative
		{
			negative_length = clamped_length;
			positive_length = max_width;
		}
		else if (clamped_length > 0) //intersection is positive
		{
			positive_length = clamped_length;
			negative_length = max_width * -1;
		}

		V.row(added_vertices * 2)     = rulings_controller->sampled_curve.row(i) + negative_length * rulings.row(i);
		V.row(added_vertices * 2 + 1) = rulings_controller->sampled_curve.row(i) + positive_length * rulings.row(i);
		added_vertices++;
	}

	if (added_vertices < 3)
		return Mesh();


	V.conservativeResize(added_vertices*2, Eigen::NoChange);


	int number_faces = (added_vertices*2) - 2;
	Eigen::MatrixXi F(number_faces, 3);

	for (int fi = 0; fi < number_faces; fi++)
	{
		if (fi % 2 == 0)
		{
			F(fi, 0) = fi;
			F(fi, 1) = fi + 1;
			F(fi, 2) = fi + 2;
		}
		else
		{
			F(fi, 0) = fi + 1;
			F(fi, 1) = fi;
			F(fi, 2) = fi + 2;
		}
	}
	//write_log(4) << "ruled surface F: " << std::endl << F << std::endl << std::endl;

	Mesh surface;
	surface.V = V;
	surface.F = F;

	igl::per_vertex_normals(surface.V, surface.F, surface.normals_vertices);
	igl::per_face_normals(surface.V, surface.F, surface.normals_faces);

	return surface;
}

//Eigen::Vector3d RuledDevelopableSurface::normalize_rulings_side(Eigen::MatrixXd& rulings, int max_k_index)
//{
//	//normalize ruling directions to a stable ruling
//	direction_normalizer = rulings.row(max_k_index);
//
//	for (int i = 0; i < rulings.rows(); i++)
//	{
//		Eigen::RowVector3d ruling_current = rulings.row(i);
//
//		if (ruling_current.dot(direction_normalizer) < 0)
//			rulings.row(i) = ruling_current * -1;
//	}
//
//	return direction_normalizer;
//}

//Eigen::VectorXd RuledDevelopableSurface::compute_ruling_lengths(Eigen::MatrixXd& rulings)
//{
//	const int number_rulings = rulings.rows();
//	Eigen::VectorXd ruling_lengths(number_rulings);
//	ruling_lengths.setZero();
//
//	const double curve_length = curve::compute_length(rulings_controller->sampled_curve);
//	const double segment = curve_length / rulings_controller->sampled_curve.rows();
//
//	////int max_k_index;
//	////k_div_torsion.cwiseAbs().maxCoeff(&max_k_index);
//	////direction_normalizer = rulings.row(max_k_index);
//
//	//Eigen::RowVector3d ruling_current = rulings.row(0);
//	//
//	//for (int i = 0; i < number_rulings - 1; i++)
//	//{
//	//	Eigen::RowVector3d tangent = rulings_controller->tangents.row(i);
//	//	Eigen::RowVector3d ruling_next = rulings.row(i + 1);
//	//
//	//	double length = compute_ruling_length(tangent, ruling_current, ruling_next, direction_normalizer, segment);
//	//	ruling_lengths(i) = length;
//	//	ruling_current = ruling_next;
//	//}
//	//
//	////add last ruling length, compute between last two rulings, but from the opposite direction
//	//Eigen::RowVector3d tangent = rulings_controller->tangents.row(number_rulings - 1) * -1;
//	//ruling_current = rulings.row(number_rulings - 1);
//	//Eigen::RowVector3d ruling_next = rulings.row(number_rulings - 2);
//	//
//	//double length = compute_ruling_length(tangent, ruling_current, ruling_next, direction_normalizer, segment);
//	//ruling_lengths(number_rulings - 1) = length;
//
//
//	/*
//	int previous_index = 0;
//	for (int i = 1; i < number_rulings; i++)
//	{
//		if (rulings.row(i).isZero())
//			continue;
//
//		Eigen::RowVector3d tangent = rulings_controller->tangents.row(previous_index);
//		Eigen::RowVector3d ruling = rulings.row(i);
//		Eigen::RowVector3d ruling_previous = rulings.row(previous_index);
//
//		double length = compute_ruling_length(tangent, ruling_previous, ruling, direction_normalizer, segment*(i-previous_index));
//		ruling_lengths(previous_index) = length;
//		write_log(0) << "    ruling_lengths:: between " << previous_index << " and " << i << " --> length = " << length << endl;
//
//		if(i < number_rulings-1) // dont update last b/c i need to go backwards
//			previous_index = i;
//	}
//
//	//add last ruling length, compute between last two rulings, but from the opposite direction
//	Eigen::RowVector3d tangent_last = rulings_controller->tangents.row(number_rulings-1) * -1;
//	Eigen::RowVector3d ruling_last = rulings.row(number_rulings - 1);
//	Eigen::RowVector3d ruling_previous = rulings.row(previous_index);
//
//	double length = compute_ruling_length(tangent_last, ruling_last, ruling_previous, direction_normalizer, segment * (number_rulings - 1 - previous_index));
//	ruling_lengths(number_rulings - 1) = length;
//	write_log(0) << "    ruling_lengths:: LAST between " << number_rulings - 1 << " and " << previous_index << " --> length = " << length << endl;
//	*/
//
//
//	//compute lengths forward
//	int current_index = 0;
//	for (int forward_index = 1; forward_index < number_rulings; forward_index++)
//	{
//		if (rulings.row(forward_index).isZero())
//			continue;
//	
//		Eigen::RowVector3d tangent = rulings_controller->tangents.row(current_index);
//		double length = compute_ruling_length(current_index, forward_index, rulings, tangent, segment * (forward_index - current_index));
//		ruling_lengths(current_index) = length;
//		
//		//write_log(0) << "    ruling_lengths forward:: between " << current_index << " and " << forward_index << " --> length = " << length << endl;
//		current_index = forward_index;
//	}
//	//write_log(0) << "ruling_lengths = " << linebreak << ruling_lengths << endl << endl;
//
//	//compute lengths backward
//	current_index = number_rulings-1;
//	for (int backward_index = number_rulings - 2; backward_index >= 0; backward_index--)
//	{
//		if (rulings.row(backward_index).isZero())
//			continue;
//
//		Eigen::RowVector3d tangent = rulings_controller->tangents.row(current_index)*-1;
//		double length = compute_ruling_length(current_index, backward_index, rulings, tangent, segment * (current_index - backward_index));
//		if(abs(length) < abs(ruling_lengths(current_index)))
//			ruling_lengths(current_index) = length;
//
//		//write_log(0) << "    ruling_lengths backward:: between " << current_index << " and " << backward_index << " --> length = " << length << endl;
//		current_index = backward_index;
//	}
//	//write_log(0) << "ruling_lengths = " << linebreak << ruling_lengths << endl << endl;
//
//	////TODO remove!
//	//ruling_lengths(number_rulings - 1) = -5;
//
//
//	/*
//	int previous_index = 0;
//	for (int i = 1; i < number_rulings; i++)
//	{
//		if (rulings.row(i).isZero())
//			continue;
//
//
//
//		Eigen::RowVector3d tangent = rulings_controller->tangents.row(previous_index);
//		Eigen::RowVector3d ruling = rulings.row(i);
//		Eigen::RowVector3d ruling_previous = rulings.row(previous_index);
//
//		double length = compute_ruling_length(tangent, ruling_previous, ruling, direction_normalizer, segment*(i-previous_index));
//		ruling_lengths(previous_index) = length;
//		write_log(0) << "    ruling_lengths:: between " << previous_index << " and " << i << " --> length = " << length << endl;
//
//		if(i < number_rulings-1) // dont update last b/c i need to go backwards
//			previous_index = i;
//	}
//
//	//add last ruling length, compute between last two rulings, but from the opposite direction
//	Eigen::RowVector3d tangent_last = rulings_controller->tangents.row(number_rulings-1) * -1;
//	Eigen::RowVector3d ruling_last = rulings.row(number_rulings - 1);
//	Eigen::RowVector3d ruling_previous = rulings.row(previous_index);
//
//	double length = compute_ruling_length(tangent_last, ruling_last, ruling_previous, direction_normalizer, segment * (number_rulings - 1 - previous_index));
//	ruling_lengths(number_rulings - 1) = length;
//	write_log(0) << "    ruling_lengths:: LAST between " << number_rulings - 1 << " and " << previous_index << " --> length = " << length << endl;
//	*/
//
//	return ruling_lengths;
//}

//double RuledDevelopableSurface::compute_ruling_length_bidirectional(const int index, const int previous, const int next, const Eigen::MatrixXd& rulings, const double step)
//{
//	double length_forward = compute_ruling_length(rulings_controller->tangents.row(index), rulings.row(index), rulings.row(next), direction_normalizer, step * (next - index));
//	double length_backward= compute_ruling_length(rulings_controller->tangents.row(index)*-1, rulings.row(index), rulings.row(previous), direction_normalizer, step * (index - previous));
//	
//	double length = abs(length_backward) < abs(length_forward) ? length_backward : length_forward;
//	if (!isfinite(length))
//		return 0;
//
//	return length;
//}

//double RuledDevelopableSurface::compute_ruling_length(const int index, const int other, const Eigen::MatrixXd& rulings, const Eigen::RowVector3d& tangent, const double step)
//{
//	double length = compute_ruling_length(tangent, rulings.row(index), rulings.row(other), direction_normalizer, step * abs(other - index));
//	return length;
//}
//
//double RuledDevelopableSurface::compute_ruling_length(const Eigen::RowVector3d& tangent, const Eigen::RowVector3d& ruling_current, const Eigen::RowVector3d& ruling_next, const Eigen::RowVector3d& ruling_direction, const double step)
//{
//	double angle_current = acos(ruling_current.dot(tangent));
//
//	//if the ruling is to close to the tangent, skip it 
//	if (angle_current < 0.1 || M_PI - angle_current < 0.1)
//		return 0;
//	
//	//compute length until intersection with next ruling (singularity)
//	double angle_next = acos(ruling_next.dot(tangent));
//	double length = (step * sin(angle_next)) / sin(angle_next - angle_current);
//
//	if (!isfinite(length))
//		return 0;
//
//	//if (abs(angle_current - angle_next) < 0.0001)
//	//	length = 0;
//
//	return length;
//}

Eigen::MatrixXd RuledDevelopableSurface::generator_curve()
{
	return rulings_controller->sampled_curve;
}
int RuledDevelopableSurface::get_width()
{
	return width;
}

void RuledDevelopableSurface::print_mathematica_data()
{
	std::cout << "(*vertex: " << vertex << "*)" << std::endl;

	rulings_controller->print_mathematica_data();

	std::cout << "V="  << to_mathematica(developable.V) << std::endl;
	std::cout << "nV=" << developable.V.rows() << ";" << std::endl;
}
