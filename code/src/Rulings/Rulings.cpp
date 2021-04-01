#include "Rulings.h"

#include "CurveHelper.h"
#include "Logger.h"


void Rulings::compute()
{
	//write_log(0) << "Rulings::compute()" << std::endl;

	compute_curve();
	ruling_directions = compute_ruling_directions();
	direction_normalizer = normalize_rulings_side(ruling_directions);
	ruling_lengths = compute_ruling_lengths(ruling_directions);
}

Eigen::Vector3d Rulings::normalize_rulings_side(Eigen::MatrixXd& rulings)
{
	//write_log(0) << "Rulings::normalize_rulings_side()" << std::endl;

	//normalize ruling directions to a stable ruling
	int max_k_index;
	curvature.maxCoeff(&max_k_index);
	direction_normalizer = rulings.row(max_k_index);

	for (int i = 0; i < rulings.rows(); i++)
	{
		Eigen::RowVector3d ruling_current = rulings.row(i);

		if (ruling_current.dot(direction_normalizer) < 0)
			rulings.row(i) = ruling_current * -1;
	}

	return direction_normalizer;
}

Eigen::VectorXd Rulings::compute_ruling_lengths(Eigen::MatrixXd& rulings)
{
	//write_log(4) << "Rulings::compute_ruling_lengths()" << std::endl;
	//write_log(5) << "rulings: " << linebreak << rulings << std::endl;



	const int number_rulings = rulings.rows();
	Eigen::VectorXd ruling_lengths(number_rulings);
	ruling_lengths.setZero();

	const double curve_length = curve::compute_length(sampled_curve);
	const double segment = curve_length / sampled_curve.rows();

	//int max_k_index;
	//k_div_torsion.cwiseAbs().maxCoeff(&max_k_index);
	//direction_normalizer = rulings.row(max_k_index);
	
	Eigen::RowVector3d ruling_current = rulings.row(0);
	
	for (int i = 0; i < number_rulings - 1; i++)
	{
		Eigen::RowVector3d tangent = tangents.row(i);
		Eigen::RowVector3d ruling_next = rulings.row(i + 1);
	
		double length = compute_ruling_length(tangent, ruling_current, ruling_next, direction_normalizer, segment);
		ruling_lengths(i) = length;
		ruling_current = ruling_next;
	}
	
	//add last ruling length, compute between last two rulings, but from the opposite direction
	Eigen::RowVector3d tangent = tangents.row(number_rulings - 1) * -1;
	ruling_current = rulings.row(number_rulings - 1);
	Eigen::RowVector3d ruling_next = rulings.row(number_rulings - 2);
	
	double length = compute_ruling_length(tangent, ruling_current, ruling_next, direction_normalizer, segment);
	ruling_lengths(number_rulings - 1) = length;


	/*
	//compute lengths forward
	int current_index = 0;
	for (int forward_index = 1; forward_index < number_rulings; forward_index++)
	{
		if (rulings.row(forward_index).isZero())
			continue;

		Eigen::RowVector3d tangent = tangents.row(current_index);
		double length = compute_ruling_length(current_index, forward_index, rulings, tangent, segment * (forward_index - current_index));
		ruling_lengths(current_index) = length;

		//write_log(0) << "    ruling_lengths forward:: between " << current_index << " and " << forward_index << " --> length = " << length << endl;
		current_index = forward_index;
	}
	//write_log(0) << "ruling_lengths = " << linebreak << ruling_lengths << endl << endl;

	//compute lengths backward
	current_index = number_rulings - 1;
	for (int backward_index = number_rulings - 2; backward_index >= 0; backward_index--)
	{
		if (rulings.row(backward_index).isZero())
			continue;

		Eigen::RowVector3d tangent = tangents.row(current_index) * -1;
		double length = compute_ruling_length(current_index, backward_index, rulings, tangent, segment * (current_index - backward_index));
		if (abs(length) < abs(ruling_lengths(current_index)))
			ruling_lengths(current_index) = length;

		//write_log(0) << "    ruling_lengths backward:: between " << current_index << " and " << backward_index << " --> length = " << length << endl;
		current_index = backward_index;
	}
	//write_log(0) << "ruling_lengths = " << linebreak << ruling_lengths << endl << endl;
	*/
	return ruling_lengths;
}


double Rulings::compute_ruling_length(const int index, const int other, const Eigen::MatrixXd& rulings, const Eigen::RowVector3d& tangent, const double step)
{
	double length = compute_ruling_length(tangent, rulings.row(index), rulings.row(other), direction_normalizer, step * abs(other - index));
	return length;
}

double Rulings::compute_ruling_length(const Eigen::RowVector3d& tangent, const Eigen::RowVector3d& ruling_current, const Eigen::RowVector3d& ruling_next, const Eigen::RowVector3d& ruling_direction, const double step)
{
	//write_log(0) << "ruling_current: " << ruling_current << std::endl;
	//write_log(0) << "ruling_next:    " << ruling_next << std::endl;
	//write_log(0) << "tangent:        " << tangent << std::endl;

	double angle_current = acos(ruling_current.dot(tangent));

	//if the ruling is to close to the tangent, skip it 
	if (angle_current < 0.1 || M_PI - angle_current < 0.1)
		return 0;

	//compute length until intersection with next ruling (singularity)
	double angle_next = acos(ruling_next.dot(tangent));
	double length = (step * sin(angle_next)) / sin(angle_next - angle_current);

	if (isinf(length))
		return 1e6;
	if (!isfinite(length))
		return 0;

	//if (abs(angle_current - angle_next) < 0.0001)
	//	length = 0;

	return length;
}

void Rulings::print_mathematica_data()
{
	std::cout << "geodesic=" << to_mathematica(geodesic) << std::endl;
	std::cout << "ng=" << geodesic.rows() << ";" << std::endl << std::endl;
	std::cout << "curve=" << to_mathematica(sampled_curve) << std::endl;
	std::cout << "n=" << sampled_curve.rows() << ";" << std::endl << std::endl;

	std::cout << "tangents=" << to_mathematica(tangents) << std::endl;
	std::cout << "principalNormals=" << to_mathematica(principal_normals) << std::endl;
	std::cout << "binormals=" << to_mathematica(binormals) << std::endl;

	std::cout << "curvature=" << to_mathematica(curvature) << std::endl;
	std::cout << "torsion=" << to_mathematica(torsion) << std::endl;
	std::cout << "rulings=" << to_mathematica(ruling_directions) << std::endl;
	std::cout << "rulingsLenghts=" << to_mathematica(ruling_lengths) << std::endl;
}
