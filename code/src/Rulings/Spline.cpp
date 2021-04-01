#include "Spline.h"

#include <iostream>


const double cut_edges = 0.0;

Spline::Spline(const Eigen::MatrixXd& geodesic) : geodesic(geodesic)
{
	this->smoothness = 3.0;
	this->degree = 5;

	create(smoothness, degree);
}

Spline::Spline(const Eigen::MatrixXd& geodesic, const double smoothness, const int degree) : geodesic(geodesic), smoothness(smoothness), degree(degree)
{
	this->smoothness = smoothness;
	this->degree = degree;

	create(smoothness, degree);
}

Spline::~Spline()
{
	delete x_spline;
	delete y_spline;
	delete z_spline;
}

void Spline::create(const double& smoothness, const int& degree)
{
	if (geodesic.rows() < 1)
		return;

	calculate_geodesic_length();

	std::vector<double> x_samples;
	std::vector<double> y_samples;
	std::vector<double> z_samples;
	std::vector<double> t_samples;

	double t = 0.0;
	double current_length = 0.0;

	for (int i = 0; i < geodesic.rows(); i++)
	{
		if (i > 0)
		{
			current_length += (geodesic.row(i) - geodesic.row(i - 1)).norm();
			t = current_length / geodesic_length;
		}

		t_samples.push_back(t);
		x_samples.push_back(geodesic(i, 0));
		y_samples.push_back(geodesic(i, 1));
		z_samples.push_back(geodesic(i, 2));
	}

	x_spline = new fitpackpp::BSplineCurve(t_samples, x_samples, degree, smoothness);
	y_spline = new fitpackpp::BSplineCurve(t_samples, y_samples, degree, smoothness);
	z_spline = new fitpackpp::BSplineCurve(t_samples, z_samples, degree, smoothness);

	is_initialized = true;
}

Eigen::MatrixXd Spline::sample_curve(const double sampling_step)
{
	if (!is_initialized)
		return Eigen::MatrixXd();

	const double query_start = cut_edges;
	const double query_end = 1.0 - cut_edges;

	const int number_samples = (query_end - query_start) / sampling_step;
	Eigen::MatrixXd spline(number_samples + 1, 3);

	//std::cout << "spline: query_start=" << query_start << ", query_end=" << query_end << ", number_samples=" << number_samples << ", reserved rows=" << spline.rows() << std::endl;

	for (int i = 0; i <= number_samples; i++)
	{
		double query_t = i * sampling_step;
		//std::cout << "  spline(t = " << query_t << ")" << std::endl;

		double query_x = x_spline->eval(query_t);
		double query_y = y_spline->eval(query_t);
		double query_z = z_spline->eval(query_t);

		spline(i, 0) = query_x;
		spline(i, 1) = query_y;
		spline(i, 2) = query_z;
	}

	//std::cout << "spline(t=1.0) = " << x_spline->eval(1.0) << ", " << y_spline->eval(1.0) << ", " << z_spline->eval(1.0) << std::endl;

	return spline;
}

Eigen::MatrixXd Spline::sample_derivative(const int order, const double sampling_step)
{
	if (!is_initialized)
		return Eigen::MatrixXd();

	const double query_start = cut_edges;
	const double query_end = 1.0 - cut_edges;

	const int number_samples = (query_end - query_start) / sampling_step;
	Eigen::MatrixXd derivative(number_samples + 1, 3);

	//write_log(4) << "derivative " << order << ": query_start=" << query_start << ", query_end=" << query_end << ", number_samples=" << number_samples << ", reserved rows=" << derivative.rows() << std::endl;

	for (int i = 0; i <= number_samples; i++)
	{
		double query_t = i * sampling_step;
		//std::cout << "  derivative of spline(t = " << query_t << ")" << std::endl;

		double query_x = x_spline->der(query_t, order);
		double query_y = y_spline->der(query_t, order);
		double query_z = z_spline->der(query_t, order);

		derivative(i, 0) = query_x;
		derivative(i, 1) = query_y;
		derivative(i, 2) = query_z;
	}

	//std::cout << "derivative " << order << " of spline(t=1.0) = " << x_spline->der(1.0, order) << ", " << y_spline->der(1.0, order) << ", " << z_spline->der(1.0, order) << std::endl;

	return derivative;
}

void Spline::calculate_geodesic_length()
{
	double t = 0.0;

	for (int i = 1; i < geodesic.rows(); i++)
		t += (geodesic.row(i) - geodesic.row(i - 1)).norm();

	this->geodesic_length = t;
}