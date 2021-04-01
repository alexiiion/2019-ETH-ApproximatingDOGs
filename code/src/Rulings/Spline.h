#pragma once

#include <Eigen/Core>
#include "BSplineCurve.h"

class Spline
{
public:

	Spline(const Eigen::MatrixXd& geodesic);
	Spline(const Eigen::MatrixXd& geodesic, const double smoothness, const int degree);
	//Spline(Eigen::MatrixXd& geodesic) : geodesic(geodesic) {};
	//Spline(Eigen::MatrixXd& geodesic, double& smoothness, int& degree) : geodesic(geodesic), smoothness(smoothness), degree(degree) {};
	~Spline();

	Eigen::MatrixXd sample_curve(const double sampling_step = 0.1);
	Eigen::MatrixXd sample_derivative(const int order, const double sampling_step = 0.1);


private:
	bool is_initialized = false;

	int degree;
	double smoothness;
	double sampling_step;

	const Eigen::MatrixXd& geodesic;
	double geodesic_length;

	fitpackpp::BSplineCurve* x_spline;
	fitpackpp::BSplineCurve* y_spline;
	fitpackpp::BSplineCurve* z_spline;

	void create(const double& smoothness, const int& degree);
	void calculate_geodesic_length();

};
