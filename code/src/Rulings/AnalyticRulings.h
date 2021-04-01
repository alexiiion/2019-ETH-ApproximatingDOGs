#pragma once

#include <Eigen/Core>

#include "Rulings.h"
#include "Spline.h"

class AnalyticRulings : public Rulings
{
public:
	AnalyticRulings(Mesh& target, Eigen::MatrixXd& geodesic);
	~AnalyticRulings();

	virtual void print_mathematica_data() override;

protected:
	double smoothness;
	double sample_step;
	int degree;

	Spline* spline;
	Eigen::MatrixXd derivative_1;
	Eigen::MatrixXd derivative_2;
	Eigen::MatrixXd derivative_3;

	virtual void compute_curve() override;
	virtual Eigen::MatrixXd compute_ruling_directions() override;

	Eigen::MatrixXd cross_rowwise(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2);
};
