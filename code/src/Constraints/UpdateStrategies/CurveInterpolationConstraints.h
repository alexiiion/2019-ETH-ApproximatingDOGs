#pragma once

#include "ConstraintsStrategy.h"

class CurveInterpolationConstraints : public ConstraintsStrategy
{
public:

	CurveInterpolationConstraints(Patch* patch, PositionConstraints& constraints) : ConstraintsStrategy(patch, constraints) {};

	virtual void configure() override;

	virtual bool can_update() override;
	virtual void update() override;


private:
	double timestep_increment = 0.02;
	double current_timestep = 0.0;

	bool previous_should_converge;

	void interpolate_geodesic_constraints(Eigen::MatrixXd& out_constraints_from);

};