#pragma once

#include "ConstraintsStrategy.h"

class CoverTargetConstraints : public ConstraintsStrategy
{

public:
	CoverTargetConstraints(Patch* patch, PositionConstraints& constraints) : ConstraintsStrategy(patch, constraints) {};

	virtual void configure() override;

	virtual bool can_update() override;
	virtual void update() override;

private:
	bool do_update;

	int current_iteration;
	int max_iterations;

	double previous_error = 0.0;
	double previous_derivative = 1000000000.0;
	double min_derivative = 1.0;

	int minimum_iterations = 0;
	int minimum_window = 0;

	OptimizationSettings previous_weights;
	
	double current_error_squared();
};