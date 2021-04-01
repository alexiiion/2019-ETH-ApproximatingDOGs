#pragma once

#include "ConstraintsStrategy.h"

class GrowTargetConstraints : public ConstraintsStrategy
{

public:
	GrowTargetConstraints(Patch* patch, PositionConstraints& constraints) : ConstraintsStrategy(patch, constraints) {};

	virtual void configure() override;

	virtual bool can_update() override;
	virtual void update() override;

private:
	std::vector<int> constrainted_vertices;
	std::vector<int> last_added_vertices;
	
	bool do_update;
	bool is_first_update;

	float previous_weight_fitting_energy;
	bool previous_should_converge;


	int current_grow_iteration;
	int max_number_grow_iterations;

};