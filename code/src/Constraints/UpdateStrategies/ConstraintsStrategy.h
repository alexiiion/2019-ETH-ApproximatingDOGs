#pragma once

#include "PositionConstraintsModel.h"
#include "PatchModel.h"

class ConstraintsStrategy
{
public:
	ConstraintsStrategy(Patch* patch, PositionConstraints& constraints) : constraints(constraints)
	{
		this->patch = patch;
	};

	virtual void configure() = 0;

	virtual bool can_update() = 0;
	virtual void update() = 0;

protected:
	Patch* patch;
	PositionConstraints& constraints;
};

