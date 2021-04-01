#pragma once

#include <Eigen/Core>

#include "GeodesicsStopping.h"
#include "WalkedData.h"
class Mesh;

struct LengthGeodesicsStopping : public GeodesicsStopping
{
public:
	LengthGeodesicsStopping() {};

	virtual bool should_stop() override
	{
		return current_length >= max_length;
	};

	virtual void update(const WalkedData& data, const Mesh& mesh) override
	{
		current_length = data.length;
	};

	virtual void reset() override
	{
		current_length = 0.0;
	};

	double max_length = 1e3;


private:
	double current_length = 0.0;

};
