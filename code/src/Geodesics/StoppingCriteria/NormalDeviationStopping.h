#pragma once

#include <vector>
#include "GeodesicsStopping.h"
#include "MeanCluster.h"

struct NormalDeviationStopping : public GeodesicsStopping
{

public:
	NormalDeviationStopping();

	virtual bool should_stop() override;
	virtual void update(const WalkedData& data, const Mesh& mesh) override;
	virtual void reset() override;

	int window_size = 5;
	double outlier_factor = 3.0;

	double absolute_threshold = 1e12;
	double bounds_min_absolute = 0.01;
	bool use_upper_bound_only = true;
	bool use_adaptive_bounds_update = true;

	double winding_factor = 1.25;

private:
	MeanCluster cluster;

	double last_deviation = 0.0;
	double angle_sum = 0.0;
};