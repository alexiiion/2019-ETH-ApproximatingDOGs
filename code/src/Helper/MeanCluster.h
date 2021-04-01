#pragma once

#include <vector>

struct MeanCluster
{
public:

	MeanCluster()
		: window_size(5), outlier_factor(2.0), use_adaptive_bounds_update(true)
	{ /* empty on purpose*/ };

	MeanCluster(const int window_size, const double outlier_factor)
		: window_size(window_size), outlier_factor(outlier_factor)
	{ /* empty on purpose*/ };


	/* Methods */

	bool is_value_within(double value);
	void reset();


	/* Members */

	int window_size;
	double outlier_factor;

	double absolute_threshold = 1e12;
	double bounds_min_absolute = 0.01;
	bool use_upper_bound_only = false;
	bool use_adaptive_bounds_update = false;

private:
	/* Members */

	std::vector<double> data;

	double sum;
	double sum_squared;
	
	double mean;
	double variance;
	double stdev;
	
	double bounds;


	/* Methods */

	bool is_within_bounds(double value);
	void update_statistics();
	void reset_statistics();
};
