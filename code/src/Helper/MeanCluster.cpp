#include "MeanCluster.h"

#include "Logger.h"

bool MeanCluster::is_value_within(double value)
{
	//check for absolute threshold
	if (data.size() < 2 && value >= absolute_threshold)
	{
		write_log(6) << "value >= absolute_threshold" << std::endl;
		return false;
	}

	if (data.size() < 2)
	{
		data.push_back(value);
		update_statistics();

		return true;
	}

	bool is_within = is_within_bounds(value);
	if (!is_within)
	{
		reset();
		return false;
	}

	if (data.size() >= window_size)
		data.erase(data.begin());

	data.push_back(value);
	update_statistics();

	write_log(6) << "current value: " << value << ", bounds: " << mean - std::abs(bounds) << "  -  " << mean + std::abs(bounds) << ", is in? " << std::boolalpha << is_within << std::endl;
	write_log(6) << "   (mean: " << mean << ", varicance: " << variance << ", stdev: " << stdev << ")" << std::endl;
	
	return is_within;
}

void MeanCluster::update_statistics()
{
	reset_statistics();

	for (double value : data)
	{
		sum += value;
		sum_squared += value * value;
	}

	mean = sum / data.size();
	variance = (sum_squared / data.size()) - (mean * mean);
	stdev = std::sqrt(variance);

	if (use_adaptive_bounds_update)
	{
		double factor = std::pow(outlier_factor, (window_size + 1) - data.size());
		bounds = stdev * factor;
	}
	else
	{
		bounds = stdev * outlier_factor;
	}

	if (bounds < bounds_min_absolute)
		bounds = bounds_min_absolute;
}

bool MeanCluster::is_within_bounds(double value)
{
	if (abs(value) > absolute_threshold)
		return false;

	double upper_bound = mean + std::abs(bounds);
	double lower_bound = use_upper_bound_only ? 0 : mean - std::abs(bounds);

	bool is_within = value >= lower_bound && value <= upper_bound;
	return is_within;
}

void MeanCluster::reset()
{
	reset_statistics();
	data.clear();
}

void MeanCluster::reset_statistics()
{
	sum = 0.0;
	sum_squared = 0.0;

	mean = 0.0;
	variance = 0.0;
	stdev = 0.0;

	bounds = 0.0;
}
