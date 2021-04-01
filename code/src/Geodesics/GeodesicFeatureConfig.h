#pragma once

#include <vector>

struct GeodesicFeatureConfig
{
	int number_geodesic_features;
	int dynamic_geodesic_features_index;
	std::vector<float> geodesic_feature_weights;
};
