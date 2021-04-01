#pragma once

#include "MeshModel.h"
#include "GeodesicWalker.h"

namespace temp
{
	void get_geodesic(const int vertex, Mesh& target, GeodesicWalker& walker, const std::vector<int>& search_directions, Eigen::MatrixXd& out_geodesic, double& out_length);
}
