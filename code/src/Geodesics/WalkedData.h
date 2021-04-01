#pragma once

#include <vector>
#include <eigen/Core>

struct WalkedData
{
	std::vector<Eigen::RowVector3d> path;
	std::vector<int> faces;
	std::vector<std::pair<int, int>> vertices;
	double length = 0.0;
};
