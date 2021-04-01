#pragma once

#include <vector>
#include "Eigen/Core"//#include "MeshModel.h";

namespace label
{

	//bool try_get_flat_area(const int& source_vertex, const Mesh& mesh, std::vector<int>& flat_vertices, int& center_vertex);

	std::vector<int> dfs_flat_area(const int& start_index, const std::vector<std::vector<int>>& connectivity_list, const Eigen::MatrixXd& N);

}