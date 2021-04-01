#pragma once

#include <Eigen/Dense>
#include "MeshModel.h"

namespace utils
{
	Eigen::MatrixXd compute_surface_normals_vertices(const Eigen::MatrixXd& curve, const Mesh& mesh);
	Eigen::MatrixXd compute_surface_normals_faces(const Eigen::MatrixXd& curve, const Mesh& mesh);
}