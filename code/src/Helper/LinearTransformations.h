#pragma once
#include <Eigen/Dense>

Eigen::RowVector3d rotate_vec(const Eigen::RowVector3d& pt, const Eigen::RowVector3d& center, const Eigen::Vector3d& axis, double angle);