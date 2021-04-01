#include "LinearTransformations.h"

#include <Eigen/Geometry>

Eigen::RowVector3d rotate_vec(const Eigen::RowVector3d& pt, const Eigen::RowVector3d& center, const Eigen::Vector3d& axis,
														 double angle) {
	Eigen::Matrix3d rot; rot =  Eigen::AngleAxis<double>(angle, axis);
    return (rot*((pt-center).transpose())).transpose()+center;
}