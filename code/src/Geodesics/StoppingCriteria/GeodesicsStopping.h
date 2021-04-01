#pragma once

#include <Eigen/Core>
#include "WalkedData.h"
class Mesh;

struct GeodesicsStopping
{

	virtual bool should_stop() = 0;
	virtual void update(const WalkedData& data, const Mesh& mesh) = 0;
	virtual void reset() = 0;

	//virtual void update(const Eigen::MatrixXd& path, const Mesh& mesh, const int last_index) = 0;
};
