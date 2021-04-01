#pragma once

#include <Eigen/Core>

//class Mesh;
#include "MeshModel.h"


class GaussCurvature
{
public:
	//GaussCurvature() {};
	GaussCurvature(Mesh& mesh) : mesh(mesh) 
	{ reset(); };

	//void update(bool exclude_border = true);
	void update(Eigen::VectorXd mask_vertices);
	void reset();

	Eigen::VectorXd get();
	double average();
	double min();
	double max();

private:
	Mesh& mesh;
	Eigen::VectorXd K;

};