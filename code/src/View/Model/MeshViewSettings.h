#pragma once

#include <Eigen/Core>
#include "MeshVisualization.h"
#include "Utils/Colors.h"

struct MeshViewSettings
{
	//MeshViewSettings() {};
	//Mesh& mesh;

	int view_index = -1;
	bool has_changed = true;

	bool show_mesh = true;
	bool show_faces = true;
	bool show_boundary = false;
	bool show_wireframe = true;

	bool label_faces = false;
	bool label_vertices = false;

	Eigen::MatrixXd cached_V;

	MeshVisualization visualization = MeshVisualization::Standard;
	Eigen::MatrixXd fill;

	Eigen::Vector3d standard_fill_color = Colors::GRAY_ULTRALIGHT;
	
	//TODO delete
	Eigen::RowVectorXd dimensions;
};

struct PatchViewSettings : MeshViewSettings
{
	bool label_quad_faces = false;
	bool show_position_constraints = true;
	bool show_target_constraints = true;
};
