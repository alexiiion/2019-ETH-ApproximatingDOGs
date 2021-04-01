#pragma once

#include "MeshModel.h"
#include "DevelopableOptimizationSettings.h"
#include "GaussianCurvature.h"


class Developable
{
public:
	Developable() : gaussian_curvature(GaussCurvature(concatenated))
	{};

	Mesh non_developable; //target mesh with eventual upsamling
	Mesh result; //result mesh with merged seams

	Mesh concatenated; //concatenated developable parts
	std::vector<Mesh> developable_parts;
	std::vector<Mesh> developable_parts_2D;

	DevelopableOptimizationSettings optimization_settings;
	GaussCurvature gaussian_curvature;

	
	std::vector<int> assigned_labels;
	//std::vector<int> vertex_patch_assignment;			// # non_developable.V -> patch index
	std::vector<int> face_patch_assignment;				// # non_developable.F -> patch index
	std::vector<Eigen::VectorXd> patch_face_distances;  // # non_developable.F -> distance to each patch 
	std::vector<Eigen::MatrixXd> patch_closest_points;  // # non_developable.F -> closest point on each patch 


	//data for mapping between submeshes and concatenated mesh
	std::vector<std::vector<int>> global_patch_correspondance;
	std::vector<int> cumulative_offset_vertices;
	int number_patches;
	int number_concatenated_vertices;

	//all indices refer to concatenated mesh
	Eigen::MatrixXi seam_pairs;
	Eigen::VectorXi boundary_indices;
	Eigen::MatrixXd boundary_vertices_original;

	std::vector<int> inner_indices;
	std::vector<std::vector<int>> inner_indices_adjacency;
	Eigen::VectorXd inner_indices_mask; //0-1 mask marking only inner indices (for K calculation)


	int upsampling_step = 0;
	//void upsample();


	//output from objectives, should be somewhere else!
	double output_angle_defect;
	double output_stiching;


	void reset()
	{
		assigned_labels.clear();
		face_patch_assignment.clear();
		patch_face_distances.clear();
		patch_closest_points.clear();

		global_patch_correspondance.clear();
		cumulative_offset_vertices.clear();
		number_patches = 0;
		number_concatenated_vertices = 0;


		seam_pairs.resize(0, 0);
		boundary_indices.resize(0, 0);
		inner_indices.clear();
		inner_indices_adjacency.clear();
		inner_indices_mask.resize(0, 0);

		gaussian_curvature.reset();
		output_angle_defect = -1;
		output_stiching = -1;
	}

};
