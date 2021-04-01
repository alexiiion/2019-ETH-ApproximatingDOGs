#pragma once

#include <Eigen/Core>
#include <vector>

//#include "GeodesicCandidates.h"
struct GeodesicCandidates;
#include "MeshModel.h"
//class Mesh;


class MeritSelector
{
public:
	MeritSelector(std::vector<int>& pre_selected_indices);
	MeritSelector(std::vector<int>& pre_selected_indices, GeodesicCandidates& geodesics_candidates, std::vector<Mesh> associated_submeshes);

	void update_covering_mesh(Mesh covering_mesh);
	int try_get_next();

private:
	std::vector<int>& pre_selected_indices;
	GeodesicCandidates& geodesics_candidates; //debug: store longest of geodesics_directions (so 1 geodesic per vertex)
	Mesh covering_mesh;
	std::vector<Mesh> associated_submeshes;

	Eigen::VectorXd merit; //unsorted (i.e. index is geodesic path index), column entries are properties
	Eigen::MatrixXd feature_matrix;
	std::vector<double> feature_weights;
	int dynamic_feature_index;
	int number_features;

	//int current_merit = -1;
	//int current_selected_index = - 1;


	void initialize();
	void update_merit();

	void compute_dynamic_features();
	void compute_constant_features();
	void add_normalized_feature(const Eigen::VectorXd& feature, const int current_feature_index, Eigen::MatrixXd& out_features);

	Eigen::VectorXd get_lengths();
	Eigen::VectorXd get_areas();

};