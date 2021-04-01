#pragma once

#include <Eigen/Core>
#include <vector>

class Mesh;
//#include "MeshModel.h"


class Coverage
{
public:
	//Coverage(const Mesh& mesh, const std::vector<Mesh>& covering_meshes, float* coverage_theshold);
	Coverage(const Mesh& mesh, const std::vector<Mesh>& covering_meshes, float* coverage_theshold);
	Coverage(const Mesh& mesh, float* coverage_theshold);

	//metrics
	double hausdorff_distance() { return _hausdorff_distance; };
	double root_mean_square_error() { return sqrt(distances_squared().mean()); };
	double distance_sum() { return _distance_sum; };
	double distance_percent() { return _distance_percent; };

	// #V *mesh*, distance from each *mesh* vertex to the closest point on any surface *covering_meshes*
	Eigen::VectorXd distances();
	Eigen::VectorXd distances_squared();

	Eigen::VectorXd covered_vertices_mask() { return _covered_vertices_mask; };
	
	// #V *mesh*, for each *from_mesh* vertex the closest point on any surface in *covering_meshes*
	Eigen::MatrixXd covering_points() { return _covering_points; }; 
	Eigen::VectorXi covering_faces() { return _covering_faces; };

	std::vector<int> covered_vertices() { return _covered_vertices; };
	std::vector<int> uncovered_vertices() { return _uncovered_vertices; };


	void compute(const std::vector<Mesh>& covering_meshes);
	void compute(const Mesh& covering_mesh);
	bool is_covered() { return is_computed && _uncovered_vertices.size() < 1; };

	//indices refer to *mesh*
	//void get_uncovered_areas(std::vector<std::vector<int>>& uncovered_components, int& projected_centroid);
	//std::vector<Mesh> get_uncovered_areas();
	//
	//void get_uncovered_areas(std::vector<std::vector<int>>& uncovered_components, int& projected_centroid, std::vector<Mesh>& uncovered_areas);


private:
	const Mesh& mesh;
	float* coverage_threshold;

	bool is_computed = false;
	double _hausdorff_distance = 1e6;
	double _distance_sum = 1e6;
	double _distance_mean = 1e6;
	double _distance_percent = 1e6;

	std::vector<int> _covered_vertices;
	std::vector<int> _uncovered_vertices;

	Eigen::VectorXd _distances_squared;
	Eigen::MatrixXd _covering_points; // #V *from_mesh*, from each target vertex to the closest point on any surface (*to_meshes*)
	Eigen::VectorXi _covering_faces;
	Eigen::VectorXd _covered_vertices_mask;

	Eigen::VectorXd vertex_mass;

	void reset();
	void compute_mass();

	void compute_hausdorff_distance(const Mesh& covering_mesh);
	void compute_distance_sum();
	void compute_distance_percent();

};