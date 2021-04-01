#pragma once

#include <Eigen/Core>
#include <vector>

struct GeodesicCandidates
{
public:
	GeodesicCandidates() {};
	GeodesicCandidates(int capacity) 
	{
		paths.reserve(capacity);

		lengths.resize(capacity);
		//gauss_curvatures.resize(capacity);
		//curvatures.resize(capacity);

		number_geodesics = 0;
		//current_geodesic_index = -1;
		//current_merit = -1;

		//geodesics_merit.resize(capacity);
	};

	void clear()
	{
		paths.clear();
		flat_geodesics.clear();

		lengths.resize(0);
		//gauss_curvatures.resize(0);
		//curvatures.resize(0);

		number_geodesics = 0;
		//current_geodesic_index = -1;
		//current_merit = -1;

		//geodesics_merit.resize(0);
	}

	bool is_empty()
	{
		return paths.size() < 1;
	}

	int number_geodesics;

	//all geodesics to select from
	std::vector<Eigen::MatrixXd> paths;
	//Eigen::VectorXd coverage;

	std::vector<int> flat_geodesics;

	//geodesic features
	Eigen::VectorXd lengths;
	Eigen::VectorXi lengths_lookup; //in descending order

	//Eigen::VectorXd gauss_curvatures;
	//Eigen::VectorXi gauss_curvatures_lookup; //in descending order

	//Eigen::VectorXd curvatures;
	//Eigen::VectorXi curvatures_lookup; //in descending order


	////feature vectors for all geodesics
	//int current_merit;
	//int current_geodesic_index;
	//Eigen::VectorXd geodesics_merit; //unsorted (i.e. index is geodesic path index), column entries are properties
	//Eigen::MatrixXd geodesic_features;
	////std::vector<int> used_geodesics;


	//std::vector<std::vector<double>> pairwise_similarity;
	//double min_similarity = 1e6;
	//int min_similarity_vertex = -1;
};
