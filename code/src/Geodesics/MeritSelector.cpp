#include "MeritSelector.h"

#include <igl/point_mesh_squared_distance.h>

#include "GeodesicCandidates.h"
#include "Coverage.h"

#include "MeshController.h"
#include "CurveHelper.h"
#include "Utils.h"
#include "Logger.h"

using namespace std;


MeritSelector::MeritSelector(std::vector<int>& pre_selected_indices, GeodesicCandidates& geodesics_candidates, std::vector<Mesh> associated_submeshes)
	: pre_selected_indices(pre_selected_indices), geodesics_candidates(geodesics_candidates), associated_submeshes(associated_submeshes)
{
	feature_weights =
	{
		0.3,	//length
		0.1,	//area
		0.6		//coverage
	};
	number_features = feature_weights.size();
	dynamic_feature_index = 2;

	initialize();
}

MeritSelector::MeritSelector(std::vector<int>& pre_selected_indices)
	: pre_selected_indices(pre_selected_indices), geodesics_candidates(GeodesicCandidates())
{
	number_features = 0;
	initialize();
}

void MeritSelector::update_covering_mesh(Mesh covering_mesh)
{
	this->covering_mesh = covering_mesh;
}

//returns next index referring to *pre_selected_indices*, -1 if there is none to select anymore
int MeritSelector::try_get_next()
{
	update_merit();

	int best_index;
	double max_merit = merit.maxCoeff(&best_index);
	write_log(4) << endl << "MERIT selected index: " << best_index << " (merit = " << max_merit << ") from: " << list_to_string(pre_selected_indices) << endl;

	if (max_merit == -1)
		return -1;

	if(feature_matrix.hasNaN())
		write_log(1) << "ERROR: feature matrix has NaN!" << endl;
	if(merit.hasNaN())
		write_log(1) << "ERROR: MERIT has NaN!" << endl;

	merit(best_index) = -1; //invalidate used geodesic
	return pre_selected_indices[best_index];
}

void MeritSelector::initialize()
{
	feature_matrix.resize(pre_selected_indices.size(), number_features);
	merit.resize(pre_selected_indices.size());

	//compute constant features: length, area
	compute_constant_features();
	update_merit();
}

void MeritSelector::update_merit()
{
	write_log(3) << "updating merit..." << endl;

	//dynamic: coverage of associated area? --> coverage.sum() / area
	compute_dynamic_features();

	double max_merit = 0.0;

	//update merit of geodesics
	for (int i = 0; i < pre_selected_indices.size(); i++)
	{
		if (merit(i) == -1)
			continue;

		double current_merit = 0;

		for (int j = 0; j < number_features; j++)
		{
			double weight = feature_weights[j];
			double value = feature_matrix(i, j);

			if (weight < 0)
			{
				value = 1 - value;
				weight *= -1;
			}

			current_merit += weight * value;
		}

		merit(i) = current_merit;

		if (current_merit > max_merit)
			max_merit = current_merit;
	}

	write_log(4) << "features: " << endl << feature_matrix << endl;
	write_log(4) << "merit: " << endl << merit << endl;
	write_log(4) << "max_merit: " << max_merit << endl;

	write_log(3) << "  ...done updating merit." << endl;
}

void MeritSelector::compute_dynamic_features()
{
	if (number_features < 1)
		return;

	Eigen::VectorXd coverage(pre_selected_indices.size());
	coverage.setZero();
	feature_matrix.col(dynamic_feature_index) = coverage;

	if (covering_mesh.V.rows() < 1)
		return;

	for (int i = 0; i < pre_selected_indices.size(); i++)
	{
		Mesh& mesh = associated_submeshes[i];

		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd covering_points; // #V *from_mesh*, from each target vertex to the closest point on any surface (*to_meshes*)
		Eigen::VectorXi covering_faces;
		igl::point_mesh_squared_distance(mesh.V, covering_mesh.V, covering_mesh.F, distances_squared, covering_faces, covering_points);

		double area = meshhelper::compute_face_area(mesh.V, mesh.F);
		if (abs(area) < 1e-6)
			area = 1e-6;

		coverage(i) = sqrt(distances_squared.sum()) / area;

		if(std::isnan(coverage(i)))
			write_log(0) << "check merit coverage computation: coverage = " << coverage(i) << " = " << sqrt(distances_squared.sum()) << " / " << area << "  (sqrt(distances_squared.sum()) / area)" << endl;
	}

	add_normalized_feature(coverage, dynamic_feature_index, feature_matrix);

	//write_log(0) << "feature coverage: " << linebreak << coverage << linebreak << endl;
	//write_log(0) << "FULL feature matrix: " << linebreak << feature_matrix << linebreak << endl;
}

void MeritSelector::compute_constant_features()
{
	if (number_features < 1)
		return;

	int current_feature_index = 0;

	add_normalized_feature(get_lengths(), current_feature_index++, feature_matrix);
	add_normalized_feature(get_areas(), current_feature_index++, feature_matrix);

	//write_log(0) << "geodesic_constant_features " << endl << feature_matrix << endl;
}

void MeritSelector::add_normalized_feature(const Eigen::VectorXd& feature, const int current_feature_index, Eigen::MatrixXd& out_features)
{
	double max = feature.maxCoeff();

	if (abs(max) < 1e-6) //avoid division by zero
		max = 1e-6;

	for (int i = 0; i < feature.rows(); i++)
		out_features(i, current_feature_index) = feature(i) / max;
}

Eigen::VectorXd MeritSelector::get_lengths()
{
	Eigen::VectorXd lengths(pre_selected_indices.size());
	lengths.setZero();

	for (int i = 0; i < pre_selected_indices.size(); i++)
	{
		int vi = pre_selected_indices[i];
		double length = geodesics_candidates.lengths[vi];
		lengths(i) = length;
	}

	return lengths;
}

Eigen::VectorXd MeritSelector::get_areas()
{
	Eigen::VectorXd areas(pre_selected_indices.size());
	areas.setZero();

	for (int i = 0; i < pre_selected_indices.size(); i++)
	{
		Mesh& mesh = associated_submeshes[i];
		double area = meshhelper::compute_face_area(mesh.V, mesh.F);
		areas(i) = area;
	}

	return areas;
}