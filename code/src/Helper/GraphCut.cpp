#include "GraphCut.h"

#include <iostream>
#include <igl/readOFF.h>
#include <igl/edges.h>
#include <tuple>
#include <fstream>
#include <iterator>
#include <ctime>
//#include <MeshModel.h>

#undef GCO_ENERGYTYPE32
#include "GCoptimization.h"

#include "Utils.h"


using std::get;
using std::cout;
using std::endl;

using namespace std; 


Eigen::MatrixXd D;
Eigen::MatrixXd N;

double scale = 10000000;
//double scale = 10000;
double cover_smoothness = 50;
double max_value;
double max_distance = 30;

int smoothnessAdaptiveTerm(const int n1, const int n2, const int l1, const int l2)
{
	//currently UNUSED

	if (l1 == l2)
		return 0;

	double distance1 = D(l1, n1);
	double distance2 = D(l2, n2);
	int data = cover_smoothness + pow(max(distance1, distance2), 4);

	if (data < 0 || data > GCO_MAX_ENERGYTERM)
		return GCO_MAX_ENERGYTERM;

	return data;
}

int smoothnessTerm(const int n1, const int n2, const int l1, const int l2)
{
	if (l1 == l2)
		return 0;

	return cover_smoothness;
}

int dataTerm(const int node, const int label)
{
	double distance = D(label, node);
	//if(distance > max_distance)
	//	return GCO_MAX_ENERGYTERM;

	//double angle = N(label, node);
	//cout << "(" << label << "," << node <<") graph cut:: data term:: angle = " << angle << ", distance = " << distance << endl;

	int data = distance * scale;
	//int data = (distance*distance)/max_value * scale;

	if (data < 0 || data > GCO_MAX_ENERGYTERM)
		return GCO_MAX_ENERGYTERM;

	return data;
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness)
{
	std::vector<double> resulting_distances;
	return compute_graph_cut_labeling(mesh, patch_vertex_distances, smoothness, resulting_distances);
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness, std::vector<double>& out_resulting_distances)
{
	std::vector<Eigen::MatrixXd> empty_normals;
	return compute_graph_cut_labeling(mesh, patch_vertex_distances, empty_normals, smoothness, out_resulting_distances);
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, const std::vector<Eigen::MatrixXd>& patch_vertex_normals, double smoothness, std::vector<double>& out_resulting_distances)
{
	const int number_assignments = mesh.V.rows();
	const int number_labels = patch_vertex_distances.size();

	D.resize(number_labels, number_assignments);

	for (int i = 0; i < patch_vertex_distances.size(); ++i)
	{
		Eigen::VectorXd distances = patch_vertex_distances[i];
		for (int j = 0; j < distances.rows(); j++)
			D(i, j) = distances(j);
	}

	max_value = D.maxCoeff();
	D /= max_value;
	std::cout << "GraphCut:: first max_value: " << max_value << std::endl;

	for (int i = 0; i < D.rows(); ++i)
	{
		for (int j = 0; j < D.cols(); j++)
		{
			D(i, j) = std::min(std::exp(D(i, j) * 10), 100000.0);
		}
	}

	cover_smoothness = smoothness;
	scale = 10000000;
	//max_distance = std::exp(1.50 * 10);

	max_value = D.maxCoeff();
	scale /= max_value;
	std::cout << "GraphCut:: final max_value: " << max_value << ", scale: " << scale << std::endl;



	N.resize(number_labels, number_assignments);

	for (int i = 0; i < patch_vertex_normals.size(); ++i)
	{
		Eigen::MatrixXd normals = patch_vertex_normals[i];

		for (int j = 0; j < normals.rows(); j++)
		{
			Eigen::VectorXd normal_mesh = mesh.normals_vertices.row(j);
			Eigen::VectorXd normal_patch = normals.row(j);

			N(i, j) = angle(normal_mesh, normal_patch);
			//cout << "(" << i << "," << j << ")  normal_mesh: " << normal_mesh.transpose() << "\nnormal_patch: " << normal_patch.transpose() << " --> angle = " << N(i, j) << endl;
		}
	}

	//cout << "D: \n" << D << "\n\nN: \n" << N << endl;


	GCoptimizationGeneralGraph mrf(number_assignments, number_labels);
	//GCoptimizationGeneralGraph mrf((int)mesh.V.rows(), (int)patch_vertex_distances.size());

	Eigen::MatrixXi E;
	igl::edges(mesh.F, E);

	//Eigen::MatrixXi connectivity = E;

	for (int i = 0; i < E.rows(); ++i)
	{
		mrf.setNeighbors(E(i, 0), E(i, 1));
	}

	std::srand(unsigned(std::time(0)));

	mrf.setDataCost(dataTerm);
	mrf.setSmoothCost(smoothnessTerm);
	mrf.setLabelCost(0);
	mrf.setVerbosity(1);
	mrf.setLabelOrder(true);

	cout << "    start optimization" << endl;

	auto energy = mrf.expansion();

	auto en_data = mrf.giveDataEnergy();
	auto en_smooth = mrf.giveSmoothEnergy();
	auto en_label = mrf.giveLabelEnergy();

	cout << "    Energy: " << energy << ", Data Energy: " << en_data << ", Smoothness Energy: " << en_smooth << ", Label Energy: " << en_label << endl;

	std::vector<int> new_labels(mesh.V.rows());

	mrf.whatLabel(0, mesh.V.rows(), new_labels.data());

	//resulting_distances.resize(mesh.V.rows());
	out_resulting_distances.clear();
	for (int i = 0; i < mesh.V.rows(); i++)
		out_resulting_distances.push_back(D(new_labels[i], i));

	return new_labels;
}

//UNUSED
std::vector<int> graphcut::compute_graph_cut_face_labeling(const Eigen::MatrixXd& face_midpoints, const Eigen::MatrixXi& connectivity, const std::vector<Eigen::VectorXd>& patch_face_distances, double smoothness, std::vector<double>& out_resulting_distances)
{
	const int number_assignments = face_midpoints.rows();
	const int number_labels = patch_face_distances.size();

	if (number_labels == 1)
	{
		out_resulting_distances.clear();
		out_resulting_distances.resize(number_assignments, 0);

		std::vector<int> new_labels(number_assignments, 0);
		return new_labels;
	}

	D.resize(number_labels, number_assignments);

	for (int i = 0; i < patch_face_distances.size(); ++i)
	{
		Eigen::VectorXd distances = patch_face_distances[i];
		for (int j = 0; j < distances.rows(); j++)
			D(i, j) = distances(j);
	}

	max_value = D.maxCoeff();
	D /= max_value;
	std::cout << "GraphCut:: first max_value: " << max_value << std::endl;

	//for (int i = 0; i < D.rows(); ++i)
	//{
	//	for (int j = 0; j < D.cols(); j++)
	//		D(i, j) = std::min(std::exp(D(i, j) * 10), 100000.0);
	//}

	cover_smoothness = smoothness;
	scale = 10000000;

	max_value = D.maxCoeff();
	scale /= max_value;
	std::cout << "GraphCut:: final max_value: " << max_value << ", scale: " << scale << std::endl;

	GCoptimizationGeneralGraph mrf(number_assignments, number_labels);

	//Eigen::MatrixXi E;
	//igl::edges(F, E);

	//for (int i = 0; i < E.rows(); ++i)
	//	mrf.setNeighbors(E(i, 0), E(i, 1));


	for (int i = 0; i < connectivity.rows(); ++i)
		mrf.setNeighbors(connectivity(i, 0), connectivity(i, 1));


	//int dim = adjacency_FF.cols();
	//for (int r = 0; r < adjacency_FF.rows(); r++)
	//{
	//	for (int c = 0; c < dim; c++)
	//		if(adjacency_FF(r, c) >= 0 && adjacency_FF(r, (c + 1) % dim) >= 0)
	//			mrf.setNeighbors(adjacency_FF(r, c), adjacency_FF(r, (c+1)%dim));
	//}

	std::srand(unsigned(std::time(0)));

	mrf.setDataCost(dataTerm);
	mrf.setSmoothCost(smoothnessAdaptiveTerm);
	mrf.setLabelCost(0);
	mrf.setVerbosity(1);
	mrf.setLabelOrder(true);

	cout << "    start optimization" << endl;

	auto energy = mrf.expansion();

	auto en_data = mrf.giveDataEnergy();
	auto en_smooth = mrf.giveSmoothEnergy();
	auto en_label = mrf.giveLabelEnergy();

	cout << "    Energy: " << energy << ", Data Energy: " << en_data << ", Smoothness Energy: " << en_smooth << ", Label Energy: " << en_label << endl;

	std::vector<int> new_labels(number_assignments);
	mrf.whatLabel(0, number_assignments, new_labels.data());

	//resulting_distances.resize(number_assignments);
	out_resulting_distances.clear();
	for (int i = 0; i < number_assignments; i++)
		out_resulting_distances.push_back(D(new_labels[i], i));

	return new_labels;
}



/*
Eigen::MatrixXd D;
Eigen::MatrixXd N;
Eigen::MatrixXd DN;

double scale = 10000000;
//double scale = 10000;
double cover_smoothness = 50;
double max_value;
double max_distance = 30;

int smoothnessAdaptiveTerm(const int n1, const int n2, const int l1, const int l2)
{
	//currently UNUSED

	if (l1 == l2) 
		return 0;
	
	double distance1 = D(l1, n1);
	double distance2 = D(l2, n2);
	int data = cover_smoothness + pow(max(distance1, distance2), 4);

	if (data < 0 || data > GCO_MAX_ENERGYTERM)
		return GCO_MAX_ENERGYTERM;

	return data;
}

int smoothnessTerm(const int n1, const int n2, const int l1, const int l2)
{
	if (l1 == l2) 
		return 0;
	
	return cover_smoothness;
}

int dataTerm(const int node, const int label)
{
	double distance = DN(label, node);
	int data = distance * scale;

	if (data < 0 || data > GCO_MAX_ENERGYTERM)
		return GCO_MAX_ENERGYTERM;

	return data;
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness)
{
	std::vector<double> resulting_distances;
	return compute_graph_cut_labeling(mesh, patch_vertex_distances, smoothness, resulting_distances);
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, double smoothness, std::vector<double>& out_resulting_distances)
{
	std::vector<Eigen::VectorXd> empty_normals;
	return compute_graph_cut_labeling(mesh, patch_vertex_distances, empty_normals, smoothness, out_resulting_distances);
}

std::vector<int> graphcut::compute_graph_cut_labeling(const Mesh& mesh, const std::vector<Eigen::VectorXd>& patch_vertex_distances, const std::vector<Eigen::VectorXd>& patch_vertex_normals, double smoothness, std::vector<double>& out_resulting_distances)
{
	const double multiplier = 10;
	double normalizer = 1;

	const int number_labels = patch_vertex_distances.size();
	const int number_assignments = patch_vertex_distances[0].rows();


	data_to_matrix(patch_vertex_distances, D);
	double distanceThreshold = 1.5 / D.maxCoeff() * normalizer;
	distanceThreshold = std::exp(distanceThreshold * multiplier);
	std::cout << "GraphCut:: distanceThreshold: " << distanceThreshold << std::endl;

	normalize_data(D, normalizer);
	transform_data_exponentially(D, multiplier, distanceThreshold);


	if (patch_vertex_normals.size() > 0)
	{
		data_to_matrix(patch_vertex_normals, N);
		double deviationThreshold = 0.3 * M_PI / N.maxCoeff() * normalizer;
		deviationThreshold = std::exp(deviationThreshold * multiplier);
		std::cout << "GraphCut:: deviationThreshold: " << deviationThreshold << std::endl;
	
		normalize_data(N, normalizer);
		transform_data_exponentially(N, multiplier, deviationThreshold);

		normalizer = D.maxCoeff();
		normalize_data(N, normalizer);

		//combine
		DN.resize(number_labels, number_assignments);
		for (int i = 0; i < DN.rows(); ++i)
			for (int j = 0; j < DN.cols(); j++)
				DN(i, j) = std::max(D(i, j), N(i, j));
	}
	else
	{
		DN = D;
	}

	max_value = DN.maxCoeff();
	scale = 1000000;
	scale /= max_value;
	std::cout << "GraphCut:: final max_value: " << max_value << ", scale: " << scale << std::endl;

	cover_smoothness = smoothness;


	GCoptimizationGeneralGraph mrf(number_assignments, number_labels);
	//GCoptimizationGeneralGraph mrf((int)mesh.V.rows(), (int)patch_vertex_distances.size());

	Eigen::MatrixXi E;
	igl::edges(mesh.F, E);

	//Eigen::MatrixXi connectivity = E;

	for (int i = 0; i < E.rows(); ++i)
	{
		mrf.setNeighbors(E(i, 0), E(i, 1));
	}

	std::srand(unsigned(std::time(0)));

	mrf.setDataCost(dataTerm);
	mrf.setSmoothCost(smoothnessTerm);
	mrf.setLabelCost(0);
	mrf.setVerbosity(1);
	mrf.setLabelOrder(true);

	cout << "    start optimization" << endl;

	auto energy = mrf.expansion();

	auto en_data = mrf.giveDataEnergy();
	auto en_smooth = mrf.giveSmoothEnergy();
	auto en_label = mrf.giveLabelEnergy();

	cout << "    Energy: " << energy << ", Data Energy: " << en_data << ", Smoothness Energy: " << en_smooth << ", Label Energy: " << en_label << endl;

	std::vector<int> new_labels(mesh.V.rows());

	mrf.whatLabel(0, mesh.V.rows(), new_labels.data());

	//resulting_distances.resize(mesh.V.rows());
	out_resulting_distances.clear();
	for (int i = 0; i < mesh.V.rows(); i++)
		out_resulting_distances.push_back(D(new_labels[i], i));
		
	return new_labels;
}

std::vector<int> graphcut::compute_graph_cut_face_labeling(const Eigen::MatrixXd& face_midpoints, const Eigen::MatrixXi& connectivity, const std::vector<Eigen::VectorXd>& patch_face_distances, double smoothness, std::vector<double>& out_resulting_distances)
{
	const int number_assignments = face_midpoints.rows();
	const int number_labels = patch_face_distances.size();

	if (number_labels == 1)
	{
		out_resulting_distances.clear();
		out_resulting_distances.resize(number_assignments, 0);

		std::vector<int> new_labels(number_assignments, 0);
		return new_labels;
	}

	D.resize(number_labels, number_assignments);

	for (int i = 0; i < patch_face_distances.size(); ++i)
	{
		Eigen::VectorXd distances = patch_face_distances[i];
		for (int j = 0; j < distances.rows(); j++)
			D(i, j) = distances(j);
	}

	max_value = D.maxCoeff();
	D /= max_value;
	std::cout << "GraphCut:: first max_value: " << max_value << std::endl;

	//for (int i = 0; i < D.rows(); ++i)
	//{
	//	for (int j = 0; j < D.cols(); j++)
	//		D(i, j) = std::min(std::exp(D(i, j) * 10), 100000.0);
	//}

	cover_smoothness = smoothness;
	scale = 10000000;

	max_value = D.maxCoeff();
	scale /= max_value;
	std::cout << "GraphCut:: final max_value: " << max_value << ", scale: " << scale << std::endl;

	GCoptimizationGeneralGraph mrf(number_assignments, number_labels);

	//Eigen::MatrixXi E;
	//igl::edges(F, E);

	//for (int i = 0; i < E.rows(); ++i)
	//	mrf.setNeighbors(E(i, 0), E(i, 1));


	for (int i = 0; i < connectivity.rows(); ++i)
		mrf.setNeighbors(connectivity(i, 0), connectivity(i, 1));


	//int dim = adjacency_FF.cols();
	//for (int r = 0; r < adjacency_FF.rows(); r++)
	//{
	//	for (int c = 0; c < dim; c++)
	//		if(adjacency_FF(r, c) >= 0 && adjacency_FF(r, (c + 1) % dim) >= 0)
	//			mrf.setNeighbors(adjacency_FF(r, c), adjacency_FF(r, (c+1)%dim));
	//}

	std::srand(unsigned(std::time(0)));

	mrf.setDataCost(dataTerm);
	mrf.setSmoothCost(smoothnessAdaptiveTerm);
	mrf.setLabelCost(0);
	mrf.setVerbosity(1);
	mrf.setLabelOrder(true);

	cout << "    start optimization" << endl;

	auto energy = mrf.expansion();

	auto en_data = mrf.giveDataEnergy();
	auto en_smooth = mrf.giveSmoothEnergy();
	auto en_label = mrf.giveLabelEnergy();

	cout << "    Energy: " << energy << ", Data Energy: " << en_data << ", Smoothness Energy: " << en_smooth << ", Label Energy: " << en_label << endl;

	std::vector<int> new_labels(number_assignments);
	mrf.whatLabel(0, number_assignments, new_labels.data());

	//resulting_distances.resize(number_assignments);
	out_resulting_distances.clear();
	for (int i = 0; i < number_assignments; i++)
		out_resulting_distances.push_back(D(new_labels[i], i));

	return new_labels;

}

void graphcut::data_to_matrix(const std::vector<Eigen::VectorXd>& data, Eigen::MatrixXd& matrix)
{
	if (data.size() < 1)
		return;

	const int number_labels = data.size();
	const int number_assignments = data[0].rows();

	matrix.resize(number_labels, number_assignments);

	for (int i = 0; i < data.size(); ++i)
	{
		Eigen::VectorXd distances = data[i];
		matrix.row(i) = distances.transpose();
		//for (int j = 0; j < distances.rows(); j++)
		//	matrix(i, j) = distances(j);
	}
}

void graphcut::normalize_data(Eigen::MatrixXd& data, const double& max)
{
	double max_value = data.maxCoeff();
	data = data / max_value * max;
	std::cout << "GraphCut:: normalize_data to " << max << ", old max = " << max_value  << ", new max = " << data.maxCoeff() << std::endl;
}

void graphcut::transform_data_exponentially(Eigen::MatrixXd& data, const double& multiplier, const double& cutoff)
{
	for (int i = 0; i < data.rows(); ++i)
	{
		for (int j = 0; j < data.cols(); j++)
		{
			data(i, j) = std::min(std::exp(data(i, j) * multiplier), cutoff);
		}
	}
}
*/