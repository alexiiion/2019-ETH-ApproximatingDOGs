#include "Coverage.h"

#include <igl/point_mesh_squared_distance.h>
#include <igl/hausdorff.h>
#include <igl/massmatrix.h>

#include "MeshModel.h"
#include "MeshController.h"


Coverage::Coverage(const Mesh& mesh, const std::vector<Mesh>& covering_meshes, float* coverage_theshold)
	: mesh(mesh)
{
	this->coverage_threshold = coverage_theshold;
	compute_mass();
	compute(covering_meshes);
}
Coverage::Coverage(const Mesh& mesh, float* coverage_theshold)
	: mesh(mesh)
{
	this->coverage_threshold = coverage_theshold;
	compute_mass();
}

void Coverage::reset()
{
	if (!is_computed)
		return;

	_covered_vertices.clear();
	_uncovered_vertices.clear();

	_distances_squared.resize(0);
	_covering_points.resize(0, Eigen::NoChange);
	_covering_faces.resize(0);
	_covered_vertices_mask.resize(0);

	is_computed = false;
}

Eigen::VectorXd Coverage::distances()
{
	if (_distances_squared.rows() > 0)
		return _distances_squared.cwiseSqrt();

	Eigen::VectorXd empty(this->mesh.V.rows());
	empty.setConstant(100000);
	return empty;
}

Eigen::VectorXd Coverage::distances_squared()
{
	if (_distances_squared.rows() > 0)
		return _distances_squared;

	Eigen::VectorXd empty(this->mesh.V.rows());
	empty.setConstant(100000);
	return empty;
}

void Coverage::compute(const std::vector<Mesh>& covering_meshes)
{
	if (covering_meshes.size() < 1)
		return;

	Mesh covering_mesh;
	if (covering_meshes.size() == 1)
		covering_mesh = covering_meshes[0];
	else
		covering_mesh = meshhelper::concatenate_meshes(covering_meshes);

	compute(covering_mesh);
}

void Coverage::compute(const Mesh& covering_mesh)
{
	if (covering_mesh.V.rows() < 1)
		return;

	reset();

	//find closest point on _covering_ faces (as defined by update_covering_wrappers()) from all target vertices
	igl::point_mesh_squared_distance(mesh.V, covering_mesh.V, covering_mesh.F, _distances_squared, _covering_faces, _covering_points);

	//set coverage
	const float max_distance = (*coverage_threshold) * (*coverage_threshold);
	const int number_vertices = mesh.V.rows();

	_covered_vertices_mask.resize(number_vertices);
	_covered_vertices_mask.setOnes();

	for (int i = 0; i < number_vertices; i++)
	{
		if (_distances_squared(i) > max_distance)
		{
			_uncovered_vertices.push_back(i);
			_covered_vertices_mask(i) = 0;
		}
		else
		{
			_covered_vertices.push_back(i);
		}
	}

	compute_hausdorff_distance(covering_mesh);
	compute_distance_sum();
	compute_distance_percent();

	is_computed = true;
}

void Coverage::compute_hausdorff_distance(const Mesh& covering_mesh)
{
	igl::hausdorff(mesh.V, mesh.F, covering_mesh.V, covering_mesh.F, _hausdorff_distance);
}

void Coverage::compute_distance_sum()
{
	Eigen::VectorXd weighted_distances = (_distances_squared.cwiseSqrt() * vertex_mass);
	_distance_sum = weighted_distances.sum();
	_distance_mean = weighted_distances.mean();
}

void Coverage::compute_distance_percent()
{
	double area = vertex_mass.sum();
	if (abs(area) < 1e-6)
		return;

	_distance_percent = _distance_sum / area;
}

void Coverage::compute_mass()
{
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(mesh.V, mesh.F, igl::MASSMATRIX_TYPE_VORONOI, M);
	
	vertex_mass = M.diagonal();
}
