#include "NormalsHelper.h"

#include <igl/point_mesh_squared_distance.h>

#include "CoordinateConverter.h"
//#include "GlobalSettings.h"
//#include "Utils.h"
//#include "Logger.h"


using namespace std;

Eigen::MatrixXd utils::compute_surface_normals_vertices(const Eigen::MatrixXd& curve, const Mesh& mesh)
{
	Eigen::MatrixXd normals(curve.rows(), 3);

	Eigen::MatrixXi bary_indices;
	Eigen::MatrixXd bary_weigths;
	CoordinateConverter::barycentric_coords_from_points(curve, mesh.V, mesh.F, bary_indices, bary_weigths);

	for (int i = 0; i < curve.rows(); i++)
	{
		Eigen::Vector3d normal(0, 0, 0);

		for (int bary_i = 0; bary_i < 3; bary_i++)
		{
			double weight = bary_weigths(i, bary_i);
			int vertex_index = bary_indices(i, bary_i);

			normal += weight * mesh.normals_vertices.row(vertex_index);
		}

		normals.row(i) = normal.normalized();
	}

	return normals;
}

Eigen::MatrixXd utils::compute_surface_normals_faces(const Eigen::MatrixXd& curve, const Mesh& mesh)
{
	Eigen::VectorXi closest_faces;
	Eigen::VectorXd square_distance_unused;
	Eigen::MatrixXd closest_points_unused;
	igl::point_mesh_squared_distance(curve, mesh.V, mesh.F, square_distance_unused, closest_faces, closest_points_unused);

	Eigen::MatrixXd normals(curve.rows(), 3);
	for (int i = 0; i < curve.rows(); i++)
		normals.row(i) = mesh.normals_faces.row(closest_faces(i));

	return normals;
}
