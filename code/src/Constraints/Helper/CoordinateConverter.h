#pragma once

#include <igl/barycentric_coordinates.h>
#include <igl/point_mesh_squared_distance.h>

#include <vector>
#include "Logger.h"

using namespace std;

class CoordinateConverter {

public:

	static void barycentric_coords_from_points(const Eigen::MatrixXd& points, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& out_barycentric_indices, Eigen::MatrixXd& out_barycentric_weights)
	{
		Eigen::VectorXi closest_faces;
		Eigen::VectorXd square_distance_unused;
		Eigen::MatrixXd closest_points_unused;
		igl::point_mesh_squared_distance(points, V, F, square_distance_unused, closest_faces, closest_points_unused);

		barycentric_coords_from_points(points, closest_faces, V, F, out_barycentric_indices, out_barycentric_weights);
	}

	static void barycentric_coords_from_points(const Eigen::MatrixXd& points, const Eigen::VectorXi& closest_faces, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& out_barycentric_indices, Eigen::MatrixXd& out_barycentric_weights)
	{
		if (!points.rows())
		{
			write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> :: " << endl << "   tried to create barycentric coordinates from non-existent points!" << endl << endl;
			return;
		}


		const int log = 6;
		write_log(log) << "closest faces: " << endl << closest_faces << endl << endl;

		const int size = points.rows();
		out_barycentric_indices.resize(size, 3);
		out_barycentric_weights.resize(size, 3);

		const int number_vertices = 3;
		const int number_faces = closest_faces.rows();

		vector<Eigen::MatrixXd> barycentric_vertices;

		for (int i = 0; i < number_vertices; i++) {
			Eigen::MatrixXd vertices(number_faces, 3);
			barycentric_vertices.push_back(vertices);
		}

		for (int i = 0; i < number_faces; i++) {
			auto face_index = closest_faces(i);
			auto face = F.row(face_index);
			write_log(log) << i << ": face[" << face_index << "]: " << face << endl;

			for (int j = 0; j < number_vertices; j++) {
				out_barycentric_indices(i, j) = face(j);
				barycentric_vertices[j].row(i) = V.row(face(j));
			}
		}

		igl::barycentric_coordinates(points, barycentric_vertices[0], barycentric_vertices[1], barycentric_vertices[2], out_barycentric_weights);

		write_log(log) << endl << "barycentric -- vertex weights: " << endl << out_barycentric_weights << endl << endl;
		write_log(log) << "barycentric -- vertex indices: " << endl << out_barycentric_indices << endl << endl;
	};

	static void points_from_barycentric_coords(const Eigen::MatrixXi& barycentric_indices, const Eigen::MatrixXd& barycentric_weights, const Eigen::MatrixXd& V, Eigen::MatrixXd& out_points)
	{
		if (!barycentric_indices.rows())
		{
			write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> :: " << endl << "   tried to create points from non-existent barycentric coordinates!" << endl << endl;
			return;
		}

		const int log = 6;

		const int number_vertices = 3;
		const int size = barycentric_indices.rows();
		out_points.resize(size, 3);

		for (int i = 0; i < size; i++) {
			Eigen::RowVector3d reconstructed_point;
			reconstructed_point.setZero();

			for (int j = 0; j < number_vertices; j++) {
				auto weight = barycentric_weights(i, j);
				auto vertex_index = barycentric_indices(i, j);
				auto vertex = V.row(vertex_index);
				write_log(log) << "   " << i << ": weight=" << weight << ", vi=" << vertex_index << endl;

				reconstructed_point += weight * vertex;
			}

			out_points.row(i) = reconstructed_point;
			write_log(log) << i << ": reconstructed point: " << reconstructed_point << endl;
		}
	};
};
