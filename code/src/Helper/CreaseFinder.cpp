#include "CreaseFinder.h"

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

#include "Logger.h"

using namespace std;

void get_ring_neighborhood(int number_rings);
void get_geodesic_neighborhood(double geodesic_distance);

void crease::creases_ring_neighbors(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<std::vector<int>>& adjacency_VF, int number_rings, Eigen::VectorXd& out_creases)
{
	const int loglevel = 6;

	out_creases.resize(V.rows());
	
	Eigen::MatrixXd N;
	igl::per_face_normals(V, F, N);

	//out_normal = N.row(closest_faces(0)) + point;

	for (int i = 0; i < V.rows(); i++)
	{
		vector<int> ring_neighbors = adjacency_VF[i];
		if (ring_neighbors.size() < 2)
			continue;

		Eigen::VectorXd angles(ring_neighbors.size() - 1);

		write_log(loglevel) << "[" << i << "]: "; log_list(loglevel, ring_neighbors, "neighbor faces: ", false);

		for (int j = 0; j < ring_neighbors.size() - 1; j++)
		{
			Eigen::RowVector3d current_normal = N.row(ring_neighbors[j]);
			Eigen::RowVector3d next_normal = N.row(ring_neighbors[j+1]);

			angles(j) = acos(current_normal.dot(next_normal));
			write_log(loglevel) << "    angle: " << angles(j) << endl;
		}

		int max_index;
		angles.cwiseAbs().maxCoeff(&max_index);
		out_creases(i) = angles[max_index];

		write_log(loglevel) << "assigned angle: " << out_creases(i) << endl << endl;
		//write_log(loglevel) << "assigned angle: " << out_creases(i) << " from: " << endl << angles << endl;
	}
	
	//igl::dihedral_angles is for tetrahedals, not for triangles!
}

void crease::creases_geodesic_neighbors(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double geodesic_distance, Eigen::VectorXd& out_creases)
{
}

//void get_ring_neighborhood(const int vertex_index, const std::vector<std::vector<int>>& adjacency_VF, const int number_rings, ve)
//{
//
//}
