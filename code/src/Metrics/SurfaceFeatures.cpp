#include "SurfaceFeatures.h"


#include <igl/principal_curvature.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/cotmatrix.h>
#include <igl/vertex_triangle_adjacency.h>


#include "MeshModel.h"
#include "Utils.h"
#include "MeshController.h"
#include "Logger.h"

void SurfaceFeatures::build(Mesh& mesh)
{
	if (are_built)
		return; 

	compute_principal_k(mesh);
	compute_K(mesh);
	compute_H(mesh);
	compute_creases(mesh);

	are_built = true;
}

void SurfaceFeatures::reset()
{
	are_built = false;
	//TODO check if I need to delete all contents 
}

void SurfaceFeatures::compute_principal_k(Mesh& mesh)
{
	igl::principal_curvature(mesh.V, mesh.F, principal_k_min_direction, principal_k_max_direction, principal_k_min, principal_k_max);
}

void SurfaceFeatures::compute_K(Mesh& mesh)
{
	get_gaussian_curvature(mesh, K);


	Kp = principal_k_min * principal_k_max; //this includes boundary
	Eigen::VectorXd Kp_abs = Kp.cwiseAbs();
	create_lookup(Kp_abs, Kp_abs_lookup);

	////exclude boundary
	//std::vector<int> boundary = meshhelper::boundary_loop_flat(mesh.F);
	//for (int vi : boundary)
	//{
	//	K(vi) = 0;
	//	Kp(vi) = 0;
	//}
	////////


	K_abs = K.cwiseAbs();
	create_lookup(K_abs, K_abs_lookup);
	normalize_vector_entries(K_abs, K_abs_normalized);
}

void SurfaceFeatures::compute_H(Mesh& mesh)
{
	get_mean_curvature(mesh, H);
	create_lookup(H, H_lookup);
	normalize_vector_entries(H, H_normalized);
}

void SurfaceFeatures::compute_creases(Mesh& mesh)
{
	dihedral_triangle_angles(mesh, crease_vertices);
	create_lookup(crease_vertices, crease_vertices_lookup);
	normalize_vector_entries(crease_vertices, crease_vertices_normalized);
}

void SurfaceFeatures::normalize_vector_entries(const Eigen::VectorXd& curvature, Eigen::VectorXd& out_normalized_curvature)
{
	double min = curvature.minCoeff();
	double max = curvature.maxCoeff();
	double range = max - min;
	out_normalized_curvature = (curvature - min * Eigen::VectorXd::Ones(curvature.rows())) / range;
}

void SurfaceFeatures::get_gaussian_curvature(Mesh& target, Eigen::VectorXd& out_K)
{
	//get gaussian curvature
	igl::gaussian_curvature(target.V, target.F, out_K); // Compute integral of Gaussian curvature

	// Compute mass matrix
	Eigen::SparseMatrix<double> M, Minv;
	igl::massmatrix(target.V, target.F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);

	out_K = (Minv*out_K).eval(); // Divide by area to get integral average

	
	
	
	Eigen::VectorXd ones(target.V.rows());
	ones.setOnes();
	Eigen::VectorXd V_mass = M * ones;
	write_log(0) << "Mass mean = " << V_mass.mean() << std::endl;
}

void SurfaceFeatures::get_mean_curvature(Mesh& target, Eigen::VectorXd& out_H)
{
	Eigen::MatrixXd HN;
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(target.V, target.F, L);
	igl::massmatrix(target.V, target.F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);

	HN = -Minv * (L * target.V);
	out_H = HN.rowwise().norm(); //up to sign
}

void SurfaceFeatures::dihedral_triangle_angles(const Mesh& mesh, Eigen::VectorXd& out_creases)
{
	const int loglevel = 6;

	out_creases.resize(mesh.V.rows());

	Eigen::MatrixXd N;
	igl::per_face_normals(mesh.V, mesh.F, N);
	
	std::vector<std::vector<int>> VFi_unused;
	std::vector<std::vector<int>> adjacency_VF;
	igl::vertex_triangle_adjacency(mesh.V, mesh.F, adjacency_VF, VFi_unused);

	//out_normal = N.row(closest_faces(0)) + point;

	for (int i = 0; i < mesh.V.rows(); i++)
	{
		std::vector<int> ring_neighbors = adjacency_VF[i];
		if (ring_neighbors.size() < 2)
			continue;

		Eigen::VectorXd angles(ring_neighbors.size() - 1);

		write_log(loglevel) << "[" << i << "]: "; log_list(loglevel, ring_neighbors, "neighbor faces: ", false);

		for (int j = 0; j < ring_neighbors.size() - 1; j++)
		{
			Eigen::RowVector3d current_normal = N.row(ring_neighbors[j]);
			Eigen::RowVector3d next_normal = N.row(ring_neighbors[j + 1]);

			angles(j) = acos(current_normal.dot(next_normal));
			write_log(loglevel) << "    angle: " << angles(j) << std::endl;
		}

		int max_index;
		angles.cwiseAbs().maxCoeff(&max_index);
		out_creases(i) = angles[max_index];

		write_log(loglevel) << "assigned angle: " << out_creases(i) << std::endl << std::endl;
	}

	//igl::dihedral_angles is for tetrahedals, not for triangles!
}
