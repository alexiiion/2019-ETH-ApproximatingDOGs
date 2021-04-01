#include "GaussianCurvature.h"

#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

//#include "MeshModel.h"

//void GaussianCurvature::update(bool exclude_border = true)
//{
//	//TODO implement
//}

void GaussCurvature::reset()
{
	K.resize(mesh.V.rows());
}


void GaussCurvature::update(Eigen::VectorXd vertex_mask)
{
	K.resize(mesh.V.rows());
	igl::gaussian_curvature(mesh.V, mesh.F, K); // Compute integral of Gaussian curvature

	// Compute mass matrix
	Eigen::SparseMatrix<double> M, Minv;
	igl::massmatrix(mesh.V, mesh.F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);

	K = (Minv*K).eval(); // Divide by area to get integral average
	K = K.array() * vertex_mask.array();
}

Eigen::VectorXd GaussCurvature::get()
{
	return K;
}

double GaussCurvature::average()
{
	if (K.rows() < 1)
		return -1;

	return K.cwiseAbs().sum() / mesh.V.rows();
}

double GaussCurvature::min()
{
	if (K.rows() < 1)
		return -1e6;

	return K.minCoeff();
}

double GaussCurvature::max()
{
	if (K.rows() < 1)
		return 1e6;

	return K.maxCoeff();
}

//double GaussCurvature::min()
//{
//
//}
//
//double GaussCurvature::max()
//{
//
//}


/*
	Eigen::VectorXd a(3);
	a << 1, 2, 3;
	Eigen::VectorXd b(3);
	b << 1, 0, 1;
	Eigen::VectorXd elementwise_product = a.array() * b.array();
	write_log(4) << endl << "      elementwise_product: " << elementwise_product.transpose() << endl;
*/