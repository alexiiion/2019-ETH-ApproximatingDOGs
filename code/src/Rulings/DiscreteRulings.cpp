#include "DiscreteRulings.h"

#include "NormalsHelper.h"
#include "CurveHelper.h"
#include "CurveInterpolation.h"
#include "CoordinateConverter.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <igl/Timer.h>

DiscreteRulings::DiscreteRulings(Mesh& target, Eigen::MatrixXd& geodesic) : Rulings(target, geodesic)
{
	//sample_step = target.average_edge_length * 0.5;
	sample_step = target.average_edge_length / GlobalSettings::rulings_spline_sample_factor;
}

DiscreteRulings::~DiscreteRulings() { }

void DiscreteRulings::compute_curve()
{
	//write_log(0) << "DiscreteRulings::compute_curve()" << std::endl;

	//resample curve to stepsize ~ average edge length
	double resampled_length;
	curve::resample_uniformly(geodesic, sample_step, sampled_curve, resampled_length);

	compute_frames(sampled_curve, tangents, principal_normals, binormals);

	//compute curve paramteters
	CurveInterpolation curve(sampled_curve);
	const int rows = curve.k.size();
	curvature = Eigen::Map<Eigen::VectorXd>(curve.k.data(), curve.k.size());
	torsion = Eigen::Map<Eigen::VectorXd>(curve.t.data(), curve.t.size());

	is_flat = curve::is_flat(curve.k);
}

Eigen::MatrixXd DiscreteRulings::compute_ruling_directions()
{
	write_log(0) << "DiscreteRulings::compute_ruling_directions()" << std::endl;

	//if flat: use middle binormal and add to each row as rulings (could be fancier...)
	if (is_flat)
		return get_binormals_as_rulings();


	Eigen::MatrixXd rulings(sampled_curve.rows(), 3);

	for (int i = 0; i < sampled_curve.rows() - 1; i++)
	{
		Eigen::Vector3d normal_current = principal_normals.row(i);
		Eigen::Vector3d normal_next = principal_normals.row(i + 1);
		rulings.row(i) = normal_current.cross(normal_next).normalized();
	}

	return rulings;
}

Eigen::MatrixXd DiscreteRulings::compute_surface_normals(const Eigen::MatrixXd& sampled_curve)
{
	return utils::compute_surface_normals_vertices(sampled_curve, target);

	/*
	Eigen::MatrixXd normals(sampled_curve.rows(), 3);

	Eigen::MatrixXi bary_indices;
	Eigen::MatrixXd bary_weigths;
	CoordinateConverter::barycentric_coords_from_points(sampled_curve, target.V, target.F, bary_indices, bary_weigths);

	for (int i = 0; i < sampled_curve.rows(); i++)
	{
		Eigen::Vector3d normal(0, 0, 0);

		for (int bary_i = 0; bary_i < 3; bary_i++)
		{
			double weight = bary_weigths(i, bary_i);
			int vertex_index = bary_indices(i, bary_i);

			normal += weight * target.normals_vertices.row(vertex_index);
		}

		normals.row(i) = normal.normalized();
	}

	return normals;
	*/
}

void DiscreteRulings::compute_frames(const Eigen::MatrixXd& sampled_curve, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals)
{
	const int number_samples = sampled_curve.rows();
	const int columns = 3;

	igl::Timer timer;
	double t = timer.getElapsedTime();

	//compute normals (barycentric interpolation)
	out_normals.resize(number_samples, columns);
	out_normals = compute_surface_normals(sampled_curve);
	write_log(5) << "      compute_surface_normals in " << timer.getElapsedTime() - t << std::endl;

	//compute tangents & binormals
	out_tangents.resize(number_samples, columns);
	out_binormals.resize(number_samples, columns);
	
	for (int i = 0; i < number_samples-1; i++)
	{
		Eigen::Vector3d tangent = (sampled_curve.row(i+1) - sampled_curve.row(i)).normalized();
		Eigen::Vector3d normal = out_normals.row(i);

		out_binormals.row(i) = tangent.cross(normal).normalized();
		out_tangents.row(i) = tangent;
	}

	//add last tangent & binormal
	Eigen::Vector3d tangent = (sampled_curve.row(number_samples - 1) - sampled_curve.row(number_samples - 2)).normalized();
	Eigen::Vector3d normal = out_normals.row(number_samples - 1);
	out_binormals.row(number_samples - 1) = tangent.cross(normal).normalized();
	out_tangents.row(number_samples - 1) = tangent;
}

Eigen::MatrixXd DiscreteRulings::get_binormals_as_rulings()
{
	//use middle binormal and add to each row as rulings (could be fancier...)
	int mid_index = binormals.rows() % 2 == 0 ? binormals.rows() / 2 - 1 : binormals.rows() / 2;
	Eigen::RowVectorXd mid_binormal = binormals.row(mid_index);

	Eigen::MatrixXd rulings(sampled_curve.rows(), 3);

	for (int i = 0; i < sampled_curve.rows(); i++)
		rulings.row(i) = mid_binormal;

	return rulings;
}

void DiscreteRulings::print_mathematica_data()
{
	std::cout << "(* DISCRETE rulings :: sample_step: " << sample_step << "*)" << std::endl;

	auto SN = utils::compute_surface_normals_vertices(sampled_curve, target);
	std::cout << "surfaceNormals=" << to_mathematica(SN) << std::endl;

	Rulings::print_mathematica_data();
}
