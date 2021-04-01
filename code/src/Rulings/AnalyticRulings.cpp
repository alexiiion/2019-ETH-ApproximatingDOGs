#include "AnalyticRulings.h"

#include "GlobalSettings.h"
#include "NormalsHelper.h"
#include "CurveHelper.h"
#include "Logger.h"

AnalyticRulings::AnalyticRulings(Mesh& target, Eigen::MatrixXd& geodesic) : Rulings(target, geodesic)
{
	//double length = curve::compute_length(geodesic);
	//sample_step = target.average_edge_length * 0.5 / length;
	//smoothness = 0.05;
	//degree = 5;

	double length = curve::compute_length(geodesic);
	sample_step = target.average_edge_length / GlobalSettings::rulings_spline_sample_factor / length;
	smoothness = GlobalSettings::rulings_spline_smoothness;
	degree = GlobalSettings::rulings_spline_degree;
}

AnalyticRulings::~AnalyticRulings() 
{ 
	delete spline; 
}

void AnalyticRulings::compute_curve()
{
	//write_log(0) << "AnalyticRulings::compute_curve()" << std::endl;

	//double length = curve::compute_length(geodesic);
	//double segment_length = length / geodesic.rows();
	//segment_length *= GlobalSettings::rulings_spline_sample_factor;

	//Eigen::MatrixXd sparse_geodesic;
	//double resampled_length;
	//curve::resample_uniformly(geodesic, segment_length, sparse_geodesic, resampled_length);
	//
	//geodesic = sparse_geodesic;

	//create spline
	spline = new Spline(geodesic, smoothness, degree);
	sampled_curve = spline->sample_curve(sample_step);

	derivative_1 = spline->sample_derivative(1, sample_step);
	derivative_2 = spline->sample_derivative(2, sample_step);
	derivative_3 = spline->sample_derivative(3, sample_step);


	//compute frames
	const int number_samples = sampled_curve.rows();

	tangents.resize(number_samples, 3);
	for (int i = 0; i < number_samples; i++)
		tangents.row(i) = derivative_1.row(i) / derivative_1.row(i).norm();

	principal_normals.resize(number_samples, 3);
	for (int i = 0; i < number_samples; i++)
		principal_normals.row(i) = (derivative_2.row(i) - derivative_2.row(i).dot(tangents.row(i)) * tangents.row(i)).normalized();

	binormals = cross_rowwise(tangents, principal_normals);


	Eigen::MatrixXd der1_der2_cross = cross_rowwise(derivative_1, derivative_2);
	Eigen::VectorXd der1_der2_cross_norm = der1_der2_cross.rowwise().norm();
	Eigen::VectorXd der1_norm = derivative_1.rowwise().norm();
	Eigen::VectorXd der1_norm_pow3 = der1_norm.array().pow(3);


	//compute curve paramteters
	curvature.resize(der1_der2_cross_norm.rows());
	for (int i = 0; i < der1_der2_cross_norm.rows(); i++)
		curvature(i) = der1_der2_cross_norm(i) / der1_norm_pow3(i);


	Eigen::VectorXd torsion_up(der1_der2_cross.rows());
	for (int i = 0; i < der1_der2_cross.rows(); i++)
		torsion_up(i) = der1_der2_cross.row(i).dot(derivative_3.row(i));

	Eigen::VectorXd torsion_down(der1_der2_cross.rows());
	for (int i = 0; i < der1_der2_cross.rows(); i++)
		torsion_down(i) = der1_der2_cross.row(i).dot(der1_der2_cross.row(i));

	torsion.resize(torsion_up.rows());
	for (int i = 0; i < torsion_up.rows(); i++)
		torsion(i) = torsion_up(i) / torsion_down(i);
}

Eigen::MatrixXd AnalyticRulings::compute_ruling_directions()
{
	//write_log(0) << "AnalyticRulings::compute_ruling_directions()" << std::endl;

	//darboux
	Eigen::MatrixXd rulings(torsion.rows(), 3);
	for (int i = 0; i < torsion.rows(); i++)
		rulings.row(i) = (torsion(i)*tangents.row(i) + curvature(i)*binormals.row(i)).normalized();

	return rulings;
}

void AnalyticRulings::print_mathematica_data()
{
	std::cout << "(* ANALYTIC rulings :: smooth: " << smoothness << ", degree: " << degree << ", sample_step: " << sample_step << "*)" << std::endl;
	
	std::cout << "spline1Derivative=" << to_mathematica(derivative_1) << std::endl;
	std::cout << "spline2Derivative=" << to_mathematica(derivative_2) << std::endl;
	std::cout << "spline3Derivative=" << to_mathematica(derivative_3) << std::endl;

	auto SN = utils::compute_surface_normals_vertices(sampled_curve, target);
	std::cout << "surfaceNormals="    << to_mathematica(SN) << std::endl;

	Rulings::print_mathematica_data();
}

Eigen::MatrixXd AnalyticRulings::cross_rowwise(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2)
{
	//if(m1.rows() != m2.rows())

	Eigen::MatrixXd result(m1.rows(), 3);
	for (int i = 0; i < m1.rows(); i++)
	{
		Eigen::Vector3d v1 = m1.row(i);
		Eigen::Vector3d v2 = m2.row(i);
		result.row(i) = v1.cross(v2);
	}

	return result;
}
