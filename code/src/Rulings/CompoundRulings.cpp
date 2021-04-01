#include "CompoundRulings.h"

#include "GlobalSettings.h"
#include "CurveInterpolation.h"
#include "NormalsHelper.h"
#include "CurveHelper.h"
#include "Utils.h"
#include "Logger.h"

#include <igl/Timer.h>

using namespace std;

void print_bins(const std::vector<std::pair<int, int>>& bins)
{
	std::cout << "bins={";
	for (const auto& bin : bins)
		std::cout << "{" << bin.first << ", " << bin.second << "},";
	
	std::cout << "};" << std::endl;
}


CompoundRulings::CompoundRulings(Mesh& target, Eigen::MatrixXd& geodesic) : AnalyticRulings(target, geodesic)
{
	//sample_step = target.average_edge_length / GlobalSettings::rulings_spline_sample_factor;
	//analytic_rulings = new AnalyticRulings(target, geodesic);
	////discrete_rulings = new DiscreteRulings(target, geodesic);
}

CompoundRulings::~CompoundRulings()
{
	//delete analytic_rulings;
}

void CompoundRulings::compute_curve()
{
	//write_log(4) << "CompoundRulings::compute_curve()" << std::endl;

	AnalyticRulings::compute_curve();
	analytic_ruling_directions = AnalyticRulings::compute_ruling_directions();
	analytic_ruling_lengths = AnalyticRulings::compute_ruling_lengths(analytic_ruling_directions);
	write_log(5) << "analytic_ruling_lengths: " << linebreak << analytic_ruling_lengths << endl;

	//compute_discrete_descriptors();
	discrete_normals = utils::compute_surface_normals_vertices(sampled_curve, target);
}

Eigen::MatrixXd CompoundRulings::compute_ruling_directions()
{
	//write_log(4) << "CompoundRulings::compute_ruling_directions()" << std::endl;


	//if the entire curve is flat, then use binormals as rulings
	bool is_flat_geodesic = curve::is_flat(geodesic);
	if (is_flat_geodesic)
		return get_discrete_rulings();
		//return discrete_binormals;


	//build filters
	std::vector<bool> curvature_filter = build_discrete_curvature_filter();
	std::vector<bool> ruling_lengths_filter = build_ruling_lengths_filter();
	std::vector<bool> tangent_proximity_filter = build_tangent_proximity_filter();

	vector<bool> filter(analytic_ruling_directions.rows(), true);
	for (int i = 0; i < analytic_ruling_directions.rows(); i++)
		if (!curvature_filter[i] || !ruling_lengths_filter[i] || !tangent_proximity_filter[i])
			filter[i] = false;

	write_log(5) << "curvature_filter" << linebreak << list_to_string(curvature_filter) << endl;
	write_log(5) << "ruling_lengths_filter" << linebreak << list_to_string(ruling_lengths_filter) << endl;
	write_log(5) << "tangent_proximity_filter" << linebreak << list_to_string(tangent_proximity_filter) << endl;
	write_log(5) << "filter" << linebreak << list_to_string(filter) << endl;


	//binning
	filter_bins.clear();
	const int min_bin_size = 2;
	int bin_start = filter[0] ? -1 : 0;

	for (int i = 0; i < filter.size(); i++)
	{
		if (filter[i]) //is a valid ruling
		{
			if (bin_start >= 0) //was in a bin which just ended
			{
				filter_bins.push_back(make_pair(bin_start, i - 1));
				bin_start = -1;
			}
			
			continue;
		}

		if (bin_start < 0) //bin just started
			bin_start = i;
	}

	if(bin_start >= 0) //in bin but was not closed because it goes to the end
		filter_bins.push_back(make_pair(bin_start, filter.size() - 1));

	if(filter_bins.size() < 1)
		write_log(4) << "...nothing to filter" << endl;

	//write_log(5) << "filter bins: "; print_bins(filter_bins);


	//interpolate rulings within bins
	Eigen::MatrixXd rulings = analytic_ruling_directions;

	for (const auto& bin : filter_bins)
	{
		int start_bin = bin.first;
		int end_bin = bin.second;
		//write_log(0) << "   interpolate rulings: " << (end_bin - start_bin + 1) << " rulings: " << start_bin << " - " << end_bin << endl;

		Eigen::MatrixXd interpolated = get_interpolated_rulings(start_bin, end_bin, analytic_ruling_directions);
		rulings.block(start_bin, 0, end_bin - start_bin + 1, 3) = interpolated;
	}

	//write_log(0) << "filtered & interpolated rulings: " << linebreak << rulings << endl;
	return rulings;
}


Eigen::MatrixXd CompoundRulings::get_discrete_rulings()
{
	////resample curve to stepsize ~ average edge length
	//double resampled_length;
	//curve::resample_uniformly(geodesic, sample_step, sampled_curve, resampled_length);

	double resampled_length;
	sample_step = curve::compute_length(geodesic) / (sampled_curve.rows()-2);
	curve::resample_uniformly(geodesic, sample_step, sampled_curve, resampled_length);

	compute_discrete_descriptors(sampled_curve, tangents, principal_normals, binormals);
	write_log(5) << "CompoundRulings::get_discrete_rulings tangents: " << linebreak << tangents << endl;
	write_log(5) << "CompoundRulings::get_discrete_rulings principal_normals: " << linebreak << principal_normals << endl;
	write_log(5) << "CompoundRulings::get_discrete_rulings binormals: " << linebreak << binormals << endl;


	//compute curve paramteters
	CurveInterpolation curve(sampled_curve);
	const int rows = curve.k.size();
	curvature = Eigen::Map<Eigen::VectorXd>(curve.k.data(), curve.k.size());
	torsion = Eigen::Map<Eigen::VectorXd>(curve.t.data(), curve.t.size());

	ruling_directions = get_discrete_binormals();
	return ruling_directions;
}

Eigen::VectorXd CompoundRulings::get_discrete_curvature()
{
	//the sample_step is the vertex area since the curve is sampled equidistantly
	double segment = curve::compute_length(sampled_curve) / sampled_curve.rows();

	//use surface normals to compute the curvature
	discrete_normals = utils::compute_surface_normals_vertices(sampled_curve, target);
	//after changing the normals, would need to recompute binormals, but since they aren't used afterwards anymore, we skip this

	/*
	vector<double> k; //TODO use eigen vector right away
	for (int i = 0; i < surface_normals.rows() - 1; i++)
	{
		Eigen::RowVectorXd n1 = surface_normals.row(i);
		Eigen::RowVectorXd n2 = surface_normals.row(i+1);
	
		double angle = acos(clip(n1.dot(n2)));
		double current_curvature = angle / segment;
		k.push_back(current_curvature);
	
		if(i >= surface_normals.rows() - 2) //hack: add the last curvature again to have curvature for each curve vertex
			k.push_back(current_curvature); //TODO do this outside of the loop
	
		//write_log(0) << i << ": k = " << current_curvature << ", angle = " << angle << ", n1 = " << n1 << ", n2 = " << n2 << endl;
	}
	
	//k.push_back(analytic_rulings->curvature(curvature.rows() - 1)); //add the last curvature to have curvature for each curve vertex
	//write_log(0) << "curvature: " << linebreak << list_to_string(k, " ", true, true) << endl;
	
	return to_eigen_vector(k);
	*/
	
	Eigen::VectorXd k(discrete_normals.rows()); 
	double current_curvature = -1.0;

	for (int i = 0; i < discrete_normals.rows() - 1; i++)
	{
		Eigen::RowVectorXd n1 = discrete_normals.row(i);
		Eigen::RowVectorXd n2 = discrete_normals.row(i+1);
	
		double angle = acos(clip(n1.dot(n2)));
		current_curvature = angle / segment;
		k(i) = current_curvature;
		//write_log(0) << i << ": k = " << current_curvature << ", angle = " << angle << ", n1 = " << n1 << ", n2 = " << n2 << endl;
	}

	k(discrete_normals.rows() - 1) = current_curvature;
	//write_log(0) << "discrete_curvature: " << linebreak << k << endl;
	
	return k;
}

std::vector<bool> CompoundRulings::build_discrete_curvature_filter()
{
	//compute curvature based on normals
	discrete_curvature = get_discrete_curvature();

	//label all flat indices
	double const curvature_threshold = GlobalSettings::flat_curvature_threshold / 2; //can be stricter because it is on discrete curvature

	vector<bool> filter(analytic_ruling_directions.rows(), true);
	for (int i = 0; i < discrete_curvature.rows(); i++)
		if (abs(discrete_curvature(i)) <= curvature_threshold)
			filter[i] = false;

	return filter;
}

std::vector<bool> CompoundRulings::build_ruling_lengths_filter()
{
	const double min_length = GlobalSettings::min_ruling_length/2;
	vector<bool> filter(analytic_ruling_directions.rows(), true);
	write_log(5) << "    length filter: min " << min_length << endl;
	//write_log(0) << "    analytic_ruling_lengths" << linebreak << analytic_ruling_lengths << endl;

	for (int i = 0; i < analytic_ruling_lengths.rows(); i++)
		if (abs(analytic_ruling_lengths(i)) < min_length)
			filter[i] = false;

	return filter;
}

std::vector<bool> CompoundRulings::build_tangent_proximity_filter()
{
	vector<bool> filter(analytic_ruling_directions.rows(), true);

	for (int i = 0; i < analytic_ruling_directions.rows(); i++)
	{
		Eigen::RowVector3d tangent = tangents.row(i).normalized();
		Eigen::RowVector3d ruling  = analytic_ruling_directions.row(i).normalized();

		double angle = acos(clip(tangent.dot(ruling)));
		if (angle < GlobalSettings::min_tangent_ruling_deviation)
			filter[i] = false;
	}

	return filter;
}

void CompoundRulings::print_mathematica_data()
{
	std::cout << "(* COMPOUND rulings :: sample_step: " << sample_step << "*)" << std::endl;

	std::cout << "filter="; print_bins(filter_bins);
	std::cout << "nf=" << filter_bins.size() << ";" << std::endl;

	std::cout << "discreteCurvature=" << to_mathematica(discrete_curvature) << std::endl;
	std::cout << "rulingsAnalytic=" << to_mathematica(analytic_ruling_directions) << std::endl;
	std::cout << "rulingsLengthsAnalytic=" << to_mathematica(analytic_ruling_lengths) << std::endl;

	AnalyticRulings::print_mathematica_data();
}

Eigen::MatrixXd CompoundRulings::get_interpolated_rulings(const int start_index, const int end_index, const Eigen::MatrixXd& analytic_ruling_directions)
{
	int number_rulings = end_index - start_index + 1;

	Eigen::RowVector3d start_ruling = analytic_ruling_directions.row(start_index);
	Eigen::RowVector3d end_ruling   = analytic_ruling_directions.row(end_index);
	
	//check if this flat area is at the start or end of the curve, because then the rulings are typically bad so we might want to use the binormal instead
	if(start_index == 0) 
		start_ruling = get_discrete_ruling_at(start_index, analytic_ruling_directions);
	if(end_index == analytic_ruling_directions.rows()-1)
		end_ruling = get_discrete_ruling_at(end_index, analytic_ruling_directions);

	if (start_ruling.dot(end_ruling) < 0)
		end_ruling = -1 * end_ruling;


	//write_log(0) << "interpolate " << number_rulings << " rulings between " << start_index << " - " << end_index << ": " << linebreak 
	//			 << "  start: " << start_ruling << linebreak 
	//			 << "  end:   " << end_ruling << endl;


	Eigen::MatrixXd interpolated(number_rulings, 3);
	for (int i = 0; i < number_rulings; i++)
	{
		double weight = i / (double) (end_index - start_index);
		interpolated.row(i) = ((1-weight) * start_ruling + (weight) * end_ruling).normalized();

		//write_log(0) << i + start_index << ": interpolated ruling = " << interpolated.row(i) << ", weight=" << weight << endl;
	}

	return interpolated;
}

Eigen::RowVector3d CompoundRulings::get_discrete_ruling_at(const int curve_index, const Eigen::MatrixXd& analytic_ruling_directions)
{
	Eigen::RowVector3d tangent = tangents.row(curve_index).normalized();
	Eigen::RowVector3d ruling = analytic_ruling_directions.row(curve_index).normalized();
	
	double angle = acos(clip(tangent.dot(ruling)));
	if (angle > GlobalSettings::min_tangent_ruling_deviation)
		return ruling;

	write_log(5) << "     -->  use binormal as ruling at " << curve_index << ", angle = " << angle << endl;

	Eigen::RowVector3d surface_normal = discrete_normals.row(curve_index);
	Eigen::RowVector3d binormal = tangent.cross(surface_normal);
	return binormal;
}

////Eigen::MatrixXd CompoundRulings::get_interpolated_rulings(const int number_rulings, const Eigen::RowVector3d& start_ruling, const Eigen::RowVector3d& end_ruling)
////{
////	Eigen::MatrixXd interpolated(number_rulings, 3);
////	for (int i = 0; i < number_rulings; i++)
////		interpolated.row(i) = (((number_rulings - i) / number_rulings) * start_ruling + (i / number_rulings) * end_ruling).normalized();
////
////	return interpolated;
////}

//std::vector<std::pair<int, int>> CompoundRulings::get_flat_parts(const Eigen::VectorXd& discrete_curvature)
//{
//	//label all flat indices
//	//double const curvature_threshold = 0.001;
//	//double const curvature_threshold = GlobalSettings::flat_curvature_threshold;
//	double const curvature_threshold = GlobalSettings::flat_curvature_threshold / 2; //can be stricter because it is on discrete curvature
//
//	vector<int> flat_indices;
//	for (int i = 0; i < discrete_curvature.rows(); i++)
//		if (abs(discrete_curvature(i)) <= curvature_threshold)
//			flat_indices.push_back(i);
//
//
//	//find bins of flat incides. "flat_parts" contains start & end indices, both are inclusive.
//	std::vector<std::pair<int, int>> flat_bins;
//	const int min_bin_size = 2;
//	int bin_start = 0;
//
//	for (int i = 1; i < flat_indices.size(); i++)
//	{
//		int diff = flat_indices[i] - flat_indices[i - 1];
//
//		if (i == flat_indices.size() - 1 && diff == 1 && i - bin_start > min_bin_size) //at last position, if still within bucket and bucket is large enough, add
//		{
//			flat_bins.push_back(make_pair(flat_indices[bin_start], flat_indices[i]));
//			continue;
//		}
//
//
//		if (diff == 1) //still in flat area 
//			continue;
//
//		if (i - bin_start > min_bin_size) //flat area done and it is large enough
//			flat_bins.push_back(make_pair(flat_indices[bin_start], flat_indices[i - 1]));
//
//		//write_log(0) << "i=" << i << ", diff=" << diff << ", last bin_start=" << bin_start << endl;
//		bin_start = i;
//	}
//
//	if (LOG_LEVEL <= 4)
//	{
//		write_log(0) << "curvature: " << linebreak << list_to_string(to_std_vector(discrete_curvature), " ", true, true) << endl;
//		write_log(0) << "flat_indices: " << list_to_string(flat_indices) << endl;
//
//		write_log(0) << "flat areas: " << endl;
//		for (const auto& bin : flat_bins)
//			write_log(0) << bin.first << " - " << bin.second << std::endl;
//	}
//
//	return flat_bins;
//}

////Eigen::RowVector3d CompoundRulings::get_discrete_binormal_at(int curve_index)
////{
////	Eigen::RowVector3d tangent = tangents.row(curve_index).normalized();
////	Eigen::RowVector3d surface_normal = surface_normals.row(curve_index);
////
////	Eigen::RowVector3d binormal = tangent.cross(surface_normal);
////	return binormal;
////}

////Eigen::MatrixXd CompoundRulings::get_binormals_as_rulings(int start, int end)
////{
////	//use middle binormal and add to each row as rulings (could be fancier...)
////	int number_flat_rulings = end - start + 1;
////	int mid_index = number_flat_rulings % 2 == 0 ? start + number_flat_rulings / 2 - 1 : start + number_flat_rulings / 2;
////	Eigen::RowVectorXd mid_binormal = discrete_rulings->binormals.row(mid_index);
////
////	Eigen::MatrixXd flat_rulings(number_flat_rulings, 3);
////
////	for (int i = 0; i < number_flat_rulings; i++)
////		flat_rulings.row(i) = mid_binormal;
////
////	return flat_rulings;
////}


void CompoundRulings::compute_discrete_descriptors()
{
	//double resampled_length;
	//double segment = curve::compute_length(geodesic) / (sampled_curve.rows()-2);
	//curve::resample_uniformly(geodesic, segment, discrete_sampled_curve, resampled_length);
	//curve::compute_frames(sampled_curve, discrete_tangents, discrete_normals, discrete_binormals);

	//write_log(0) << linebreak << "compute_discrete_descriptors:: " << linebreak << endl;
	//write_log(0) << "discrete_tangents: "  << linebreak << discrete_tangents << std::endl;
	//write_log(0) << "discrete_normals: "   << linebreak << discrete_normals << std::endl;
	//write_log(0) << "discrete_binormals: " << linebreak << discrete_binormals << std::endl;


	/*
	int number_samples = sampled_curve.rows();
	int columns = 3;

	//compute normals
	bool is_flat_geodesic = curve::is_flat(geodesic);
	if(is_flat_geodesic)
		discrete_normals = utils::compute_surface_normals_faces(sampled_curve, target);
	else
		discrete_normals = utils::compute_surface_normals_vertices(sampled_curve, target);


	//compute tangents & binormals
	discrete_tangents.resize(number_samples, columns);
	discrete_binormals.resize(number_samples, columns);

	for (int i = 0; i < number_samples - 1; i++)
	{
		Eigen::Vector3d tangent = (sampled_curve.row(i + 1) - sampled_curve.row(i)).normalized();
		Eigen::Vector3d normal = discrete_normals.row(i);

		discrete_tangents.row(i) = tangent;
		discrete_binormals.row(i) = tangent.cross(normal).normalized();
	}

	//add last tangent & binormal
	Eigen::Vector3d tangent = (sampled_curve.row(number_samples - 1) - sampled_curve.row(number_samples - 2)).normalized();
	Eigen::Vector3d normal = discrete_normals.row(number_samples - 1);

	discrete_tangents.row(number_samples - 1) = tangent;
	discrete_binormals.row(number_samples - 1) = tangent.cross(normal).normalized();
	*/
}


void CompoundRulings::compute_discrete_descriptors(const Eigen::MatrixXd& sampled_curve, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals)
{
	const int number_samples = sampled_curve.rows();
	const int columns = 3;

	igl::Timer timer;
	double t = timer.getElapsedTime();

	//compute normals (barycentric interpolation)
	out_normals.resize(number_samples, columns);
	out_normals = utils::compute_surface_normals_vertices(sampled_curve, target);
	write_log(5) << "      compute_surface_normals in " << timer.getElapsedTime() - t << std::endl;

	//compute tangents & binormals
	out_tangents.resize(number_samples, columns);
	out_binormals.resize(number_samples, columns);

	for (int i = 0; i < number_samples - 1; i++)
	{
		Eigen::Vector3d tangent = (sampled_curve.row(i + 1) - sampled_curve.row(i)).normalized();
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

Eigen::MatrixXd CompoundRulings::get_discrete_binormals()
{
	//use middle binormal and add to each row as rulings (could be fancier...)
	int mid_index = binormals.rows() % 2 == 0 ? binormals.rows() / 2 - 1 : binormals.rows() / 2;
	Eigen::RowVectorXd mid_binormal = binormals.row(mid_index);
	write_log(5) << "  mid index = " << mid_index << ", mid_binormal: " << mid_binormal << endl;

	Eigen::MatrixXd binormal_rulings(sampled_curve.rows(), 3);

	for (int i = 0; i < sampled_curve.rows(); i++)
		binormal_rulings.row(i) = mid_binormal;

	write_log(5) << "get_discrete_binormals as rulings: " << linebreak << binormal_rulings << std::endl;


	return binormal_rulings;
	/*
		double resampled_length;
		double segment = curve::compute_length(geodesic) / sampled_curve.rows();
		curve::resample_uniformly(geodesic, segment, sampled_curve, resampled_length);

		int number_samples = sampled_curve.rows();
		Eigen::MatrixXd binormals(number_samples, 3);

		// TODO make into function (is from DiscreteRulings::compute_frames())
		for (int i = 0; i < number_samples - 1; i++)
		{
			Eigen::Vector3d tangent = (sampled_curve.row(i + 1) - sampled_curve.row(i)).normalized();
			Eigen::Vector3d normal = surface_normals.row(i);

			binormals.row(i) = tangent.cross(normal).normalized();
		}

		//add last tangent & binormal
		Eigen::Vector3d tangent = (sampled_curve.row(number_samples - 1) - sampled_curve.row(number_samples - 2)).normalized();
		Eigen::Vector3d normal = surface_normals.row(number_samples - 1);
		binormals.row(number_samples - 1) = tangent.cross(normal).normalized();
	*/
}
