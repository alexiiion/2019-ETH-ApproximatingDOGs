#include "CurveHelper.h"

#include "CurveInterpolation.h"

#include "CoordinateConverter.h"
#include "GlobalSettings.h"
#include "Utils.h"
#include "NormalsHelper.h"
#include "Logger.h"


using namespace std;

Eigen::VectorXd laplace_on_angles(const Eigen::MatrixXd points, const int i)
{
	Eigen::VectorXd direction_back = (points.row(i - 1) - points.row(i)).normalized();
	Eigen::VectorXd direction_front = (points.row(i + 1) - points.row(i)).normalized();

	//double dot_back = (points.row(i).dot(points.row(i - 1)));
	//double dot_front = (points.row(i).dot(points.row(i + 1)));
	double angle_back = acos(points.row(i - 1).dot(points.row(i)));
	double angle_front = acos(points.row(i + 1).dot(points.row(i)));

	Eigen::VectorXd Lp = direction_back  * angle_back  / (angle_back + angle_front) +
						 direction_front * angle_front / (angle_back + angle_front);

	return Lp;
}

Eigen::MatrixXd curve::remove_short_segments(const Eigen::MatrixXd& path, const double min_length)
{
	const int loglevel = 5;

	const int number_points = path.rows();
	Eigen::MatrixXd path_processed(number_points, 3);

	if (number_points < 2)
		return path_processed;

	path_processed.row(0) = path.row(0); //add first point of geodesic
	int count = 1;

	double min_squared_length = min_length * min_length;
	//double length = 0.0;

	for (int j = 1; j < number_points - 1; j++)
	{
		double segment_length = (path_processed.row(count-1) - path.row(j)).squaredNorm();
		if (segment_length < min_length)
			continue;
		
		path_processed.row(count) = path.row(j);
		//length += (path_processed.row(count) - path_processed.row(count - 1)).norm();

		count++;
		write_log(loglevel + 1) << "  --> add vector!" << endl;
	}

	path_processed.row(count) = path.row(number_points - 1); //add last point of geodesic
	path_processed.conservativeResize(count + 1, Eigen::NoChange);

	return path_processed;
}

Eigen::MatrixXd curve::remove_straight_segments(const Eigen::MatrixXd& path, const double epsilon_angle, const double min_length, const double max_length)
{
	const int loglevel = 5;

	const int number_points = path.rows();
	Eigen::MatrixXd path_processed(number_points, 3);

	if (number_points < 2)
		return path_processed;

	path_processed.row(0) = path.row(0); //add first point of geodesic
	int count = 1;

	for (int j = 1; j < number_points - 1; j++)
	{
		double segment_length = (path.row(j + 1) - path.row(j - 1)).squaredNorm();
		if (segment_length < min_length)
			continue;

		double gap_length = (path_processed.row(count - 1) - path.row(j)).squaredNorm();
		if (gap_length >= max_length)
		{
			path_processed.row(count) = path.row(j);
			count++;

			write_log(loglevel + 1) << "  --> add vector!" << endl;
			continue;
		}

		Eigen::RowVector3d e1 = (path.row(j) - path.row(j - 1)).normalized();
		Eigen::RowVector3d e2 = (path.row(j + 1) - path.row(j)).normalized();

		double dot = e1.dot(e2);
		double cos_angle = clip(dot, -1, 1);

		if (cos_angle == 1 || cos_angle == -1)
			continue;

		double angle = acos(cos_angle);

		write_log(loglevel + 1) << "  point [" << j << "]: " << endl;
		write_log(loglevel + 1) << "    cos(theta) = " << cos_angle << endl;

		if (abs(angle) < epsilon_angle)
			continue;

		path_processed.row(count) = path.row(j);
		count++;

		write_log(loglevel + 1) << "  --> add vector!" << endl;
	}

	path_processed.row(count) = path.row(number_points - 1); //add last point of geodesic
	path_processed.conservativeResize(count + 1, Eigen::NoChange);

	return path_processed;
}

double curve::compute_length(const Eigen::MatrixXd& path)
{
	const int loglevel = 6;

	const int number_points = path.rows();
	double length = 0.0;

	write_log(loglevel) << "segment lengths:  " << endl;
	for (int i = 0; i < number_points - 1; i++)
	{
		double segment_length = (path.row(i + 1) - path.row(i)).norm();
		length += segment_length;

		write_log(loglevel) << i << ": " << segment_length << endl;
	}
	return length;
}

bool curve::get_frame_at(const int i, const Eigen::MatrixXd& path, Eigen::RowVector3d& out_T, Eigen::RowVector3d& out_N, Eigen::RowVector3d& out_B)
{
	if (i < 1 || i > path.rows() - 2)
		return false;

	Eigen::RowVector3d ef = (path.row(i + 1) - path.row(i)).normalized();
	Eigen::RowVector3d eb = (path.row(i - 1) - path.row(i)).normalized();
	//write_log(0) << i << ": ef = " << ef << ", eb = " << eb << endl;

	out_T = (ef - eb).normalized();

	if ((out_T - ef).norm() < 1e-6) //if this segment is perfectly flat, then construct normal
	{
		//write_log(0) << "  flat segment: T = " << out_T << endl;

		Eigen::RowVector3d e0 = abs(out_T(0)) > .9 ? Eigen::RowVector3d(0, 1, 0) : Eigen::RowVector3d(1, 0, 0);
		e0 -= e0.dot(out_T) * out_T;
		e0.normalize();
		Eigen::RowVector3d e1 = out_T.cross(e0);
		//write_log(0) << "  e0 = " << e0 << ", e1 = " << e1 << endl;

		out_N = e1;
		out_B = e0;
		return true;
	}


	out_N = (ef + eb).normalized();
	if (out_N.hasNaN())
		write_log(0) << "  N has nan: " << out_N << endl;

	out_B = out_T.cross(out_N);
	
	return true;
}

Eigen::Matrix3d curve::get_frame_at(int i, const Eigen::MatrixXd& path)
{
	Eigen::RowVector3d T, N, B;
	bool success = curve::get_frame_at(i, path, T, N, B);

	if (!success)
		return Eigen::Matrix3d();

	Eigen::Matrix3d frame; 
	frame.col(0) = T; 
	frame.col(1) = N; 
	frame.col(2) = B;
	
	return frame;
}

void curve::compute_frames(const Eigen::MatrixXd& path, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals)
{
	const int number_samples = path.rows();
	const int columns = 3;

	out_tangents.resize(number_samples, columns);
	out_normals.resize(number_samples, columns);
	out_binormals.resize(number_samples, columns);

	for (int i = 1; i < number_samples - 1; i++)
	{
		Eigen::RowVector3d T, N, B;
		bool success = curve::get_frame_at(i, path, T, N, B);
		if (!success)
			return;

		out_tangents.row(i)  = T;
		out_normals.row(i)   = N;
		out_binormals.row(i) = B;
	}

	//repeat at beginning
	out_tangents.row(0)  = out_tangents.row(1);
	out_normals.row(0)   = out_normals.row(1);
	out_binormals.row(0) = out_binormals.row(1);

	//repeat at end
	out_tangents.row(number_samples-1)  = out_tangents.row(number_samples-2);
	out_normals.row(number_samples-1)   = out_normals.row(number_samples-2);
	out_binormals.row(number_samples-1) = out_binormals.row(number_samples-2);
}

void curve::compute_surface_frames(const Eigen::MatrixXd& sampled_curve, const Mesh& mesh, Eigen::MatrixXd& out_tangents, Eigen::MatrixXd& out_normals, Eigen::MatrixXd& out_binormals)
{
	const int number_samples = sampled_curve.rows();
	const int columns = 3;

	//igl::Timer timer;
	//double t = timer.getElapsedTime();

	//compute normals (barycentric interpolation)
	out_normals.resize(number_samples, columns);
	out_normals = utils::compute_surface_normals_vertices(sampled_curve, mesh);
	//write_log(5) << "      compute_surface_normals in " << timer.getElapsedTime() - t << std::endl;

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

std::vector<double> curve::compute_normal_deviations(const Eigen::MatrixXd& path, const Mesh& target)
{
	//get normal deviations
	vector<double> deviations;

	for (int i = 1; i < path.rows(); i++)
	{
		auto nb = get_surface_normal(path.row(i - 1), target.V, target.F, target.normals_vertices, target.adjacency_VF);
		auto nf = get_surface_normal(path.row(i), target.V, target.F, target.normals_vertices, target.adjacency_VF);

		deviations.push_back(angle(nb, nf));
	}

	return deviations;
}

//std::vector<double> to_std_vector(const Eigen::VectorXd& eigen_vector)
//{
//	vector<double> std_vector(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
//	return std_vector;
//}
//
//Eigen::VectorXd to_eigen_vector(std::vector<double>& std_vector)
//{
//	Eigen::VectorXd eigen_vector = Eigen::Map<Eigen::VectorXd>(std_vector.data(), std_vector.size());
//	return eigen_vector;
//}

bool curve::is_flat(const Eigen::MatrixXd& path)
{
	return curve::is_flat(path, GlobalSettings::flat_curvature_threshold);
}

bool curve::is_flat(const Eigen::MatrixXd& path, const double flat_threshold)
{
	const int loglevel = 5; 

	if (path.rows() < 4)
		return false;

	CurveInterpolation curve(path);
	write_log(loglevel) << "curvature = " << to_mathematica(curve.k) << endl;

	return curve::is_flat(curve.k, flat_threshold);
}

bool curve::is_flat(std::vector<double>& path_curvature)
{
	return curve::is_flat(path_curvature, GlobalSettings::flat_curvature_threshold);
}

bool curve::is_flat(std::vector<double>& path_curvature, const double flat_threshold)
{
	const int loglevel = 6; 

	if (path_curvature.size() < 1)
		return false;


	//use absolute curvature
	for (auto& k : path_curvature)
		k = abs(k); 

	write_log(loglevel) << "curvatureAbs = " << to_mathematica(path_curvature) << endl;


	double median;
	double q1;
	double q3;
	compute_median(path_curvature, median, q1, q3);

	double iqr = q3 - q1;
	double lower_bound = q1 - 1.5*iqr;
	double upper_bound = q3 + 1.5*iqr;

	vector<double> outlier_cleaned_curvature;
	for (double k : path_curvature)
	{
		if (k >= lower_bound && k <= upper_bound)
			outlier_cleaned_curvature.push_back(k);
	}

	double max_k = *std::max_element(outlier_cleaned_curvature.begin(), outlier_cleaned_curvature.end());
	write_log(loglevel) << "outlierCleanedCurvature = " << to_mathematica(outlier_cleaned_curvature);
	write_log(loglevel-1) << "max_k = " << max_k << " --> is flat = " << boolalpha << (max_k <= flat_threshold) << endl;

	if (max_k <= flat_threshold)
		return true;

	return false;
}

double curve::compute_curvature(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2)
{
	Eigen::RowVector3d e1 = v1.normalized();
	Eigen::RowVector3d e2 = v2.normalized();

	double dot = e1.dot(e2);
	double cos_angle = clip(dot, -1, 1);

	if (cos_angle == 1 || cos_angle == -1)
		return 0.0;

	double angle = acos(cos_angle);
	double hypothenuse = (v2 - v1).norm();
	double k = 2 * sin(angle) / hypothenuse;

	if (hypothenuse < 1e-3) //too short hypothenuse creates huge curvature. but with dense meshes this happens... 
		return 0.0;

	return k;
}

Eigen::VectorXd curve::compute_curvature(const Eigen::MatrixXd& path)
{
	CurveInterpolation curve(path);
	return to_eigen_vector(curve.k);

	/*
	const int loglevel = 6;

	const int number_points = path.rows();
	Eigen::VectorXd curvature(number_points - 2);
	//double curvature = 0.0;

	write_log(loglevel) << "curvatures:  " << endl;
	for (int i = 0; i < number_points - 2; i++)
	{
		Eigen::RowVector3d e1 = (path.row(i + 1) - path.row(i)).normalized();
		Eigen::RowVector3d e2 = (path.row(i + 2) - path.row(i + 1)).normalized();

		double dot = e1.dot(e2);
		double cos_angle = clip(dot, -1, 1);

		if (cos_angle == 1 || cos_angle == -1)
			continue;

		double angle = acos(cos_angle);
		double hypothenuse = (path.row(i + 2) - path.row(i)).norm();
		double current_curvature = 2 * sin(angle) / hypothenuse;

		if (hypothenuse < 1e-3) //too short hypothenuse creates huge curvature. but with dense meshes this happens... 
			current_curvature = 0;

		curvature(i) = current_curvature;
		//curvature += current_curvature;
		write_log(loglevel) << current_curvature << ", ";
	}
	write_log(loglevel) << endl;

	return curvature;
	*/
}

Eigen::MatrixXd curve::smooth_curve(const Eigen::MatrixXd& path, const int iterations, const double smooth_rate, const double inflate_rate)
{
	const int n = path.rows();
	Eigen::MatrixXd path_smoothed = path;

	for (int t = 0; t < iterations; t++)
	{
		//write_log(0) << "smoothing:: " << t << endl;

		for (int i = 1; i < n - 1; i++)
		{
			Eigen::RowVector3d eb = path_smoothed.row(i - 1) - path_smoothed.row(i);
			Eigen::RowVector3d ef = path_smoothed.row(i + 1) - path_smoothed.row(i);

			Eigen::RowVector3d Lp = eb * 0.5 + ef * 0.5;
			Eigen::RowVector3d point = path_smoothed.row(i) + smooth_rate * Lp;

			//eb = path_smoothed.row(i - 1) - point;
			//ef = path_smoothed.row(i + 1) - point;
			point = point + inflate_rate * Lp;

			//write_log(0) << "  " << i << ": difference = " << point - path_smoothed.row(i) << endl;
			path_smoothed.row(i) = point;
		}
	}

	return path_smoothed;
}

Eigen::MatrixXd curve::smooth_normals(const Eigen::MatrixXd& normals, const int iterations, const double smooth_rate, const double inflate_rate)
{
	if (abs(smooth_rate) < 1e-6)
		return normals;

	const int n = normals.rows();
	Eigen::MatrixXd normals_smoothed = normals;

	for (int t = 0; t < iterations; t++)
	{
		//write_log(0) << "smoothing:: " << t << endl;

		for (int i = 1; i < n - 1; i++)
		{
			Eigen::VectorXd Lp = laplace_on_angles(normals_smoothed, i);
			Eigen::VectorXd point = (normals_smoothed.row(i).transpose() + smooth_rate * Lp).normalized();
			normals_smoothed.row(i) = point;

			Lp = laplace_on_angles(normals_smoothed, i);
			point = (normals_smoothed.row(i).transpose() + inflate_rate * Lp).normalized();
			normals_smoothed.row(i) = point;



			/*
			Eigen::RowVector3d direction_back = (normals_smoothed.row(i - 1) - normals_smoothed.row(i)).normalized();
			Eigen::RowVector3d direction_front = (normals_smoothed.row(i + 1) - normals_smoothed.row(i)).normalized();
			
			double dot_back = (normals_smoothed.row(i).dot(normals_smoothed.row(i - 1)));
			double dot_front = (normals_smoothed.row(i).dot(normals_smoothed.row(i + 1)));
			double angle_back = acos(normals_smoothed.row(i - 1).dot(normals_smoothed.row(i)));
			double angle_front = acos(normals_smoothed.row(i + 1).dot(normals_smoothed.row(i)));

			//write_log(0) << "  " << i << ": dot_back = " << dot_back << ", dot_front: " << dot_front << endl;


			Eigen::RowVector3d Lp = direction_back  * angle_back  / (angle_back + angle_front) + 
									direction_front * angle_front / (angle_back + angle_front);
			Eigen::RowVector3d point = (normals_smoothed.row(i) + smooth_rate * Lp).normalized();

			//eb = normals_smoothed.row(i - 1) - point;
			//ef = normals_smoothed.row(i + 1) - point;
			//Eigen::RowVector3d inflated_point = point + inflate_rate * Lp;
			//point = (point + inflate_rate * Lp).normalized();

			//write_log(0) << "  " << i << ": prev. normal = (" << normals_smoothed.row(i) << "), new normal = (" << point << "), difference = (" << point - normals_smoothed.row(i) << ")" << endl;
			normals_smoothed.row(i) = point;
			*/
		}
	}

	return normals_smoothed;
}

void curve::resample_uniformly(const Eigen::MatrixXd& original_path, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	resample_uniformly(original_path, 1.0, out_sampled, out_sampled_length);
}

/*
void curve::resample_uniformly(const Eigen::MatrixXd& original_path, const double segment_length, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	const int loglevel = 4;

	const int number_points = original_path.rows();
	write_log(loglevel) << "--> original_path: " << endl << original_path << endl;

	if (!original_path.allFinite())
	{
		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesic is not finite!" << endl << endl;
		write_log(1) << "--> original_path: " << endl << original_path << endl << endl;
		exit(EXIT_FAILURE);
	}

	out_sampled.resize(number_points, 3);
	out_sampled_length = 0;

	Eigen::RowVectorXd current_point = original_path.row(0);
	out_sampled.row(0) = current_point;
	write_log(loglevel) << "  adding point at " << 0 << ": " << current_point << endl << endl;
	
	int number_points_added = 1;
	int segment_source_index = 0;
	int segment_target_index = 0;
	double current_length = 0.0;

	while (segment_target_index < number_points)
	{
		while (current_length < segment_length && segment_target_index < number_points - 1)
		{
			segment_target_index++;
			current_length = (original_path.row(segment_target_index) - current_point).norm(); //could be squared norm
		}

		double previous_segment_length = (current_point - original_path.row(segment_target_index - 1)).norm();
		double target_length = sqrt(segment_length * segment_length - previous_segment_length * previous_segment_length);
		
		auto target_segment_vector = original_path.row(segment_target_index) - original_path.row(segment_target_index - 1);
		double scale = target_length / target_segment_vector.norm();

		Eigen::RowVectorXd point = current_point + target_segment_vector * scale;

		double check_length = (point - current_point).norm();
		double check_difference = abs(check_length - segment_length);
		//assert(check_difference < 1e-3);
		write_log(loglevel) << "  -> CHECK: difference =  " << check_difference << ", lenght = " << check_length << " (" << (point - out_sampled.row(number_points_added - 1)).norm() << ")" << endl;
		write_log(loglevel) << "  adding point at " << number_points_added << ": " << point << endl << endl;

		out_sampled_length += check_length;
		out_sampled.row(number_points_added) = point;
		number_points_added++;

		current_point = point;
	}
}
*/

void curve::resample_uniformly(const Eigen::MatrixXd& original_path, const double segment_length, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	const int loglevel = 6;

	const int number_points = original_path.rows();
	write_log(loglevel) << "--> original_path: " << endl << original_path << endl;

	if (number_points < 3)
	{
		out_sampled = original_path;
		out_sampled_length = curve::compute_length(original_path);
		return;
	}

	if (!original_path.allFinite())
	{
		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesic is not finite!" << endl << endl;
		write_log(1) << "--> original_path: " << endl << original_path << endl << endl;
		exit(EXIT_FAILURE);
	}

	out_sampled.resize(number_points, 3);
	out_sampled_length = 0;

	Eigen::RowVectorXd current_point = original_path.row(0);
	out_sampled.row(0) = current_point;
	write_log(loglevel) << "  adding point at " << 0 << ": " << current_point << endl << endl;

	//const int edge_length = 1.0;
	int segment_start_index = 0;
	int segment_end_index = 0;
	double current_length = 0.0;
	int number_points_added = 1;

	//while (segment_end_index < number_points - 1)
	while (segment_end_index < number_points)
		//while (true)
	{
		if (number_points_added >= out_sampled.rows() - 1)
			out_sampled.conservativeResize(number_points_added + number_points, Eigen::NoChange);

		Eigen::RowVectorXd next_point;
		current_length = 0.0;
		segment_end_index = segment_start_index;

		while (current_length < segment_length && segment_end_index < number_points - 1)
		{
			segment_end_index++;
			next_point = original_path.row(segment_end_index);
			current_length = (next_point - current_point).norm(); //could be squared norm
		}

		if (segment_end_index - segment_start_index == 1 && current_length >= segment_length)
		{
			write_log(loglevel) << "in ONE segment: " << current_length << endl;
			write_log(loglevel) << "  segment (" << segment_start_index << " - " << segment_end_index << ")" << endl;
			Eigen::RowVectorXd segment_direction = (original_path.row(segment_end_index) - original_path.row(segment_start_index)).normalized();
			double scale = segment_length / 1.0;
			Eigen::RowVectorXd point = current_point + segment_direction*scale;

			double check_length = (point - current_point).norm();
			double check_difference = abs(check_length - segment_length);
			//assert(check_difference < 1e-3);
			write_log(loglevel) << "  -> CHECK: difference =  " << check_difference << ", length = " << check_length << " (" << (point - out_sampled.row(number_points_added - 1)).norm() << ")" << endl;
			write_log(loglevel) << "  adding point at " << number_points_added << ": " << point << endl << endl;

			out_sampled_length += check_length;
			out_sampled.row(number_points_added) = point;
			number_points_added++;

			current_point = point;
			continue;
		}

		//if is last point to add
		if (segment_end_index >= number_points - 1 && current_length < segment_length)
		{
			write_log(loglevel) << "in LAST segment: " << current_length << endl;
			write_log(loglevel) << "  remaining current_length: " << current_length << " (end_index: " << segment_end_index << ")" << endl;
			write_log(loglevel) << "  adding point at " << number_points_added << ": " << original_path.row(segment_end_index) << endl << endl;

			out_sampled_length += current_length;
			out_sampled.row(number_points_added) = original_path.row(segment_end_index);
			number_points_added++;

			break;
		}

		write_log(loglevel) << "in MULTIPLE segments: " << current_length << endl;
		write_log(loglevel) << "  segment (" << segment_start_index << " - " << segment_end_index << ")" << endl;

		segment_start_index = segment_end_index - 1;
		if (segment_start_index < 1)
			write_log(0) << "error" << endl;

		double a = (original_path.row(segment_start_index) - current_point).norm();
		double b = (original_path.row(segment_end_index) - original_path.row(segment_start_index)).norm();
		double c = (original_path.row(segment_end_index) - current_point).norm();
		write_log(loglevel) << "  length span: a = " << a << ", b = " << b << ", c = " << c << endl;

		double proportional_b = (segment_length - a) / (c - a);
		Eigen::RowVectorXd segment_direction = (original_path.row(segment_end_index) - original_path.row(segment_start_index)).normalized();
		Eigen::RowVectorXd point = original_path.row(segment_start_index) + segment_direction * (b * proportional_b);
		//auto point = current_point + segment_direction * (b * proportional_b);
		write_log(loglevel) << "  proportional_b: " << proportional_b << endl;

		double check_length = (point - current_point).norm();
		double check_difference = abs(check_length - segment_length);
		//assert(check_difference < 1e-3);
		write_log(loglevel) << "  -> CHECK: difference =  " << check_difference << ", length = " << check_length << " (" << (point - out_sampled.row(number_points_added - 1)).norm() << ")" << endl;
		write_log(loglevel) << "  adding point at " << number_points_added << ": " << point << endl << endl;

		out_sampled_length += check_length;
		out_sampled.row(number_points_added) = point;
		number_points_added++;

		current_point = point;
	}

	out_sampled.conservativeResize(number_points_added, Eigen::NoChange);
	write_log(loglevel) << "--> out_sampled: " << endl << out_sampled << endl;
	write_log(loglevel) << "--> target resampled_length: " << out_sampled_length << ", number points: " << number_points_added << endl << endl;


	//double geodesic_length = geodesics_candidates.geodesic_lengths(index);

	////double-check resampling
	//double off = geodesic_length - out_sampled_length;
	//assert(abs(off) > 5e-2);
	//write_log(2) << endl << endl << "ERROR: resampled geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(off) << endl << endl;

	//geodesic_length = resampled_length;
	//target_geodesic = resampled_geodesic;
	//data_model.target_curve = target_geodesic;

	//debug_show_curve(original_geodesic, data_model, true);
}


void curve::resample_matching_curve(const Eigen::MatrixXd& to_sample, const int center_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	Eigen::MatrixXd sampled_first_half;
	double sampled_length_first_half;
	sample_backward(to_sample, center_index, direction*(-1), sampled_first_half, sampled_length_first_half);
	write_log(6) << "sampled_first_half: " << endl << sampled_first_half << endl;

	Eigen::MatrixXd sampled_second_half;
	double sampled_length_second_half;
	sample_forward(to_sample, center_index, direction, sampled_second_half, sampled_length_second_half);
	write_log(6) << "sampled_second_half: " << endl << sampled_second_half << endl;

	out_sampled_length = sampled_length_first_half + sampled_length_second_half;

	out_sampled.resize(to_sample.rows(), 3);
	out_sampled <<
		sampled_first_half,
		Eigen::RowVector3d::Zero(),
		sampled_second_half;


	if (LOG_LEVEL == 6) //CHECK
	{
		write_log(0) << endl << endl;

		double to_sample_length = 0;
		double sampled_length = 0;
		for (int i = 1; i < to_sample.rows(); i++)
		{
			double length_original = (to_sample.row(i) - to_sample.row(i - 1)).norm();
			double length_resampled = (out_sampled.row(i) - out_sampled.row(i - 1)).norm();

			to_sample_length += length_original;
			sampled_length += length_resampled;

			write_log(0) << "length_original:  " << std::to_string(length_original) << endl;
			write_log(0) << "length_resampled: " << std::to_string(length_resampled) << endl << endl;

			if (abs(length_original - length_resampled) > 1e-3)
				write_log(0) << "   it's off at [" << i << "] by " << std::to_string(length_original - length_resampled) << endl << endl;
		}
		write_log(0) << "TOTAL length_original:  " << std::to_string(to_sample_length) << endl;
		write_log(0) << "TOTAL length_resampled: " << std::to_string(sampled_length) << endl << endl;
		write_log(0) << endl << endl;
	}

	write_log(6) << endl;
	write_log(4) << "wrapper sampled_length: " << out_sampled_length << endl;
	write_log(6) << "sampled_length_first_half: " << sampled_length_first_half << endl;
	write_log(6) << "sampled_length_second_half: " << sampled_length_second_half << endl;
	write_log(6) << "out_sampled: " << endl << out_sampled << endl;
	write_log(6) << endl;
}

void curve::sample_forward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	out_sampled.resize(to_sample.rows() - index - 1, 3);
	out_sampled_length = 0;

	int out_i = 0;
	for (int i = index + 1; i < to_sample.rows(); i++)
	{
		double segment_length = (to_sample.row(i - 1) - to_sample.row(i)).norm();
		out_sampled_length += segment_length;
		out_sampled.row(out_i) = direction * out_sampled_length;

		write_log(6) << "   geodesic: wrapper resampled[" << out_i << "] = " << out_sampled.row(out_i) << "   (segment_l=" << segment_length << ", summed_l=" << out_sampled_length << ")" << endl;

		out_i++;
	}
}

void curve::sample_backward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
{
	out_sampled.resize(index, 3);
	out_sampled_length = 0;

	for (int i = index - 1; i >= 0; i--)
	{
		double segment_length = (to_sample.row(i) - to_sample.row(i + 1)).norm();
		out_sampled_length += segment_length;
		out_sampled.row(i) = direction * out_sampled_length;

		write_log(6) << "   geodesic: wrapper resampled[" << i << "] = " << out_sampled.row(i) << "   (segment_l=" << segment_length << ", summed_l=" << out_sampled_length << ")" << endl;
	}
}

//TODO remove! it's legacy
double curve::compute_curvature_sum(const Eigen::MatrixXd& path)
{
	const int loglevel = 6;

	const int number_points = path.rows();
	double curvature = 0.0;

	write_log(loglevel) << "curvatures:  " << endl;
	for (int i = 0; i < number_points - 2; i++)
	{
		Eigen::RowVector3d e1 = (path.row(i + 1) - path.row(i)).normalized();
		Eigen::RowVector3d e2 = (path.row(i + 2) - path.row(i + 1)).normalized();

		double dot = e1.dot(e2);
		double cos_angle = clip(dot, -1, 1);

		if (cos_angle == 1 || cos_angle == -1)
			continue;

		double angle = acos(cos_angle);
		double hypothenuse = (path.row(i + 2) - path.row(i)).norm();
		double current_curvature = 2 * sin(angle) / hypothenuse;

		if (hypothenuse < 1e-3) //too short hypothenuse creates huge curvature. but with dense meshes this happens... 
			current_curvature = 0;

		curvature += current_curvature;
		write_log(loglevel) << current_curvature << ", ";
	}
	write_log(loglevel) << endl;

	return curvature;
}

double curve::compute_gauss_curvature_sum(Mesh& target, const Eigen::MatrixXd& path)
{
	const int loglevel = 6;

	const int number_points = path.rows();
	if (number_points < 1)
		return 0.0;

	Eigen::MatrixXi bary_indices;
	Eigen::MatrixXd bary_weigths;
	CoordinateConverter::barycentric_coords_from_points(path, target.V, target.F, bary_indices, bary_weigths);

	write_log(loglevel) << endl <<
		"bary_indices: " << endl << bary_indices << endl <<
		"bary_weigths: " << endl << bary_weigths << endl <<
		"gaussian_curvature_target: " << endl << target.surface_features().K << endl <<
		endl;

	//get all K for geodesic path points
	Eigen::VectorXd K_per_point(number_points);
	for (int j = 0; j < number_points; j++)
	{
		double k = 0;

		for (int bary_i = 0; bary_i < 3; bary_i++)
		{
			double weight = bary_weigths(j, bary_i);
			int vertex_index = bary_indices(j, bary_i);
			double vertex_K = abs(target.surface_features().K(vertex_index));

			k += weight * vertex_K;
		}
		K_per_point(j) = k;
	}
	write_log(loglevel) << "K_per_point:  " << endl << K_per_point << endl << endl;


	double curvature = 0.0;

	for (int j = 1; j < number_points - 1; j++)
	{
		//TODO somehow normalize for lengths
		//double prev_K = K_per_point(j - 1);
		//double next_K = K_per_point(j + 1);

		curvature += K_per_point(j);
	}

	// add K for first point
	curvature += K_per_point(0);

	// add K for last point
	curvature += K_per_point(number_points - 1);

	////apply this on sum
	//if (geodesic_lengths(i) > 1e-6)
	//{
	//	const double length_normalization = 1 / geodesic_lengths(i);
	//	curvature *= length_normalization;
	//}
	//else
	//{
	//	curvature = 0;
	//}

	write_log(loglevel) << " gauss curvature:  " << curvature << endl;

	return curvature;
}
