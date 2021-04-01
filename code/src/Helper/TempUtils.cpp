#include "TempUtils.h"

void temp::get_geodesic(const int vertex, Mesh& target, GeodesicWalker& walker, const std::vector<int>& search_directions, Eigen::MatrixXd& out_geodesic, double& out_length)
{
	out_geodesic.resize(0,0);
	out_length = 0;

	Eigen::RowVector3d k_max = target.surface_features().principal_k_max_direction.row(vertex);
	Eigen::RowVector3d k_min = target.surface_features().principal_k_min_direction.row(vertex);
	Eigen::RowVector3d axis = k_max.cross(k_min).normalized();

	//double k_max_value = data_model.target.surface_features().principal_k_max(vertex);
	//double k_min_value = data_model.target.surface_features().principal_k_min(vertex);
	////write_log(0) << i << ": k_max = " << k_max << " k_min = " << k_min << " axis = " << axis << endl;
	//write_log(5) << "[" << vertex << "]: " << endl;
	//write_log(5) << "  k_max: " << " | " << k_max_value << " | , direction: " << k_max << endl;
	//write_log(5) << "  k_min: " << " | " << k_min_value << " | , direction: " << k_min << endl;
	//write_log(5) << "  axis:  " << axis << "  is flat? " << boolalpha << (is_equal(k_max_value, 0.0) && is_equal(k_min_value, 0.0)) << linebreak << endl;

	//create vectors for all search directions
	std::vector<Eigen::RowVector3d> directions;
	for (int angle_deg : search_directions)
	{
		Eigen::RowVector3d direction;

		if (angle_deg == 0)
		{
			direction = k_max;
		}
		else if (angle_deg == 90)
		{
			direction = k_min;
		}
		else
		{
			double angle = angle_deg / 180.0 * igl::PI;

			Eigen::Matrix3d R;
			R = Eigen::AngleAxisd(angle, axis);
			direction = k_max * R;
		}

		directions.push_back(direction);
	}

	//from all search directions, select longest geodesic
	for (const Eigen::RowVector3d& direction : directions)
	{
		Eigen::MatrixXd current_geodesic;
		double current_length;
		walker.get_geodesic_at(vertex, direction, current_geodesic, current_length);

		if (current_length > out_length)
		{
			out_geodesic = current_geodesic;
			out_length = current_length;
		}
	}

	//write_log(4) << "Geodesic at v" << vertex << ": length = " << length << ", rows = " << geodesic.rows() << ", is flat ? " << boolalpha << curve::is_flat(geodesic) << endl << endl;
}