#include "NormalDeviationStopping.h"

#include <igl/PI.h>
#include <igl/barycentric_coordinates.h>

#include "MeshModel.h"
#include "GlobalSettings.h"
#include "Utils.h"
#include "Logger.h"

//int get_face_index(const Eigen::RowVectorXd& point, const Mesh& mesh);

NormalDeviationStopping::NormalDeviationStopping()
	: cluster(window_size, outlier_factor)
{ 
	absolute_threshold = GlobalSettings::crease_angle_threshold;
	bounds_min_absolute = GlobalSettings::flat_angle_threshold;
};


bool NormalDeviationStopping::should_stop()
{
	angle_sum += last_deviation;
	bool is_done_winding = angle_sum >= igl::PI * winding_factor;
	bool is_done_clustering = cluster.is_value_within(last_deviation) == false;
	
	//bool is_done(is_done_winding || is_done_clustering);
	//if (is_done)
	//	write_log(0) << "  STOPPING:: should stop" << linebreak << std::endl;

	return is_done_winding || is_done_clustering;
}

void NormalDeviationStopping::update(const WalkedData& data, const Mesh& mesh)
{
	if (data.faces.size() < 2)
		return;


	//update settings... should actually be done with pointers 
	cluster.window_size = window_size;
	cluster.outlier_factor = outlier_factor;

	cluster.absolute_threshold = absolute_threshold;
	cluster.bounds_min_absolute = bounds_min_absolute;
	cluster.use_upper_bound_only = use_upper_bound_only;
	cluster.use_adaptive_bounds_update = use_adaptive_bounds_update;


	//update last normal deviation
	int index = data.faces.size() - 1;
	int current_face_index = data.faces[index];
	int previous_face_index = data.faces[index - 1];

	last_deviation = normals_deviation(current_face_index, previous_face_index, mesh.normals_faces);

	//write_log(0) << "STOPPING:: previous_face_index: " << previous_face_index << " current_face_index: " << current_face_index << std::endl;
	//write_log(0) << "STOPPING:: last_deviation: " << last_deviation << std::endl;
}

void NormalDeviationStopping::reset()
{
	cluster.reset();
	last_deviation = 0.0;
	angle_sum = 0.0;
}

/*
void NormalDeviationStopping::update(const Eigen::MatrixXd& path, const Mesh& mesh, const int last_index)
{
	if (last_index < 2)
		return;


	//update settings... should actually be done with pointers 
	cluster.window_size = window_size;
	cluster.outlier_factor = outlier_factor;

	cluster.absolute_threshold = absolute_threshold;
	cluster.bounds_min_absolute = bounds_min_absolute;
	cluster.use_adaptive_bounds_update = use_adaptive_bounds_update;


	//update last normal deviation
	Eigen::RowVectorXd point = path.row(last_index);
	Eigen::RowVectorXd segment_direction = (point - path.row(last_index - 1)).normalized();

	int previous_face_index = get_face_index((point - segment_direction*0.001), mesh);
	int current_face_index  = get_face_index((point + segment_direction*0.001), mesh);

	write_log(0) << "STOPPING:: previous_face_index: " << previous_face_index << " current_face_index: " << current_face_index << std::endl;


	if (previous_face_index < 0 || current_face_index < 0)
	{
		last_deviation = cluster.absolute_threshold + 1;
		return;
	}

	last_deviation = normals_deviation(current_face_index, previous_face_index, mesh.normals_faces);
	write_log(0) << "STOPPING:: last_deviation: " << last_deviation << std::endl;
	//last_deviation = normals_deviation(fid_current, fid_previous, target.normals_faces);
}
*/

/*
int get_face_index(const Eigen::RowVectorXd& point, const Mesh& mesh)
{
	int face_index = -1;
	int closest_vertex = get_closest_vertex(mesh.V, point);

	for (int fi : mesh.adjacency_VF[closest_vertex])
	{
		Eigen::RowVectorXi face = mesh.F.row(fi);

		Eigen::RowVector3d barycentric_weights;
		igl::barycentric_coordinates(point, mesh.V.row(face(0)), mesh.V.row(face(1)), mesh.V.row(face(2)), barycentric_weights);

		//write_log(0) << "barycentric_weights: " << barycentric_weights << " at face: " << face << std::endl;

		if (barycentric_weights(0) >= 0.0 && barycentric_weights(0) <= 1.0 &&
			barycentric_weights(1) >= 0.0 && barycentric_weights(1) <= 1.0 &&
			barycentric_weights(2) >= 0.0 && barycentric_weights(2) <= 1.0)
		{
			face_index = fi;
			break;
		}
	}

	return face_index;
}
*/