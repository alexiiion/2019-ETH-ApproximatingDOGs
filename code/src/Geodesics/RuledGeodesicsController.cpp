#include "RuledGeodesicsController.h"

#include <experimental/filesystem>

#include <Eigen/Geometry> 
#include <igl/Timer.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
//#include <igl/vertex_triangle_adjacency.h>
#include <igl/barycentric_coordinates.h>

#include "PatchModel.h"
#include "MeshController.h"

#include "GeodesicWalker.h"
#include "LengthGeodesicsStopping.h"
#include "RuledDevelopableSurface.h"
#include "MeshLabelling.h"
#include "GraphCut.h"

#include "CurveHelper.h"
#include "NormalsHelper.h"
#include "CurveInterpolation.h"
#include "CoordinateConverter.h"
#include "Utils.h"
#include "Logger.h"

#include "GlobalSettings.h"

using namespace std;

//const Eigen::RowVector3d sample_direction = Eigen::RowVector3d::UnitX();

std::vector<int> get_faces(const std::vector<int>& component_vertices, const std::vector<std::vector<int>>& adjacency_VF);
void resample_constraint_curves(const Eigen::MatrixXd& geodesic, const double length, Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve);

//const double area_scale = 2.0;
//const double angle_threshold = 0.02;
void dfs_flat_on_faces(const std::vector<std::vector<int>>& connectivity_list, const Eigen::MatrixXd& N, const int& start_index, std::vector<int>& out_component);

bool RuledGeodesicsController::initialize()
{
	if (is_initialzed)
		return true;

	if (data_model.geodesics_candidates.is_empty())
		build_geodesics();


	//data_model.ruled_vertex_indices.clear();
	//data_model.ruled_vertex_indices.push_back(2838);

	if (data_model.ruled_vertex_indices.size() < 1 && data_model.ruled_developables.size() < 1)
		build_random_developables(data_model.number_random_points);
	else if (data_model.ruled_vertex_indices.size() >= 1 && data_model.ruled_developables.size() < 1) //loading indices from file, generating surfaces for those
		data_model.ruled_developables = build_developables_at(data_model.ruled_vertex_indices);
	else
		write_log(3) << endl << "INFO at initializing Geodesics Controller:: RANDOM geodesics -- already initialized, all good" << endl;

	data_model.number_random_points = data_model.ruled_developables.size();


	if (data_model.selected_ruled_vertex_indices.size() < 1 && data_model.selected_ruled_view_indices.size() < 1)
		select_ruled_developables();
	else if (data_model.selected_ruled_vertex_indices.size() >= 1 && data_model.selected_ruled_view_indices.size() < 1)
		selector = new MeritSelector(data_model.selected_ruled_vertex_indices);
		//log_list(4, data_model.selected_ruled_vertex_indices, "loading SELECTED geodesics from file: ", false, false);
	else
		write_log(3) << endl << "INFO at initializing Geodesics Controller:: SELECTED geodesics -- already initialized, all good" << endl;

	if (data_model.selected_ruled_vertex_indices <= data_model.ruled_vertex_indices)
	{
		std::vector<Mesh> selected_ruled_developables = meshhelper::get_selected_ruled_developable_meshes(data_model.selected_ruled_vertex_indices, data_model.ruled_vertex_indices, data_model.ruled_developables);
		data_model.initialization_coverage->compute(selected_ruled_developables);
	}

	is_initialzed = true;
	return is_initialzed;
}

void RuledGeodesicsController::build_geodesics()
{
	if (!data_model.geodesics_candidates.is_empty())
		data_model.geodesics_candidates.clear();

	igl::Timer timer;
	double init_time = timer.getElapsedTime();

	write_log(3) << endl << "build all walking geodesics..." << endl;

	const int number_vertices = data_model.target.V.rows();
	GeodesicCandidates candidates(number_vertices);

	Eigen::MatrixXi unused_adjacent_vertices;
	Eigen::VectorXi unused_incident_faces;

	//Eigen::MatrixXd bary_weigths;

	GeodesicWalker walker(data_model.target, data_model.geodesics_stopping);
	//GeodesicWalker walker(data_model.target, &data_model.max_frame_error, &data_model.use_weighted_frame_error);

	for (int i = 0; i < number_vertices; i++)
	{
		LOG_LEVEL = 4;
		if (i == 916)
		{
			write_log(0) << "916" << endl;
			//LOG_LEVEL = 5;
		}

		Eigen::MatrixXd geodesic;
		double length = 0;
		int deg = 0;

		Eigen::RowVector3d k_max = data_model.target.surface_features().principal_k_max_direction.row(i);
		double k_max_value = data_model.target.surface_features().principal_k_max(i);
		Eigen::RowVector3d k_min = data_model.target.surface_features().principal_k_min_direction.row(i);
		double k_min_value = data_model.target.surface_features().principal_k_min(i);

		Eigen::RowVector3d axis = k_max.cross(k_min).normalized();
		
		//write_log(0) << i << ": k_max = " << k_max << " k_min = " << k_min << " axis = " << axis << endl;
		write_log(5) << "[" << i << "]: " << endl;
		write_log(5) << "  k_max: " << " | " << k_max_value << " | , direction: " << k_max << endl;
		write_log(5) << "  k_min: " << " | " << k_min_value << " | , direction: " << k_min << endl;
		write_log(5) << "  axis:  " << axis << "  is flat? " << boolalpha << (is_equal(k_max_value, 0.0) && is_equal(k_min_value, 0.0)) << linebreak << endl;

		//TODO build all geodesics for std::vector<int> geodesics_directions; store the longest at vertex
		for (int d = 0; d < data_model.geodesics_directions.size(); d++)
		{
			int angle_deg = data_model.geodesics_directions[d];
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
			//write_log(5) << "axis: " << axis <<  ", angle_deg: " << angle_deg << ", rotated direction: " << direction << ", R: " << endl << R << endl;
			write_log(5) << "axis: " << axis <<  ", angle_deg: " << angle_deg << ", rotated direction: " << direction << endl;
			
			Eigen::MatrixXd current_geodesic;
			double current_length;

			////float max_length = 10;
			////GeodesicWalker walker(data_model.target, &max_length);
			walker.get_geodesic_at(i, direction, current_geodesic, current_length);
			write_log(5) << "v" << i << " current_geodesic.rows(): " << current_geodesic.rows() << ", current_length: " << current_length << endl << endl;

			if (current_length > length)
			{
				geodesic = current_geodesic;
				length = current_length;
				deg = angle_deg;
			}
		}

		//if (length < 0.00001)
		//	continue;

		//if (curve::is_flat(geodesic))
		//	candidates.flat_geodesics.push_back(i);


		write_log(5) << "v" << i << " selected geodesic.rows(): " << geodesic.rows() << ", length: " << length << " at " << deg << "°" << endl;

		candidates.paths.push_back(geodesic);
		candidates.lengths(i) = length;
		candidates.number_geodesics++;
	}

	data_model.geodesics_candidates = candidates;
	//data_model.do_rebuild_geodesics = false;

	double t = timer.getElapsedTime();
	//log_list(0, candidates.flat_geodesics, "flat_geodesics: ", false);
	write_log(3) << "...done building geodesics. found " << candidates.paths.size() << " geodesics -- ELAPSED TIME: " << t - init_time << endl << endl;
}

void RuledGeodesicsController::build_random_developables(int n)
{
	write_log(3) << endl << "create ~" << n << " random ruled developable surfaces..." << endl;

	igl::Timer timer;
	double init_time = timer.getElapsedTime();


	//clear 
	for (int i = 0; i < data_model.ruled_developables.size(); i++)
		delete data_model.ruled_developables[i];
	data_model.ruled_developables.clear();

	for (int i = 0; i < data_model.selected_ruled_view_indices.size(); i++)
		data_model.viewer.data_list[data_model.selected_ruled_view_indices[i]].clear();
	data_model.selected_ruled_view_indices.clear();

	if (data_model.ruled_view_index >= 0)
	{
		data_model.viewer.data_list[data_model.ruled_view_index].clear();
		data_model.ruled_view_index = -1;
	}

	data_model.ruled_vertex_indices.clear();
	patch_vertex_distances.clear();
	flat_indices.clear();
	flat_area_indices.clear();



	//get n random vertices (later: farthest point sampling)

	const int n_vertices = data_model.target.V.rows();

	bool create_all_surfaces = false;
	if (n >= n_vertices) //TODO go trough them and check length & validity of developables, remove invalid vertices
	{
		//data_model.ruled_vertex_indices = std::vector<int>(n_vertices);
		//std::iota(data_model.ruled_vertex_indices.begin(), data_model.ruled_vertex_indices.end(), 0);
		data_model.ruled_vertex_indices.clear();

		for (int i = 0; i < n_vertices; i++)
		{

			Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[i];
			double length = data_model.geodesics_candidates.lengths[i];

			if (length < 2) //TODO extract parameter!
				continue;

			/*
			RuledDevelopableSurface* surface;
			bool is_valid_surface = try_get_ruled_developable(i, surface);
			if (!is_valid_surface)
				continue;
			*/

			RuledDevelopableSurface* surface;
			std::vector<int> used_vertices;

			bool is_valid_surface = try_get_surface_at(i, surface, used_vertices);
			if (!is_valid_surface)
				continue;

			data_model.ruled_developables.push_back(surface);
			data_model.ruled_vertex_indices.push_back(i);
		}

		data_model.number_random_points = data_model.ruled_vertex_indices.size();
		return;
	}


	//get n random vertices (later: farthest point sampling)

	//keep a list of invalid_indices (flat ones, or ones that previously gave invalid surface)
	vector<int> invalid_indices;
	
	int added_indices = 0;
	int iterations = 0;
	int max_iterations = n * 3.0;
	std::srand(unsigned(std::time(0)));

	while (added_indices < n && iterations < max_iterations)
	{
		iterations++;

		int random_index = std::rand() % n_vertices;

		if (is_contained(random_index, data_model.ruled_vertex_indices))
			continue;
		if (is_contained(random_index, invalid_indices))
			continue;


		RuledDevelopableSurface* surface;
		std::vector<int> used_vertices;

		bool is_valid_surface = try_get_surface_at(random_index, surface, used_vertices);
		write_log(5) << "[" << random_index << "]: valid surface? " << boolalpha << is_valid_surface << ", used_vertices: " << list_to_string(used_vertices) << endl;

		if (!is_valid_surface)
		{
			invalid_indices.push_back(random_index);
			continue;
		}

		invalid_indices.insert(invalid_indices.end(), used_vertices.begin(), used_vertices.end());
		sort(invalid_indices.begin(), invalid_indices.end());


		/*
		//Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[random_index];
		//double length = data_model.geodesics_candidates.lengths[random_index];

		//if (length < 2) //TODO extract parameter!
		//	continue;


		//RuledDevelopableSurface* surface;
		//bool is_valid_surface = try_get_surface(random_index, surface);

		//TODO check if geodesic is flat. if flat, run dfs and label all vertices in data_model.ruled_vertex_indices, create only one simple ruled surface of correct size
		//bool is_flat = is_contained(random_index, data_model.geodesics_candidates.flat_geodesics);
		bool is_flat = curve::is_flat(data_model.geodesics_candidates.paths[random_index]);

		RuledDevelopableSurface* surface;
		bool is_valid_surface = false;

		if (is_flat)
		{
			std::vector<int> flat_vertices;
			is_valid_surface = try_get_flat_surface(random_index, surface, flat_vertices);
			if(is_valid_surface)
				invalid_indices.insert(invalid_indices.end(), flat_vertices.begin(), flat_vertices.end());
		}
		else
		{
			is_valid_surface = try_get_ruled_developable(random_index, surface);
		}

		if (!is_valid_surface)
		{
			invalid_indices.push_back(random_index);
			write_log(5) << "invalid surface at " << random_index << endl;
			continue; 
		}
		*/

		data_model.ruled_developables.push_back(surface);
		data_model.ruled_vertex_indices.push_back(random_index);
		added_indices++;

		//write_log(0) << "     ruled_vertex_indices: " << list_to_string(data_model.ruled_vertex_indices) << endl;
		//write_log(0) << "     invalid_indices: " << list_to_string(invalid_indices) << endl << endl;
	}

	data_model.number_random_points = data_model.ruled_developables.size();

	////TODO move to view!
	//Mesh concatenated_developables = meshhelper::concatenate_ruled_developables(data_model.ruled_developables);
	//meshhelper::add_mesh(data_model.viewer, concatenated_developables.V, concatenated_developables.F, false, false);
	//data_model.ruled_view_index = data_model.viewer.selected_data_index;


	write_log(4) << "  random indices for valid ruled developables: " << linebreak << list_to_string(data_model.ruled_vertex_indices, ", ") << endl;

	double t = timer.getElapsedTime();
	write_log(3) << "...done creating ruled developables -- ELAPSED TIME: " << t - init_time << endl << endl;
}

std::vector<RuledDevelopableSurface*> RuledGeodesicsController::build_developables_at(std::vector<int>& ruled_indices)
{
	log_list(3, ruled_indices, "loading geodesics from file: ", false, false);

	std::vector<RuledDevelopableSurface*> developables;
	//std::vector<int> valid_indices;
	std::vector<int> to_delete;

	for (int i = 0; i < ruled_indices.size(); i++)
	{
		int index = ruled_indices[i];

		//RuledDevelopableSurface* surface;
		//bool is_valid_surface = try_get_ruled_developable(index, surface);
	
		std::vector<int> used_vertices;

		RuledDevelopableSurface* surface;
		bool is_valid_surface = try_get_surface_at(index, surface, used_vertices);
		/*
		if (!is_valid_surface)
			to_delete.push_back(i);
		else
			developables.push_back(surface);
		*/
		
		if (!is_valid_surface)
		{
			to_delete.push_back(i);
			write_log(4) << "    ruled surface at v" << index << " is not valid!" << endl;
		}
		else
		{
			//if(index != ruled_indices[i])
			//	to_delete.push_back(ruled_indices[i]);

			ruled_indices[i] = index;
			developables.push_back(surface);
		}
	}

	for (int i = 0; i < to_delete.size(); i++)
		ruled_indices.erase(ruled_indices.begin() + to_delete[i] - i - 1);

	write_log(4) << "  created " << developables.size() << " developables." << endl;
	return developables;
}

void RuledGeodesicsController::select_ruled_developables()
{
	write_log(3) << endl << "select ruled developables to keep..." << endl;

	if (data_model.ruled_vertex_indices.size() < 1)
		return;

	have_checked_holes = false;

	if (data_model.ruled_vertex_indices.size() == 1)
	{
		int index = data_model.ruled_vertex_indices[0];
		data_model.selected_ruled_vertex_indices.push_back(index);
		data_model.vertex_ruled_assignment = vector<int> (data_model.target.V.rows(), index);

		//create_ruled_constraints_selector();
		selector = new MeritSelector(data_model.selected_ruled_vertex_indices);
		return;
	}

	igl::Timer timer;
	double init_time = timer.getElapsedTime();


	if (patch_vertex_distances.size() < 1)
	{
		write_log(3) << "  compute distance from each vertex to each patch" << endl;

		std::vector<Eigen::MatrixXd> patch_vertex_normals; //the normals on the patches for the projected points
		//std::vector<Eigen::VectorXd> patch_vertex_normal_deviations; //the normals on the patches for the projected points

		/*
		std::vector<Eigen::VectorXd> patch_ray_deviations;
		std::vector<Eigen::MatrixXd> patch_closest_points; //for DEBUG
		std::vector<Eigen::MatrixXd> patch_rays; //for DEBUG
		*/

		Eigen::VectorXi face_indices;
		Eigen::MatrixXd closest_points;

		//for (int i = 0; i < data_model.number_random_points; i++) 
		for (int i = 0; i < data_model.ruled_developables.size(); i++)
		{
			Mesh ruled_mesh = data_model.ruled_developables[i]->developable;
			if (ruled_mesh.F.rows() < 1)
			{
				write_log(1) << "ERROR:: in select_ruled_developables() where ruled_mesh.F.rows() < 1" << endl;
				continue;
			}

			Eigen::VectorXd distances_squared;
			igl::point_mesh_squared_distance(data_model.target.V, ruled_mesh.V, ruled_mesh.F, distances_squared, face_indices, closest_points);

			distances_squared = distances_squared.cwiseSqrt();
			patch_vertex_distances.push_back(distances_squared);
		}
	}

	/*
	if (patch_vertex_distances.size() < 1)
	{
		write_log(3) << "  compute distance from each vertex to each patch" << endl;

		std::vector<Eigen::MatrixXd> patch_vertex_normals; //the normals on the patches for the projected points
		//std::vector<Eigen::VectorXd> patch_vertex_normal_deviations; //the normals on the patches for the projected points

		//std::vector<Eigen::VectorXd> patch_ray_deviations;
		//std::vector<Eigen::MatrixXd> patch_closest_points; //for DEBUG
		//std::vector<Eigen::MatrixXd> patch_rays; //for DEBUG
	
		Eigen::VectorXi face_indices;
		Eigen::MatrixXd closest_points;

		//for (int i = 0; i < data_model.number_random_points; i++) 
		for (int i = 0; i < data_model.ruled_developables.size(); i++)
		{
			Mesh ruled_mesh = data_model.ruled_developables[i]->developable;
			if (ruled_mesh.F.rows() < 1) 
			{
				write_log(1) << "ERROR:: in select_ruled_developables() where ruled_mesh.F.rows() < 1" << endl;
				continue;
			}

			Eigen::VectorXd distances_squared;
			igl::point_mesh_squared_distance(data_model.target.V, ruled_mesh.V, ruled_mesh.F, distances_squared, face_indices, closest_points);

			distances_squared = distances_squared.cwiseSqrt();
			patch_vertex_distances.push_back(distances_squared);



			//Eigen::MatrixXd VN;
			//igl::per_vertex_normals(ruled_mesh.V, ruled_mesh.F, VN);
			////std::vector<std::vector<int>> VFi_unused;
			////std::vector<std::vector<int>> adjacency_VF;
			////igl::vertex_triangle_adjacency(ruled_mesh.V, ruled_mesh.F, adjacency_VF, VFi_unused);

			//Eigen::MatrixXd surface_normals(closest_points.rows(), closest_points.cols());
			//Eigen::VectorXd n_deviations(closest_points.rows());
			//for (int pi = 0; pi < closest_points.rows(); pi++)
			//{
			//	//Eigen::RowVectorXd normal = get_surface_normal(closest_points.row(pi), ruled_mesh.V, ruled_mesh.F, VN, adjacency_VF);

			//	Eigen::RowVector3i face = ruled_mesh.F.row(face_indices(pi));
			//	Eigen::RowVector3d barycentric_weights;
			//	igl::barycentric_coordinates(closest_points.row(pi).transpose(), ruled_mesh.V.row(face(0)), ruled_mesh.V.row(face(1)), ruled_mesh.V.row(face(2)), barycentric_weights);

			//	Eigen::Vector3d normal(0, 0, 0);
			//	for (int bary_i = 0; bary_i < 3; bary_i++)
			//	{
			//		double weight = barycentric_weights(bary_i);
			//		int vertex_index = face(bary_i);
			//		normal += weight * VN.row(vertex_index);
			//	}

			//	surface_normals.row(pi) = normal;


			//	
			//	n_deviations(pi) = angle(normal, data_model.target.normals_vertices.row(pi));

			//}
			//patch_vertex_normals.push_back(surface_normals);
			

			//Eigen::MatrixXd surface_normals = utils::compute_surface_normals_vertices(closest_points, ruled_mesh);
			Eigen::MatrixXd surface_normals = utils::compute_surface_normals_faces(closest_points, ruled_mesh);

			int v = data_model.ruled_vertex_indices[i];
			double a = angle(data_model.target.normals_vertices.row(v), surface_normals.row(v));
			write_log(5) << " normal deviation at ruled v" << v << ": " << a << endl;

			if (a >= igl::PI*0.75)
			{
				write_log(5) << "   invert ruled normals at v" << v << endl;
				ruled_mesh.normals_vertices *= -1;
				surface_normals *= -1;
			}

			patch_vertex_normals.push_back(surface_normals);


			Eigen::VectorXd n_deviations(closest_points.rows());

			//get angles between vertex normals and closest point on each patch, use a penalty in graph cut
			Eigen::MatrixXd r(closest_points.rows(), 3);
			Eigen::VectorXd ray_deviations(closest_points.rows());
			for (int pi = 0; pi < closest_points.rows(); pi++)
			{
				n_deviations(pi) = angle(surface_normals.row(pi), data_model.target.normals_vertices.row(pi));



				//Eigen::RowVector3d ray = (closest_points.row(pi) - data_model.target.V.row(pi));
				//double normal_deviation = angle(ray, data_model.target.normals_vertices.row(pi));
				//ray_deviations(pi) = normal_deviation;

				//r.row(pi) = ray;
			}
			patch_vertex_normal_deviations.push_back(n_deviations);

			//patch_closest_points.push_back(closest_points);
			//patch_ray_deviations.push_back(ray_deviations);
			//patch_rays.push_back(r);

		}
	}

	
	string folder = "./__data/";
	bool success = false;
	if (!std::experimental::filesystem::exists(folder))
		success = std::experimental::filesystem::create_directories(folder);

	//int number_torsal = data_model.ruled_developables.size();
	//int number_vertices = data_model.target.V.rows();

	//ofstream file_distances(folder+"1-patch_vertex_distances.txt");
	//	Eigen::MatrixXd distance_matrix(number_torsal, number_vertices);
	//	for (int i = 0; i < patch_vertex_distances.size(); ++i)
	//	{
	//		Eigen::VectorXd distances = patch_vertex_distances[i];
	//		for (int j = 0; j < distances.rows(); j++)
	//			distance_matrix(i, j) = distances(j);
	//	}
	//	file_distances << distance_matrix << endl;
	//file_distances.close();



	//ofstream file_deviations(folder + "1-patch_vertex_normal_deviations.txt");
	//	Eigen::MatrixXd deviation_matrix(number_torsal, number_vertices);
	//	for (int i = 0; i < patch_vertex_normal_deviations.size(); ++i)
	//	{
	//		Eigen::VectorXd deviations = patch_vertex_normal_deviations[i];
	//		for (int j = 0; j < deviations.rows(); j++)
	//			deviation_matrix(i, j) = deviations(j);
	//	}
	//	file_deviations << deviation_matrix << endl;
	//file_deviations.close();



	//ofstream file_raydeviations(folder + "1-patch_vertex_ray_deviations.txt");
	//	Eigen::MatrixXd raydeviation_matrix(number_torsal, number_vertices);
	//	for (int i = 0; i < patch_ray_deviations.size(); ++i)
	//	{
	//		Eigen::VectorXd deviations = patch_ray_deviations[i];
	//		for (int j = 0; j < deviations.rows(); j++)
	//			raydeviation_matrix(i, j) = deviations(j);
	//	}
	//	file_raydeviations << raydeviation_matrix << endl;
	//file_raydeviations.close();

	//
	//
	////ofstream file_rays("1-patch_vertex_rays.txt");
	//for (int i = 0; i < number_torsal; i++)
	//{
	//	//Eigen::MatrixXd points = patch_closest_points[i];
	//	//Eigen::MatrixXd rays = points - data_model.target.V;
	//	Eigen::MatrixXd rays = patch_rays[i];

	//	ofstream file_rays(folder + "1" + to_string(i) + "-patch_vertex_rays.txt");
	//	file_rays << rays << endl;
	//	file_rays.close();
	//}
	////file_rays.close();

	//
	//
	////ofstream file_rays("1-patch_vertex_rays.txt");
	//for (int i = 0; i < number_torsal; i++)
	//{
	//	//Eigen::MatrixXd points = patch_closest_points[i];
	//	//Eigen::MatrixXd rays = points - data_model.target.V;
	//	Eigen::MatrixXd normals = patch_vertex_normals[i];

	//	ofstream file_normals(folder + "1" + to_string(i) + "-patch_normals.txt");
	//	file_normals << normals << endl;
	//	file_normals.close();
	//}
	////file_rays.close();

	//
	//
	////ofstream file_rays("1-patch_vertex_rays.txt");
	//for (int i = 0; i < number_torsal; i++)
	//{
	//	Eigen::MatrixXd points = patch_closest_points[i];

	//	ofstream file_points(folder + "1" + to_string(i) + "-patch_points.txt");
	//	file_points << points << endl;
	//	file_points.close();
	//}
	////file_rays.close();

	//
	//
	//ofstream file_normals(folder + "1-vertex_normals.txt");
	//	file_normals << data_model.target.normals_vertices << endl;
	//file_normals.close();
	*/

	write_log(4) << "  optimize label assignment" << endl << endl;
	std::vector<double> resulting_distances;
	vector<int> assignment = graphcut::compute_graph_cut_labeling(data_model.target, patch_vertex_distances, data_model.label_selection_smoothness, resulting_distances);
	//vector<int> assignment = graphcut::compute_graph_cut_labeling(data_model.target, patch_vertex_distances, patch_vertex_normal_deviations, data_model.label_selection_smoothness, resulting_distances);
	//vector<int> assignment = graphcut::compute_graph_cut_labeling(data_model.target, patch_vertex_distances, patch_vertex_normals, data_model.label_selection_smoothness, resulting_distances);
	//vector<int> assignment = graphcut::compute_graph_cut_labeling(data_model.target, patch_vertex_distances, data_model.label_selection_smoothness);
	//assignment: size is #V, values are patch index {0, number_random_points} for each vertex

	//ofstream file(folder + "selected_patches_distances.txt");
	//for (double value : resulting_distances)
	//	file << value << endl;
	//file.close();


	//re-map assignment from patch indices to geodesic (vertex) indices
	//extract unique geodesic (vertex) indices
	data_model.vertex_ruled_assignment.clear();
	data_model.vertex_ruled_assignment.reserve(data_model.target.V.rows());

	data_model.selected_ruled_vertex_indices.clear();

	if (data_model.selected_ruled_view_indices.size() > 0)
	{
		for (int i = 0; i < data_model.selected_ruled_view_indices.size(); i++)
			data_model.viewer.data_list[data_model.selected_ruled_view_indices[i]].clear();

		data_model.selected_ruled_view_indices.clear();
	}

	for (int i = 0; i < assignment.size(); i++)
	{
		//re-map patch index to vertex index
		int patch_index = assignment[i];
		int geodesic_vertex_index = data_model.ruled_vertex_indices[patch_index];
		data_model.vertex_ruled_assignment.push_back(geodesic_vertex_index);

		//store unique selected indices
		if (is_contained(geodesic_vertex_index, data_model.selected_ruled_vertex_indices))
			continue;

		data_model.selected_ruled_vertex_indices.push_back(geodesic_vertex_index);
	}


	//add previously saved flat parts
	if (data_model.use_geodesics_flat_detetction)
	{
		for (int i = 0; i < flat_indices.size(); i++)
		{
			int flat_index = flat_indices[i];
			if (is_contained(flat_index, data_model.selected_ruled_vertex_indices))
				continue;

			if (data_model.geodesics_candidates.paths[flat_index].rows() < 5)
				continue;

			data_model.selected_ruled_vertex_indices.push_back(flat_index);

			for (int flat_vi : flat_area_indices[i])
				data_model.vertex_ruled_assignment[flat_vi] = flat_index;

			int index = index_of(flat_index, data_model.ruled_vertex_indices);
		}
	}



	double t = timer.getElapsedTime();
	write_log(3) << "...found " << data_model.selected_ruled_vertex_indices.size() << " ruled developables to keep -- ELAPSED TIME: " << t - init_time << endl << endl;
	write_log(4) << "...selected indices: " << list_to_string(data_model.selected_ruled_vertex_indices, ", ") << endl;


	create_ruled_constraints_selector();

	std::vector<Mesh> selected_ruled_developables = meshhelper::get_selected_ruled_developable_meshes(data_model.selected_ruled_vertex_indices, data_model.ruled_vertex_indices, data_model.ruled_developables);
	data_model.initialization_coverage->compute(selected_ruled_developables);



	//only debug!
	vector<int> selected_patch_indices = assignment;
	std::sort(selected_patch_indices.begin(), selected_patch_indices.end());
	selected_patch_indices.erase(unique(selected_patch_indices.begin(), selected_patch_indices.end()), selected_patch_indices.end());

	log_list(5, data_model.ruled_vertex_indices, "  random vertex indices: ", true, true);
	log_list(5, selected_patch_indices, "  selected_patch_indices: ", true, true);
	log_list(5, data_model.selected_ruled_vertex_indices, "  selected vertex indices: ", true, true);
}

/* NEW try_get_surface_at() */
bool RuledGeodesicsController::try_get_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	if (vertex < 0 || vertex >= data_model.geodesics_candidates.paths.size())
		return false;

	
	Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[vertex];
	if (geodesic.rows() < 4)
		return false;

	bool is_surface_valid = false;
	if (data_model.use_geodesics_flat_detetction)
	{
		bool is_flat = curve::is_flat(geodesic);
		if (is_flat)
			return try_get_flat_surface_at(vertex, out_surface, assigned_indices);

		//if (is_flat)
		//	is_surface_valid = try_get_flat_surface_at(vertex, out_surface, assigned_indices);
		//if (is_surface_valid)
		//	return true;
		//else
		//	write_log(0) << "is flat but invalid surface" << endl;
	}

	return try_get_curved_surface_at(vertex, out_surface, assigned_indices);
}

/* OLD try_get_surface_at() */
/*
bool RuledGeodesicsController::try_get_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	bool is_flat = is_contained(vertex, data_model.geodesics_candidates.flat_geodesics);
	bool is_valid_surface = false;

	if (is_flat)
	{
		std::vector<int> flat_vertices;
		is_valid_surface = try_get_flat_surface(vertex, out_surface, flat_vertices);
		//if (is_valid_surface)
		//	invalid_indices.insert(invalid_indices.end(), flat_vertices.begin(), flat_vertices.end());
	}
	else
	{
		is_valid_surface = try_get_ruled_developable(vertex, out_surface);
	}

	return is_valid_surface;
}
*/

/* NEW try_get_curved_surface_at() */
bool RuledGeodesicsController::try_get_curved_surface_at(const int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	Eigen::MatrixXd* geodesic = &data_model.geodesics_candidates.paths[vertex];
	if (geodesic->rows() < 6) //otherwise can't create spline of degree 5
		return false;

	out_surface = new RuledDevelopableSurface(data_model.target, *geodesic, vertex, RulingsType::Analytic);
	out_surface->create(data_model.ruled_width);

	assigned_indices.clear();
	assigned_indices.push_back(vertex);

	if (out_surface->developable.F.rows() < 3) //TODO extract parameter!
		return false;

	return true;
}

/* OLD try_get_curved_surface_at() */
/*
bool RuledGeodesicsController::try_get_curved_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	if (vertex < 0 || vertex >= data_model.geodesics_candidates.lengths.rows())
		return false;

	//write_log(4) << endl << "  create ruled developable at vertex "  << vertex << "..." << endl;
	//igl::Timer timer;
	//double init_time = timer.getElapsedTime();



	Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[vertex];
	if (geodesic.rows() < 6) //otherwise can't create spline of degree 5
		return false;

	out_surface = new RuledDevelopableSurface(data_model.target, data_model.geodesics_candidates.paths[vertex], vertex, RulingsType::Analytic);
	out_surface->create(data_model.ruled_width);

	if (out_surface->developable.F.rows() < 10) //TODO extract parameter!
		return false;

	//write_log(4) << endl << "  ...ruled developable successfully created at vertex " << vertex << " -- ELAPSED TIME : " << timer.getElapsedTime() - init_time << endl;

	return true;
}
*/

/* NEW try_get_flat_surface_at() */
bool RuledGeodesicsController::try_get_flat_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	const int loglevel = 5;

	assigned_indices.clear();
	assigned_indices.push_back(vertex);

	//check one ring neighborhood first, to see if starts at crease
	auto vertex_normal = data_model.target.normals_vertices.row(vertex);
	write_log(loglevel) << linebreak << "try get FLAT area at " << vertex << ", normal = " << vertex_normal << endl;

	for (int fi : data_model.target.adjacency_VF[vertex])
	{
		auto n = data_model.target.normals_faces.row(fi);
		double a = angle(n, vertex_normal);
		//double k = 2 * std::sin(a) / (n - vertex_normal).norm();
		write_log(loglevel) << a << ",";

		if (a > GlobalSettings::flat_angle_threshold)
		{
			write_log(loglevel) << "try get flat surface: not valid b/c starts at crease" << endl;
			return false;
		}
		//if (k > GlobalSettings::flat_curvature_threshold)
		//{
		//	write_log(loglevel) << "try get flat surface: not valid b/c starts at crease" << endl;
		//	return false;
		//}
	}
	write_log(loglevel) << "does not start at crease, proceed with labeling" << endl;

	vector<int> flat_vertices = label::dfs_flat_area(vertex, data_model.target.adjacency_VV, data_model.target.normals_vertices);
	if (flat_vertices.size() < 1)
	{
		write_log(loglevel) << "try get flat surface: not valid b/c no labelling found" << endl;
		return false;
	}
	write_log(loglevel) << "labeled " << flat_vertices.size() << "/" << data_model.target.V.rows() << " as 'flat'" << endl;

	
	/*
	int original_vertex = vertex; //for debug only
	vector<int> flat_vertices;
	bool found_flat = label::try_label_flat_area(original_vertex, data_model.target, flat_vertices, vertex);

	
	
	int number_labels = data_model.vertex_labels.maxCoeff();
	number_labels++;
	for (int i : flat_vertices)
		data_model.vertex_labels(i) = number_labels;

	//data_model.vertex_ruled_assignment = flat_vertices;

	//write_log(0) << "data_model.vertex_labels: " << linebreak << data_model.vertex_labels << endl;
	//write_log(0) << "flat indices: " << list_to_string(flat_vertices) << endl;


	if (vertex == 33)
		write_log(0) << "33" << endl;
	*/


	//using the center vertex of the flat area for better coverage
	int center_vertex = get_central_vertex(flat_vertices, data_model.target.V);
	int original_vertex = vertex; //for debug only
	vertex = center_vertex;

	Eigen::MatrixXd* geodesic = &data_model.geodesics_candidates.paths[vertex];
	if (geodesic->rows() < 3)
	{
		write_log(loglevel) << "try get flat surface: not valid b/c geodesic to short" << endl;
		return false;
	}

	if (!curve::is_flat(*geodesic))
	{
		write_log(4) << linebreak << "after centering geodesic isn't flat anymore!" << linebreak << endl;

		double max_length = 0.0;
		int new_center_vertex = center_vertex;
		for (int flat_i : flat_vertices)
		{
			Eigen::MatrixXd g = data_model.geodesics_candidates.paths[flat_i];
			if (!curve::is_flat(g))
				continue;

			double l = curve::compute_length(g);
			if (l > max_length)
			{
				max_length = l;
				new_center_vertex = flat_i;
			}
		}

		write_log(4) << "found new center vertex at " << new_center_vertex << endl;


		vertex = new_center_vertex;
		geodesic = &data_model.geodesics_candidates.paths[vertex];
	}

	/*
	//using the center leads to bad geodesics sometimes, therefore find the longest that is flat (hoping for the best...)
	int original_vertex = vertex; //for debug
	Eigen::MatrixXd* geodesic;

	double max_length = 0.0;
	int new_vertex = vertex;
	for (int flat_i : flat_vertices)
	{
		Eigen::MatrixXd g = data_model.geodesics_candidates.paths[flat_i];
		if (!curve::is_flat(g))
			continue;

		double l = data_model.geodesics_candidates.lengths(flat_i);
		if (l > max_length)
		{
			max_length = l;
			new_vertex = flat_i;
		}
	}

	write_log(4) << "found new center vertex at " << new_vertex << ", prev. vertex was " << vertex << endl;
	vertex = new_vertex;
	geodesic = &data_model.geodesics_candidates.paths[vertex];
	*/

	//compute width of ruled surface based on size of flat area
	vector<int> flat_faces = meshhelper::get_faces_from_vertices(flat_vertices, data_model.target.adjacency_VF);
	double area = meshhelper::compute_face_area(flat_faces, data_model.target.V, data_model.target.F);
	const double area_scale = 2.0;
	//const double area_scale = 1.3;
	double length = data_model.geodesics_candidates.lengths[vertex];
	double width = (area * area_scale) / length;

	if (area < data_model.target.surface_area * 0.01)
	{
		write_log(loglevel) << "try get flat surface: not valid b/c found area is too small (" << area << "/" << data_model.target.surface_area << "), v" << original_vertex << "->" << vertex << endl;
		return false;
	}

	//visualize flat areas
	int number_labels = data_model.vertex_labels.maxCoeff();
	number_labels++;
	for (int i : flat_vertices)
		data_model.vertex_labels(i) = number_labels;


	//create surface
	out_surface = new RuledDevelopableSurface(data_model.target, *geodesic, vertex, RulingsType::Discrete);
	out_surface->create(width);

	//check that ruled surface is valid
	if (out_surface->developable.F.rows() < 2)
	{
		write_log(loglevel) << "try get flat surface: not valid b/c ruled surface computation is bad (too few faces)" << endl;
		return false;
	}


	assigned_indices.clear();
	assigned_indices = flat_vertices;

	this->flat_indices.push_back(vertex);
	this->flat_area_indices.push_back(flat_vertices);


	//debug output
	if (loglevel <= LOG_LEVEL) out_surface->print_mathematica_data();

	write_log(loglevel+1) << "flat_vertices: " << list_to_string(flat_vertices) << endl;
	write_log(4) << "(v" << original_vertex << " to center v" << vertex << ") area: " << area << ", length: " << length << ", width = (area * 2.0) / length: " << width << endl;
	//write_log(loglevel) << "F_areas: " << endl << F_areas << endl;
	//write_log(loglevel) << "flat_faces: "; log_list(loglevel, flat_faces, "", false);

	return true;
}

/* OLD try_get_flat_surface_at() */
/*
bool RuledGeodesicsController::try_get_flat_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices)
{
	if (vertex < 0 || vertex >= data_model.geodesics_candidates.lengths.rows())
		return false;


	//double sum = 0.0;
	//double squared_sum = 0.0;
	//for (double d : normal_deviations)
	//{
	//	sum += d;
	//	squared_sum += d * d;
	//}
	//double mean = sum / normal_deviations.size();
	//double varicance = (squared_sum / normal_deviations.size()) - mean * mean;
	//double stdev = sqrt(varicance);
	//if (stdev < 1e-2)
	//	stdev = 1e-2;
	//double lower_bound = mean - abs(stdev * stdev_factor);
	//double upper_bound = mean + abs(stdev * stdev_factor);
	//bool is_done = last_deviation > upper_bound || last_deviation < lower_bound;

	//vector<double> normal_deviations;
	//vector<int> face_neighbors = data_model.target.adjacency_VF[vertex];
	//for (int i = 0; i < face_neighbors.size() + 1; i++)
	//{
	//	int fi = face_neighbors[i];
	//	int fi_compare = face_neighbors[i+1 % face_neighbors.size()];
	//	
	//	auto n = data_model.target.normals_faces.row(fi);
	//	auto n_compare = data_model.target.normals_faces.row(fi_compare);
	//	
	//	double cos_angle = clip(n.dot(n_compare), -1, 1);
	//	double angle = acos(cos_angle);
	//	normal_deviations.push_back(angle);
	//}

	assigned_indices.clear();
	assigned_indices.push_back(vertex);

	auto vertex_normal = data_model.target.normals_vertices.row(vertex);
	for (int fi : data_model.target.adjacency_VF[vertex])
	{
		auto n = data_model.target.normals_faces.row(fi);
		double cos_angle = clip(n.dot(vertex_normal), -1, 1);
		double angle = acos(cos_angle);

		if (angle > angle_threshold)
			return false;
	}



	//out_flat_vertices.clear();
	//dfs_within_list(data_model.geodesics_candidates.flat_geodesics, data_model.target.adjacency_VV, vertex, out_flat_vertices);
	//vector<int> flat_faces = get_faces(out_flat_vertices, data_model.target.adjacency_VF);
	
	//vector<int> flat_faces;
	//dfs_flat_on_faces(data_model.target.adjacency_VF, data_model.target.normals_faces, data_model.target.adjacency_VF[vertex][0], flat_faces);
	
	vector<int> out_flat_vertices;
	dfs_flat_on_faces(data_model.target.adjacency_VV, data_model.target.normals_vertices, vertex, out_flat_vertices);

	if (out_flat_vertices.size() < 1)
		return false;

	vector<int> flat_faces = get_faces(out_flat_vertices, data_model.target.adjacency_VF);

	const int loglevel = 5;
	write_log(loglevel) << "out_flat_vertices: "; log_list(loglevel, out_flat_vertices, "", false);



	//TODO select central geodesic
	Eigen::MatrixXd flat_V;
	index_to_value(out_flat_vertices, data_model.target.V, flat_V);

	auto centroid = get_centroid(flat_V);
	int center_vertex = get_closest_vertex(data_model.target.V, centroid);
	int original_vertex = vertex; //for debug only
	vertex = center_vertex;

	Eigen::VectorXd F_areas = meshhelper::compute_face_areas(data_model.target.V, data_model.target.F);
	double area = 0.0;
	for (int fi : flat_faces)
		area += F_areas(fi);

	const double area_scale = 1.5;
	double length = data_model.geodesics_candidates.lengths(vertex);
	double width = (area * area_scale) / length;

	write_log(loglevel) << "(v" << original_vertex << " to center v" << vertex <<  ") area: " << area << ", length: " << length << ", width = (area * 1.5) / length: " << width << endl;
	//write_log(loglevel) << "F_areas: " << endl << F_areas << endl;
	//write_log(loglevel) << "flat_faces: "; log_list(loglevel, flat_faces, "", false);
	
	//Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[vertex];
	//if (geodesic.rows() < 6) //otherwise can't create spline of degree 5
	//	return false;

	out_surface = new RuledDevelopableSurface(data_model.target, data_model.geodesics_candidates.paths[vertex], vertex, RulingsType::Discrete);
	out_surface->create(width);
	if(loglevel <= LOG_LEVEL)
		out_surface->print_mathematica_data();

	if (out_surface->developable.F.rows() < 2) 
		return false;

	assigned_indices.clear();
	assigned_indices = out_flat_vertices;
	
	flat_indices.push_back(vertex);

	return true;
}
*/


void RuledGeodesicsController::create_ruled_constraints_selector()
{
	int number_selected = data_model.selected_ruled_vertex_indices.size();

	//get vertices per label
	vector<vector<int>> label_vertex_assignment(number_selected);

	for (int vi = 0; vi < data_model.vertex_ruled_assignment.size(); vi++)
	{
		int label = data_model.vertex_ruled_assignment[vi];
		int selected_index = index_of(label, data_model.selected_ruled_vertex_indices);
		label_vertex_assignment[selected_index].push_back(vi);
	}


	vector<Mesh> meshes;
	for (std::vector<int>& vertex_indices : label_vertex_assignment)
		meshes.push_back(meshhelper::sub_mesh(vertex_indices, data_model.target.V, data_model.target.F, data_model.target.adjacency_VF));


	selector = new MeritSelector(data_model.selected_ruled_vertex_indices, data_model.geodesics_candidates, meshes);
}

bool RuledGeodesicsController::get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex)
{
	if(!is_initialzed)
		return false;

	//meshhelper::update_coverage_per_vertices(data_model);
	//if (!data_model.target.uncovered_V.rows())
	//	return false;

	data_model.patch_coverage->compute(meshhelper::get_wrapper_meshes(data_model.patches));
	//if (data_model.patch_coverage->is_covered())
	//	return false;

	data_model.current_candidates_index++;
 	bool success = select_constraints(out_target_curve, out_wrapper_curve, vertex);

	if(!success)
		write_log(4) << "no more geodesics to fit" << endl;


	return success;
}

bool RuledGeodesicsController::select_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex)
{
	if (data_model.current_candidates_index < 0 || data_model.current_candidates_index >= data_model.selected_ruled_vertex_indices.size())
		return select_hole(out_target_curve, out_wrapper_curve, vertex);

	return select_geodesic(out_target_curve, out_wrapper_curve, vertex);
}

bool RuledGeodesicsController::select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex)
{
	selector->update_covering_mesh(meshhelper::concatenate_wrappers(data_model.patches));
	int index = selector->try_get_next();
	if (index == -1)
		return false;

	write_log(4) << endl << "MERIT selected index: " << index << " from: " << list_to_string(data_model.selected_ruled_vertex_indices) << endl;

	Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[index];
	if (geodesic.rows() < 3) //TODO should checked before
		return false;

	//resample_constraint_curves(geodesic, data_model.geodesics_candidates.lengths(index), out_target_curve, out_wrapper_curve);
	out_target_curve = geodesic;
	vertex = index;
	return true;

	/*
	double resampled_target_length;
	curve::resample_uniformly(geodesic, 1.0, out_target_curve, resampled_target_length);

	double target_off = abs(data_model.geodesics_candidates.lengths(index) - resampled_target_length);
	//write_log(0) << "original length: " << geodesics_candidates.geodesic_lengths(best_index) << ", resampled length: " << resampled_target_length << endl;
	//assert(target_off > 1e-2);
	if (target_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;


	double resampled_wrapper_length;
	curve::resample_matching_curve(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);

	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
	//assert(wrapper_off > 1e-2);
	if (wrapper_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;

	return true;
	*/
}

/*
bool RuledGeodesicsController::select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
{
	//add merit back in to weigh between lenghts and uncovered areas?

	if (data_model.current_candidates_index < 0 || data_model.current_candidates_index >= data_model.selected_ruled_vertex_indices.size())
	{
		return false;


		data_model.is_optimizing = false;
		data_model.optimization_settings.outlier_threshold *= 4;

		check_coverage_holes();
		if(data_model.hole_average_vertex.size() < 1)
			return false;

		//found holes that need coverage, add geodesics
		//TODO check if the geodesic is not 0
		//TODO leads to error in view because there are more selected_ruled_vertex_indices than selected_ruled_view_indices
		data_model.selected_ruled_vertex_indices.insert(data_model.selected_ruled_vertex_indices.end(), data_model.hole_average_vertex.begin(), data_model.hole_average_vertex.end());
		data_model.optimization_settings.outlier_threshold /= 4;
	}



	int index = data_model.selected_ruled_vertex_indices[data_model.current_candidates_index];
	Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[index];
	if (geodesic.rows() < 5) //TODO should checked before
		return false;

	double resampled_target_length;
	curve::resample_uniformly(geodesic, 1.0, out_target_curve, resampled_target_length);

	double target_off = abs(data_model.geodesics_candidates.lengths(index) - resampled_target_length);
	//write_log(0) << "original length: " << geodesics_candidates.geodesic_lengths(best_index) << ", resampled length: " << resampled_target_length << endl;
	//assert(target_off > 1e-2);
	if (target_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;


	double resampled_wrapper_length;
	curve::resample_matching_curve(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);

	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
	//assert(wrapper_off > 1e-2);
	if (wrapper_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;

	return true;
}
*/

bool RuledGeodesicsController::select_hole(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex)
{
	write_log(3) << "checking for uncovered areas..." << endl;

	double outlier_relax_factor = 1.0;
	data_model.is_optimizing = false;

	if (data_model.hole_average_vertex.size() < 1 && have_checked_holes)
		return false;

	if (data_model.hole_average_vertex.size() < 1)
	{
		data_model.optimization_settings.outlier_threshold *= outlier_relax_factor;
		check_coverage_holes();
	}

	write_log(3) << "...found " << data_model.hole_average_vertex.size() << " holes at indices: " << list_to_string(data_model.hole_average_vertex, ", ") << endl;


	if (data_model.hole_average_vertex.size() < 1)
	{
		data_model.optimization_settings.outlier_threshold /= outlier_relax_factor;
		return false;
	}

	selector->update_covering_mesh(meshhelper::concatenate_wrappers(data_model.patches));
	vertex = selector->try_get_next();
	if (vertex == -1)
		return false;

	int hole_mesh_index = index_of(vertex, data_model.hole_average_vertex);
	write_log(4) << endl << "MERIT selected hole vertex: " << vertex << " from: " << list_to_string(data_model.hole_average_vertex) << endl;


	auto target_constraints = data_model.hole_meshes[hole_mesh_index].V;
	const auto min_point = target_constraints.colwise().minCoeff();
	const auto max_point = target_constraints.colwise().maxCoeff();
	const auto dimensions = max_point - min_point;
	write_log(4) << "min_point: " << min_point << ", max_point: " << max_point << ", dimensions: " << dimensions << endl;

	const double scale = 1.0;
	const double min_length = 5.0;

	double max_dimension = dimensions.maxCoeff();
	double geodesic_length = max(max_dimension * scale, min_length);

	LengthGeodesicsStopping length_stopping;
	length_stopping.max_length = geodesic_length;
	GeodesicWalker walker(data_model.target, length_stopping);

	//TODO compute rotated direction from data_model.geodesics_directions
	Eigen::RowVector3d k_max = data_model.target.surface_features().principal_k_max_direction.row(vertex);

	Eigen::MatrixXd geodesic;
	double length;
	walker.get_geodesic_at(vertex, dimensions, geodesic, length);

	if (length < min_length)
		walker.get_geodesic_kmax_at(vertex, geodesic, length);
	if(length < min_length)
		walker.get_geodesic_kmin_at(vertex, geodesic, length);
	if (length < min_length)
		write_log(0) << "geodesic still too short!" << endl;

	//resample_constraint_curves(geodesic, length, out_target_curve, out_wrapper_curve);
	out_target_curve = geodesic;
	return true;
}

void resample_constraint_curves(const Eigen::MatrixXd& geodesic, const double length, Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
{
	//double resampled_target_length;
	//curve::resample_uniformly(geodesic, 1.0, out_target_curve, resampled_target_length);

	//double target_off = abs(length - resampled_target_length);
	//if (target_off > 1e-2)
	//	write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;


	//double resampled_wrapper_length;
	//curve::resample_matching_curve(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);

	//double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
	//if (wrapper_off > 1e-2)
	//	write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;
}

void RuledGeodesicsController::check_coverage_holes()
{
	double relax_factor = 1.0;
	//double relax_factor = 3.0;
	data_model.optimization_settings.coverage_threshold *= relax_factor;

	data_model.hole_vertex_indices.clear();
	data_model.hole_meshes.clear();
	data_model.hole_average_vertex.clear();

	data_model.patch_coverage->compute(meshhelper::get_wrapper_meshes(data_model.patches));
	if (data_model.patch_coverage->is_covered())
		return;

	data_model.hole_vertex_indices = label_vertices_uncovered_areas();
	data_model.hole_meshes = build_uncovered_area_meshes(data_model.hole_vertex_indices);

	for (Mesh mesh : data_model.hole_meshes)
	{
		Eigen::RowVector3d average = mesh.V.colwise().sum() / mesh.V.rows();
		int average_index = get_closest_vertex(data_model.target.V, average);
		data_model.hole_average_vertex.push_back(average_index);
	}

	if (data_model.hole_average_vertex.size() < 1)
		return;

	have_checked_holes = true;
	selector = new MeritSelector(data_model.hole_average_vertex, data_model.geodesics_candidates, data_model.hole_meshes);

	data_model.is_optimizing = false;
	data_model.optimization_settings.coverage_threshold /= relax_factor;
	write_log(4) << "...found " << data_model.hole_average_vertex.size() << " holes at indices: " << list_to_string(data_model.hole_average_vertex) << endl;
}

std::vector<std::vector<int>> RuledGeodesicsController::label_vertices_uncovered_areas()
{
	//meshhelper::update_coverage_per_vertices(data_model);

	//int number_uncovered_vertices = data_model.target.uncovered_V.rows();
	int number_uncovered_vertices = data_model.patch_coverage->uncovered_vertices().size();
	int count_labeled_vertices = 0;

	std::vector<bool> are_vertices_labeled(number_uncovered_vertices, false);
	std::vector<std::vector<int>> labeled_component_vertices;
	//std::vector<int> uncovered_V(data_model.target.uncovered_V.data(), data_model.target.uncovered_V.data() + data_model.target.uncovered_V.rows() * data_model.target.uncovered_V.cols());
	std::vector<int> uncovered_V = data_model.patch_coverage->uncovered_vertices();

	//log_list(4, uncovered_V, "uncovered_V: ", false);

	int number_labels = 0;
	//int start_vertex = data_model.target.uncovered_V(0);
	int start_vertex = uncovered_V[0];
	const int min_hole_size = 3;

	while (count_labeled_vertices < number_uncovered_vertices)
	{
		vector<int> component;
		dfs_within_list(uncovered_V, data_model.target.adjacency_VV, start_vertex, component);
		//log_list(4, component, "labeled component: ", false, false);


		count_labeled_vertices += component.size();
		for each (int labeled_vertex in component)
		{
			int index = index_of(labeled_vertex, uncovered_V);

			if (index >= 0) //-1 is the error case, when the value wasn't found
				are_vertices_labeled[index] = true;
		}

		auto iterator = find(are_vertices_labeled.begin(), are_vertices_labeled.end(), false);
		int index = distance(are_vertices_labeled.begin(), iterator);
		start_vertex = uncovered_V[index];

		if (component.size() < min_hole_size)
			continue;

		labeled_component_vertices.push_back(component);
		number_labels++;


		//write_log(4) << "number_labels: " << number_labels << ", count_labeled_vertices: " << count_labeled_vertices << endl;
		//write_log(4) << "index: " << index << ", start_vertex: " << start_vertex << endl << endl;
		//log_list(4, are_vertices_labeled, "are_vertices_labeled: ", false, false);
	}

	return labeled_component_vertices;
}

std::vector<Mesh> RuledGeodesicsController::build_uncovered_area_meshes(const std::vector<std::vector<int>>& labeled_component_vertices)
{
	//re-map faces to map to local uncovered area vertices
	std::vector<Eigen::MatrixXi> labeled_component_faces;

	const double min_area = 1.0;

	for each (vector<int> component in labeled_component_vertices)
	{
		vector<int> face_list = get_faces(component, data_model.target.adjacency_VF);
		double area = meshhelper::compute_face_area(face_list, data_model.target.V, data_model.target.F);
		write_log(0) << "hole:: uncovered area: " << area << ", #faces = " << face_list.size() << endl;

		/*
		vector<int> face_list;
		for each (int labeled_vertex in component)
		{
			vector<int> faces = data_model.target.adjacency_VF[labeled_vertex];
			face_list.insert(face_list.end(), faces.begin(), faces.end());
		}

		//log_list(4, face_list, "global face_list: ", false, false);

		//std::sort(face_list.begin(), face_list.end());
		//face_list.erase(std::unique(face_list.begin(), face_list.end()), face_list.end());
		face_list = remove_duplicates(face_list);

		//log_list(4, face_list, "global unique face_list: ", false, false);
		*/

		Eigen::MatrixXi component_faces(face_list.size(), 3);
		int added_faces = 0;

		for (int i = 0; i < face_list.size(); i++)
		{
			int fi = face_list[i];
			Eigen::RowVectorXi face_vertices_global = data_model.target.F.row(fi);
			//write_log(4) << i << ": global face_vertices for fi=" << fi << ": " << face_vertices_global << endl;

			vector<int> face_vertices_local;
			for (int j = 0; j < face_vertices_global.size(); j++)
			{
				int vi = face_vertices_global(j);
				face_vertices_local.push_back(index_of(vi, component));
			}

			if (is_contained(-1, face_vertices_local))
			{
				//write_log(4) << "   NOT including this face!" << endl;
				continue;
			}

			component_faces.row(added_faces) = Eigen::Map<Eigen::RowVectorXi>(face_vertices_local.data(), face_vertices_local.size());
			added_faces++;

			//write_log(4) << "local component_faces.row(" << i << "): " << component_faces.row(i) << endl;
		}

		component_faces.conservativeResize(added_faces, 3);
		labeled_component_faces.push_back(component_faces);
		//write_log(4) << "local component_faces: " << endl << component_faces << endl;
	}

	std::vector<Mesh> uncovered_area_meshes;
	for (int i = 0; i < labeled_component_vertices.size(); i++)
	{
		write_log(0) << "hole:: component has " << labeled_component_faces[i].rows() << " faces" << endl;

		//TODO use area instead of number of Faces, for irregular meshes!
		if (labeled_component_faces[i].rows() < 3)
		{
			write_log(0) << "--> SKIP" << endl;
			continue;
		}

		Eigen::MatrixXd V;
		index_to_value(labeled_component_vertices[i], data_model.target.V, V);

		Mesh mesh;
		mesh.V = V;
		mesh.F = labeled_component_faces[i];

		uncovered_area_meshes.push_back(mesh);
	}

	return uncovered_area_meshes;
}

std::vector<int> get_faces(const std::vector<int>& component_vertices, const std::vector<std::vector<int>>& adjacency_VF)
{
	std::vector<int> face_list;
	for each (int labeled_vertex in component_vertices)
	{
		std::vector<int> faces = adjacency_VF[labeled_vertex];
		face_list.insert(face_list.end(), faces.begin(), faces.end());
	}

	//log_list(4, face_list, "global face_list: ", false, false);

	//std::sort(face_list.begin(), face_list.end());
	//face_list.erase(std::unique(face_list.begin(), face_list.end()), face_list.end());
	face_list = remove_duplicates(face_list);

	return face_list;
}


/*
bool is_flat(double current_deviation, vector<double>& normal_deviations)
{
	const int window_size = 20;
	const double stdev_factor = 1.5;

	double sum = 0.0;
	double squared_sum = 0.0;

	for (double d : normal_deviations)
	{
		sum += d;
		squared_sum += d * d;
	}

	double mean = sum / normal_deviations.size();
	double varicance = (squared_sum / normal_deviations.size()) - mean * mean;
	double stdev = sqrt(varicance);

	if (stdev < 1e-2)
		stdev = 1e-2;

	double lower_bound = mean - abs(stdev * stdev_factor);
	double upper_bound = mean + abs(stdev * stdev_factor);

	bool flat = current_deviation > upper_bound || current_deviation < lower_bound;

	if (normal_deviations.size() >= window_size)
		normal_deviations.erase(normal_deviations.begin());

	normal_deviations.push_back(current_deviation);

	return flat;
}
*/


void dfs_flat_on_faces(const std::vector<std::vector<int>>& connectivity_list, const Eigen::MatrixXd& N, const int& start_index, std::vector<int>& out_component)
{
	// Initially mark all verices as not visited 
	std::vector<bool> visited(connectivity_list.size(), false);

	// Create a stack for DFS 
	std::stack<int> stack;

	// Push the current source node. 
	stack.push(start_index);
	Eigen::RowVectorXd mean_normal = N.row(start_index);


	int index;
	out_component.clear();


	while (!stack.empty())
	{
		// Pop a vertex from stack 
		index = stack.top();
		stack.pop();

		//std::vector<int>::iterator iterator = std::find(obstacles.begin(), obstacles.end(), index);
		//bool is_contained = iterator != obstacles.end();
		//if (is_contained)
		//	continue;

		// Stack may contain same vertex twice. So we need to print the popped item only if it is not visited. 
		if (!visited[index])
		{
			visited[index] = true;

			Eigen::RowVectorXd normal = N.row(index);
			//double cos_angle = clip(normal.dot(mean_normal), -1, 1);
			//double angle = acos(cos_angle);
			////write_log(6) << "mean_normal: " << mean_normal << ", normal: " << normal << endl;
			////write_log(6) << "dot: " << normal.dot(mean_normal) << ", cos_angle: " << cos_angle << ", angle: " << angle << endl;

			if (angle(normal.transpose(), mean_normal.transpose()) <= GlobalSettings::flat_angle_threshold)
			{
				out_component.push_back(index);
				mean_normal = (mean_normal + normal).normalized();
			}
		}

		// Get all adjacent vertices of the popped vertex s 
		// If a adjacent has not been visited, then push it to the stack. 
		for (auto i = connectivity_list[index].begin(); i != connectivity_list[index].end(); ++i)
			if (!visited[*i])
				stack.push(*i);
	}
}



/*
bool RuledGeodesicsController::select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
{
	int geodesic_index = random_geodesics[data_model.patches.size()];
	Eigen::MatrixXd original_geodesic = data_model.geodesics_candidates.paths[geodesic_index];
	
	//resample geodesic on target to match wrapper vertices
	double resampled_target_length;
	curve::resample_uniformly(original_geodesic, out_target_curve, resampled_target_length);

	double target_off = abs(data_model.geodesics_candidates.lengths(geodesic_index) - resampled_target_length);

	if (target_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;


	double resampled_wrapper_length;
	curve::resample_matching_curve(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);

	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
	if(wrapper_off > 1e-2)
		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;

	return true;
}
*/
