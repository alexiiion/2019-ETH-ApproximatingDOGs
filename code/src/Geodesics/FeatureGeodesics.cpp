//#include "FeatureGeodesics.h"
//
//#include <igl/boundary_loop.h>
//#include <igl/exact_geodesic.h>
//
//#include "CoordinateConverter.h"
//#include "Logger.h"
//#include "Utils.h"
//
//using namespace std;
//
//void FeatureGeodesics::get_all_geodesics(DataModel& data_model, std::vector<int> sources, std::vector<std::vector<int>> targets)
//{
//	sources.pop_back(); //remove the last entry, which is -1 and signifies the end of the cluster
//
//	/*
//	vector<int> concatenated_targets = sources;
//	for (int cluster_index = 0; cluster_index < targets.size(); cluster_index++)
//	{
//		auto& feature_vertices = targets[cluster_index];
//		concatenated_targets.insert(concatenated_targets.end(), feature_vertices.begin(), feature_vertices.end() - 1);
//	}
//
//	for (int i = 0; i < sources.size(); i++)
//	{
//		std::vector<Eigen::MatrixXd> current_paths;
//		Eigen::VectorXd current_lengths;
//
//		find_geodesics_on_target(data_model.target, vector<int>{sources[i]}, concatenated_targets, current_paths, current_lengths);
//
//		append_vector(current_lengths, lengths);
//		paths.insert(paths.end(), current_paths.begin(), current_paths.end());
//	}
//	*/
//
//	/*
//	for (int cluster_index = 0; cluster_index < targets.size(); cluster_index++)
//	{
//		auto feature_vertices = targets[cluster_index];
//		feature_vertices.pop_back();
//
//		std::vector<Eigen::MatrixXd> current_paths;
//		Eigen::VectorXd current_lengths;
//
//		find_geodesics_on_target(data_model.target, sources, feature_vertices, current_paths, current_lengths);
//
//		append_vector(current_lengths, lengths);
//		paths.insert(paths.end(), current_paths.begin(), current_paths.end());
//	}
//
//	preprocess_target_geodesics(data_model.target);
//	compute_constant_features(data_model.feature_settings);
//	*/
//
//
//	//vector<int> concatenated_targets = sources;
//	vector<int> concatenated_targets;
//
//	for (int cluster_index = 0; cluster_index < targets.size(); cluster_index++)
//	{
//		auto& feature_vertices = targets[cluster_index];
//		concatenated_targets.insert(concatenated_targets.end(), feature_vertices.begin(), feature_vertices.end() - 1);
//	}
//
//	find_geodesics_on_target(data_model.target, sources, concatenated_targets, paths, lengths);
//	
//	number_geodesics = paths.size();
//	geodesics_merit.resize(number_geodesics);
//	geodesic_features.resize(number_geodesics, data_model.feature_settings.number_geodesic_features);
//
//	if (number_geodesics < 1)
//		return;
//
//	preprocess_target_geodesics(data_model.target);
//	compute_constant_features(data_model.feature_settings);
//
//}
//
//void FeatureGeodesics::get_all_geodesics(DataModel& data_model, std::vector<std::vector<int>>& feature_connectivity_graph)
//{
//	if (feature_connectivity_graph.size() < 1)
//		return;
//
//	write_log(3) << endl << "start finding all geodesics of feature cluster graph" << endl << endl;
//
//	for (int i = 0; i < feature_connectivity_graph.size(); i++)
//	{
//		vector<int> connected_clusters = feature_connectivity_graph[i];
//
//		std::vector<int> sources = data_model.clustered_feature_vertices[i];
//		sources.pop_back(); //remove the last entry, which is -1 and signifies the end of the cluster
//
//		vector<int> concatenated_targets;
//
//		for (int cluster_index = 0; cluster_index < connected_clusters.size(); cluster_index++)
//		{
//			auto& feature_vertices = data_model.clustered_feature_vertices[connected_clusters[cluster_index]];
//			concatenated_targets.insert(concatenated_targets.end(), feature_vertices.begin(), feature_vertices.end() - 1);
//		}
//
//		std::vector<Eigen::MatrixXd> current_paths;
//		Eigen::VectorXd current_lengths;
//
//		find_geodesics_on_target(data_model.target, sources, concatenated_targets, current_paths, current_lengths);
//		//TODO remove geodesics with length < 0.1
//
//		append_vector(current_lengths, lengths);
//		paths.insert(paths.end(), current_paths.begin(), current_paths.end());
//
//		//write_log(0) << "  [" << i << "]: " << endl;
//		//for (int j = 0; j < current_paths.size(); j++)
//		//	write_log(0) << "  path " << j << ": " << endl << current_paths[j] << endl << endl;
//		//write_log(0) << "   lengths " << endl << current_lengths << endl << endl;
//	}
//
//	write_log(3) << endl << "done. total number geodesics found: " << paths.size() << endl << endl;
//
//	preprocess_target_geodesics(data_model.target);
//	compute_constant_features(data_model.feature_settings);
//}
//
//void FeatureGeodesics::find_geodesics_on_target(const TargetMesh& target, std::vector<int>& sources, std::vector<int>& targets, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths)
//{
//	//write_log(3) << endl << "- find all GEODESICS on target" << endl << endl;
//	if (sources.size() < 1 || targets.size() < 1)
//		return;
//
//	Eigen::VectorXi source_faces, target_faces;
//	Eigen::VectorXi source_indices = Eigen::Map<Eigen::VectorXi>(sources.data(), sources.size());
//	Eigen::VectorXi target_indices = Eigen::Map<Eigen::VectorXi>(targets.data(), targets.size());
//	
//	std::vector<Eigen::MatrixXd> all_paths;
//	Eigen::VectorXd all_lengths;
//	igl::exact_geodesic(target.V, target.F, source_indices, source_faces, target_indices, target_faces, all_lengths, all_paths);
//
//	//omit very short geodesics because we cannot do curve interpolation on them
//	out_geodesic_lengths.resize(all_lengths.rows());
//
//	for (int i = 0; i < all_lengths.rows(); i++)
//	{
//		if (all_lengths(i) < 4.0)
//			continue;
//
//		out_geodesic_lengths(out_geodesic_paths.size()) = all_lengths(i);
//		out_geodesic_paths.push_back(all_paths[i]);
//	}
//	out_geodesic_lengths.conservativeResize(out_geodesic_paths.size());
//
//
//	write_log(3) << "  geodesic: number geodesics found: " << out_geodesic_paths.size() << endl;
//	write_log(6) << "geodesic: source_indices (number = " << sources.size() << "):" << endl << source_indices << endl;
//	write_log(6) << "geodesic: target_indices (number = " << targets.size() << "):" << endl << target_indices << endl;
//	write_log(6) << "geodesic: out_geodesic_distances: " << endl << out_geodesic_lengths << endl;
//
//	//for (int i = 0; i < geodesic_paths.size(); i++)
//	//	debug_show_curve(geodesic_paths[i], data_model);
//}
//
////void get_all_geodesics(DataModel& data_model, std::vector<int>& segment_indices, int source_index, FeatureGeodesics& segment_geodesics)
////{
////	if (!data_model.target.V.rows())
////	{
////		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> target mesh is not set, cannot find geodesics!" << endl << endl;
////		return;
////	}
////
////	if (segment_indices.size() < 1)
////	{
////		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> target indices are empty, cannot find geodesics!" << endl << endl;
////		return;
////	}
////
////	find_geodesics_on_target(data_model, segment_indices, source_index, segment_geodesics.geodesic_paths, segment_geodesics.geodesic_lengths);
////	preprocess_target_geodesics(data_model, segment_geodesics);
////	compute_constant_features(data_model.feature_settings, segment_geodesics);
////}
//
//void FeatureGeodesics::compute_constant_features(const GeodesicFeatureConfig& feature_settings)
//{
//	const int loglevel = 5;
//	int current_feature_index = 0;
//
//	const int number_geodesics = paths.size();
//	geodesic_features.resize(number_geodesics, feature_settings.number_geodesic_features);
//	
//	add_normalized_constant_feature(current_feature_index++, gauss_curvatures, geodesic_features);
//	add_normalized_constant_feature(current_feature_index++, lengths, geodesic_features);
//	add_normalized_constant_feature(current_feature_index++, curvatures, geodesic_features);
//
//	write_log(loglevel) << "geodesic_constant_features " << endl << geodesic_features << endl;
//}
//
//void FeatureGeodesics::add_normalized_constant_feature(const int current_feature_index, const Eigen::VectorXd& feature, Eigen::MatrixXd& out_features)
//{
//	double max = feature.maxCoeff();
//
//	if (abs(max) < 1e-6) //avoid division by zero
//		max = 1e-6;
//
//	for (int i = 0; i < feature.rows(); i++)
//		out_features(i, current_feature_index) = feature(i) / max;
//}
//
////void FeatureGeodesics::find_geodesics_on_target(DataModel& data_model, std::vector<int>& segment_indices, int& source_index, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths)
////{
////	write_log(3) << endl << "- find all GEODESICS on target" << endl << endl;
////	if (segment_indices.size() < 1)
////		return;
////
////	//source_index = segment_indices[31];
////	if(source_index == -1) //pick random vertex from target vertices
////	{
////		srand(time(NULL)); //TODO keep instance of rand
////		int random_index = rand() % segment_indices.size();
////		source_index = segment_indices[random_index];
////		write_log(4) << "geodesic: random source selected: " << source_index << endl;
////	}
////
////	Eigen::VectorXi source_indices, source_faces, target_faces;
////	source_indices.resize(1);
////	source_indices << source_index;
////
////	Eigen::VectorXi target_indices = Eigen::Map<Eigen::VectorXi>(segment_indices.data(), segment_indices.size());
////	igl::exact_geodesic(data_model.target.V, data_model.target.F, source_indices, source_faces, target_indices, target_faces, out_geodesic_lengths, out_geodesic_paths);
////
////	write_log(4) << "geodesic: source index: " << source_index << endl;
////	write_log(4) << "geodesic: number geodesics found: " << out_geodesic_paths.size() << endl;
////	write_log(6) << "geodesic: source_indices: " << endl << source_indices << endl;
////	write_log(6) << "geodesic: target_indices: " << endl << target_indices << endl;
////	write_log(6) << "geodesic: out_geodesic_distances: " << endl << out_geodesic_lengths << endl;
////
////	//for (int i = 0; i < geodesic_paths.size(); i++)
////	//	debug_show_curve(geodesic_paths[i], data_model);
////}
//
//void FeatureGeodesics::preprocess_target_geodesics(const TargetMesh& target)
//{
//	const int loglevel = 5;
//	write_log(loglevel - 1) << endl << "preprocessing geodesics..." << endl;
//
//	//coverage.resize(number_geodesics);
//	//coverage.setConstant(100000);
//
//	remove_straight_segments(lengths, paths);
//	compute_gauss_curvatures(target, paths, lengths, gauss_curvatures);
//	compute_curvatures(paths, curvatures);
//
//	create_lookup(gauss_curvatures, gauss_curvatures_lookup);
//	create_lookup(curvatures, curvatures_lookup);
//	create_lookup(lengths, lengths_lookup);
//
//#pragma region only_debug_logging
//	write_log(loglevel) << endl << "geodesic_curvatures: " << endl << curvatures << endl << endl;
//	write_log(loglevel) << endl << "geodesic_curvatures_lookup: " << endl << curvatures_lookup << endl << endl;
//
//	write_log(loglevel) << endl << "geodesic_lengths: " << endl << lengths << endl << endl;
//	write_log(loglevel) << endl << "geodesic_lengths_lookup: " << endl << lengths_lookup << endl << endl;
//
//	int max_curvature_i = curvatures_lookup(0);
//	double max_curvature = curvatures(max_curvature_i);
//	write_log(loglevel) << "max curvature: at " << max_curvature_i << " = " << max_curvature << endl;
//	//debug_show_curve(geodesic_paths_target[max_curvature_i], data_model, true);
//
//	int max_length_i = lengths_lookup(0);
//	double max_length = lengths(max_length_i);
//	write_log(loglevel) << "max distance: at " << max_length_i << " = " << max_length << endl;
//	//debug_show_curve(geodesic_paths_target[max_length_i], data_model, true);
//#pragma endregion
//
//	write_log(loglevel - 1) << "...done preprocessing geodesics" << endl << endl;
//}
//
//void FeatureGeodesics::remove_straight_segments(const Eigen::VectorXd& geodesic_lengths, std::vector<Eigen::MatrixXd>& geodesic_paths)
//{
//	const int loglevel = 5;
//
//	const double epsilon_angle = 1e-6;
//	const double min_length = 1e-3;
//	//const double min_length = 0.1;
//	const int number_geodesics = geodesic_paths.size();
//
//	for (int i = 0; i < number_geodesics; i++)
//	{
//		Eigen::MatrixXd path = geodesic_paths[i];
//		const int number_points = path.rows();
//		const double max_length = geodesic_lengths(i) / 5.0;
//
//		write_log(loglevel) << "  geodesic[" << i << "]" << endl;
//		//write_log(loglevel) << "  geodesic[" << i << "] (len=" << geodesic_lengths(i) << ", num=" << path.rows() << "): " << endl;
//
//		if (number_points < 2)
//			continue;
//
//		Eigen::MatrixXd path_processed(number_points, 3);
//		path_processed.row(0) = path.row(0); //add first point of geodesic
//
//		int count = 1;
//
//		for (int j = 1; j < number_points - 1; j++)
//		{
//			double segment_length = (path.row(j + 1) - path.row(j - 1)).squaredNorm();
//			if (segment_length < min_length)
//				continue;
//
//			double gap_length = (path_processed.row(count - 1) - path.row(j)).squaredNorm();
//			if (gap_length >= max_length)
//			{
//				path_processed.row(count) = path.row(j);
//				count++;
//
//				write_log(loglevel + 1) << "  --> add vector!" << endl;
//				continue;
//			}
//
//			Eigen::RowVector3d e1 = (path.row(j) - path.row(j - 1)).normalized();
//			Eigen::RowVector3d e2 = (path.row(j + 1) - path.row(j)).normalized();
//
//			double dot = e1.dot(e2);
//			double cos_angle = clip(dot, -1, 1);
//
//			if (cos_angle == 1 || cos_angle == -1)
//				continue;
//
//			double angle = acos(cos_angle);
//
//			write_log(loglevel + 1) << "  point [" << j << "]: " << endl;
//			write_log(loglevel + 1) << "    cos(theta) = " << cos_angle << endl;
//
//			if (abs(angle) < epsilon_angle)
//				continue;
//
//			path_processed.row(count) = path.row(j);
//			count++;
//
//			write_log(loglevel + 1) << "  --> add vector!" << endl;
//		}
//
//		path_processed.row(count) = path.row(number_points - 1); //add last point of geodesic
//
//		path_processed.conservativeResize(count + 1, Eigen::NoChange);
//		geodesic_paths[i] = path_processed;
//
//		write_log(loglevel) << "    points from " << number_points << " to " << path_processed.rows() << endl;
//	}
//}
//
//void FeatureGeodesics::compute_gauss_curvatures(const TargetMesh& target, const std::vector<Eigen::MatrixXd>& geodesic_paths, const Eigen::VectorXd& geodesic_lengths, Eigen::VectorXd& out_gaussian_curvatures)
//{
//	const int loglevel = 6;
//
//	out_gaussian_curvatures.resize(geodesic_paths.size());
//
//	for (int i = 0; i < geodesic_paths.size(); i++)
//	{
//		auto geodesic = geodesic_paths[i];
//		const int number_points = geodesic.rows();
//		if (number_points < 1)
//		{
//			out_gaussian_curvatures(i) = 0;
//			continue;
//		}
//
//		Eigen::MatrixXi bary_indices;
//		Eigen::MatrixXd bary_weigths;
//		CoordinateConverter::barycentric_coords_from_points(geodesic, target.V, target.F, bary_indices, bary_weigths);
//
//		write_log(loglevel) << endl << "geodesic_path[" << i << "]: " << endl <<
//			"bary_indices: " << endl << bary_indices << endl <<
//			"bary_weigths: " << endl << bary_weigths << endl <<
//			"gaussian_curvature_target: " << endl << target.surface_features.K << endl <<
//			endl;
//
//		//get all K for geodesic path points
//		Eigen::VectorXd K_per_point(number_points);
//		for (int j = 0; j < number_points; j++)
//		{
//			double k = 0;
//
//			for (int bary_i = 0; bary_i < 3; bary_i++)
//			{
//				double weight = bary_weigths(j, bary_i);
//				int vertex_index = bary_indices(j, bary_i);
//				double vertex_K = abs(target.surface_features.K(vertex_index));
//
//				k += weight * vertex_K;
//			}
//			K_per_point(j) = k;
//		}
//		write_log(loglevel) << "K_per_point:  " << endl << K_per_point << endl << endl;
//
//
//		double curvature = 0.0;
//
//		for (int j = 1; j < number_points - 1; j++)
//		{
//			//TODO somehow normalize for lengths
//			//double prev_K = K_per_point(j - 1);
//			//double next_K = K_per_point(j + 1);
//
//			curvature += K_per_point(j);
//		}
//
//		// add K for first point
//		curvature += K_per_point(0);
//
//		// add K for last point
//		curvature += K_per_point(number_points - 1);
//
//		////apply this on sum
//		//if (geodesic_lengths(i) > 1e-6)
//		//{
//		//	const double length_normalization = 1 / geodesic_lengths(i);
//		//	curvature *= length_normalization;
//		//}
//		//else
//		//{
//		//	curvature = 0;
//		//}
//
//		write_log(loglevel) << i << " gauss curvature:  " << curvature << endl;
//
//		out_gaussian_curvatures(i) = curvature;
//	}
//
//	write_log(loglevel-1) << "gauss curvatures:  " << endl << out_gaussian_curvatures << endl << endl;
//}
//
//void FeatureGeodesics::compute_curvatures(const std::vector<Eigen::MatrixXd>& geodesic_paths, Eigen::VectorXd & out_curvatures)
//{
//	out_curvatures.resize(geodesic_paths.size());
//
//	for (int i = 0; i < geodesic_paths.size(); i++)
//	{
//		double curvature = 0;
//
//		//Curve c(paths[i]);
//		//for (double k : c.k)
//		//	curvature += k;
//
//		compute_curvature(geodesic_paths[i], curvature);
//		out_curvatures(i) = curvature;
//
//		write_log(5) << endl << "  total curvature = " << curvature << endl;
//	}
//}
//
//void FeatureGeodesics::compute_curvature(const Eigen::MatrixXd& path, double& out_curvature)
//{
//	const int loglevel = 6;
//
//	const int number_points = path.rows();
//	out_curvature = 0.0;
//
//	write_log(loglevel) << "curvatures:  " << endl;
//	for (int i = 0; i < number_points - 2; i++)
//	{
//		Eigen::RowVector3d e1 = (path.row(i + 1) - path.row(i)).normalized();
//		Eigen::RowVector3d e2 = (path.row(i + 2) - path.row(i + 1)).normalized();
//
//		double dot = e1.dot(e2);
//		double cos_angle = clip(dot, -1, 1);
//
//		if (cos_angle == 1 || cos_angle == -1)
//			continue;
//
//		double angle = acos(cos_angle);
//		double hypothenuse = (path.row(i + 2) - path.row(i)).norm();
//		double curvature = 2 * sin(angle) / hypothenuse;
//
//		if (hypothenuse < 1e-3) //too short hypothenuse creates huge curvature. but with dense meshes this happens... 
//			curvature = 0;
//
//		out_curvature += curvature;
//		write_log(loglevel) << curvature << ", ";
//	}
//	write_log(loglevel) << endl;
//}
