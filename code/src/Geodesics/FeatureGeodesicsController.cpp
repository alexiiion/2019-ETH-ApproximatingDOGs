//#include "FeatureGeodesicsController.h"
//
//#include <igl/point_mesh_squared_distance.h>
//#include <igl/exact_geodesic.h>
//#include <igl/triangle_triangle_adjacency.h>
//#include <igl/per_vertex_normals.h>
//#include <igl/per_face_normals.h>
//
//#include "PatchModel.h"
//#include "MeshController.h"
//#include "GeodesicCandidates.h"
//
//#include "CurveHelper.h"
//#include "Utils.h"
//#include "Logger.h"
//
//using namespace std;
//
//
//const Eigen::RowVector3d sample_direction = Eigen::RowVector3d::UnitX();
//
//void scale_target(DataModel& data_model);
//bool is_probe_intersecting_features(const DataModel& data_model, const vector<int>& all_close_vertices, const int feature_index, const int candidate_feature_index);
//void get_closest_vertices(const Eigen::MatrixXd& V, const Eigen::MatrixXd& path, std::vector<int>& out_closest_vertices);
//
//bool FeatureGeodesicsController::initialize()
//{
//	if (is_initialzed)
//		return true;
//
//	if (data_model.clustered_feature_vertices.size() == 0)
//	{
//		write_log(2) << endl << "WARNING in <" << __FUNCTION__ << "> feature connectivity graph could not be built and thus optimiztaion could not be initialized. were feature clusters defined?" << endl << endl;
//		return false;
//	}
//
//	data_model.current_candidates_index = 0;
//	//data_model.current_candidates_index = 19; //TODO remove!! bunny nose
//
//	if (data_model.feature_connectivity_graph.size() == 0)
//		build_feature_graph();
//
//	if (data_model.feature_geodesics_candidates.size() == 0)
//		build_geodesics();
//	else
//		update_geodesics_merit(data_model.current_candidates_index);
//
//	if (data_model.feature_geodesics_candidates.size() == 0)
//		return false;
//
//	//if (is_equal(data_model.target.scale, 1.0))
//	//	scale_target(data_model);
//
//
//	is_initialzed = true;
//	return is_initialzed;
//}
//
//void FeatureGeodesicsController::build_feature_graph()
//{
//	write_log(3) << endl << "creating feature graph..." << endl;
//
//	const int loglevel = 0;
//	const int number_probes = 6;
//
//	for (int feature_index = 0; feature_index < data_model.clustered_feature_vertices.size(); feature_index++)
//	{
//		vector<int> connected_clusters;
//
//		vector<int> current_cluster = data_model.clustered_feature_vertices[feature_index];
//		if (current_cluster.size() < 1)
//			continue;
//
//		vector<int> source_indices = current_cluster;
//		source_indices.pop_back();
//
//		for (int candidate_feature_index = 0; candidate_feature_index < data_model.clustered_feature_vertices.size(); candidate_feature_index++)
//		{
//			write_log(loglevel) << endl << endl << "testing feature " << feature_index << " with feature " << candidate_feature_index << endl;
//
//			if (feature_index == candidate_feature_index)
//				continue;
//
//			////for already processed clusters, check if connection already exists to void duplicate connections (undirected graph)
//			//if (candidate_feature_index < feature_index)
//			//{
//			//	vector<int> graph_neighborhood = data_model.feature_connectivity_graph[candidate_feature_index];
//			//	vector<int>::iterator iterator = find(graph_neighborhood.begin(), graph_neighborhood.end(), feature_index);
//			//	
//			//	bool is_contained = iterator != graph_neighborhood.end();
//			//	if (is_contained)
//			//		continue;
//			//}
//			
//			vector<int> candidate_cluster = data_model.clustered_feature_vertices[candidate_feature_index];
//			if (candidate_cluster.size() < 1)
//				continue;
//
//
//			//get geodesic probes to test for intersections
//			write_log(loglevel) << "  find <= " << number_probes << " probe geodesics" << endl;
//
//			vector<int> target_indices;
//			if (candidate_cluster.size()-1 < number_probes)
//			{
//				target_indices = candidate_cluster;
//				target_indices.pop_back(); // remove -1
//			}
//			else
//			{
//				int step = (candidate_cluster.size()-1) / number_probes;
//				int max = step * (number_probes - 1);
//
//				for (int i = 0; i < max; i += step)
//					target_indices.push_back(candidate_cluster[i]);
//				target_indices.push_back(candidate_cluster[candidate_cluster.size() - 2]);
//			}
//
//			log_list(loglevel, target_indices, "  ...find probes to target indices: ", false);
//
//			std::vector<Eigen::MatrixXd> probe_paths;
//			Eigen::VectorXd probe_lengths;
//
//			//FeatureGeodesics geodesics;
//			//geodesics.find_geodesics_on_target(data_model.target, source_indices, target_indices, probe_paths, probe_lengths);
//			find_geodesics_on_target(data_model.target, source_indices, target_indices, probe_paths, probe_lengths);
//
//			write_log(loglevel) << "  get all closest vertices to path points" << endl;
//			int number_probes_intersecting = 0;
//
//			//get all closest vertices to path points
//			for (int probe_index = 0; probe_index < probe_paths.size(); probe_index++)
//			{
//				Eigen::MatrixXd& probe_path = probe_paths[probe_index];
//
//				vector<int> all_close_vertices;
//				get_closest_vertices(data_model.target.V, probe_path, all_close_vertices);
//				
//				//keep only unique ones
//				std::sort(all_close_vertices.begin(), all_close_vertices.end());
//				all_close_vertices.erase(std::unique(all_close_vertices.begin(), all_close_vertices.end()), all_close_vertices.end());
//				//log_list(loglevel, all_close_vertices, "  ...found all close vertices: ", false);
//
//				write_log(loglevel) << "  test intersections with all features" << endl;
//
//				bool is_intersecting = is_probe_intersecting_features(data_model, all_close_vertices, feature_index, candidate_feature_index);
//				if(is_intersecting)
//					number_probes_intersecting++;
//			}
//			write_log(loglevel) << "  ...number_probes_intersecting: " << number_probes_intersecting << endl;
//
//			if (number_probes_intersecting > 1)
//				continue;
//
//			connected_clusters.push_back(candidate_feature_index);
//			data_model.number_connections++;
//
//			write_log(loglevel) << "features " << feature_index << " & " << candidate_feature_index << " are CONNECTED!" << endl;
//		}
//
//		//if(connected_clusters.size() > 0)
//			data_model.feature_connectivity_graph.push_back(connected_clusters);
//	}
//
//	if (LOG_LEVEL <= 4)
//	{
//		write_log(0) << endl << "cluster graph: " << endl;
//		for (int i = 0; i < data_model.feature_connectivity_graph.size(); i++)
//		{
//			write_log(0) << i << ": ";
//			log_list(0, data_model.feature_connectivity_graph[i], "", false);
//		}
//	}
//
//	write_log(3) << endl << "...done creating feature graph" << endl;
//}
//
//void FeatureGeodesicsController::build_geodesics()
//{
//	if (data_model.feature_connectivity_graph.size() < 1)
//	{
//		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> feature connectivity graph not built yet!" << endl << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	for (int i = 0; i < data_model.feature_connectivity_graph.size(); i++)
//	{
//		write_log(4) << endl << endl << "cluster " << i << ": " << endl;
//		vector<int> cluster_neighbors = data_model.feature_connectivity_graph[i];
//
//		vector<vector<int>> targets;
//		for (int j = 0; j < cluster_neighbors.size(); j++)
//			targets.push_back(data_model.clustered_feature_vertices[cluster_neighbors[j]]);
//
//		//FeatureGeodesics geodesics_candidates;
//		//geodesics_candidates.get_all_geodesics(data_model, data_model.clustered_feature_vertices[i], targets);
//
//		GeodesicCandidates geodesics_candidates = get_all_geodesics(data_model, data_model.clustered_feature_vertices[i], targets);
//		data_model.geodesics_candidates = geodesics_candidates;
//	}
//
//	//TODO calc features merit: prioritize feature neighborhood with (1) little coverage, (1) few source indices, (2) few connections, (3) highest merit (in this order of importance)
//	
//
//	//Eigen::MatrixXi TT;
//	//igl::triangle_triangle_adjacency(data_model.target.F, TT);
//
//	//Eigen::MatrixXd NF;
//	//igl::per_face_normals(data_model.target.V, data_model.target.F, NF);
//
//	////get source face
//	//Eigen::MatrixXd NV;
//	//igl::per_vertex_normals(data_model.target.V, data_model.target.F, NV);
//
//	//FeatureGeodesics geodesics_max;
//	//for (int i = 0; i < data_model.target.V.rows(); i++)
//	//{
//	//	auto k_max_direction = data_model.target.V.row(i) + data_model.target.surface_features.principal_k_max_direction.row(i);
//	//	int source_face = get_source_face(data_model.target.V, data_model.target.F, i, k_max_direction, NV, data_model.target.adjacency_VV, data_model.target.adjacency_VF);
//	//	auto path = trace_geodesic(data_model.target.V, data_model.target.F, TT, NF, source_face, k_max_direction);
//
//	//	geodesics_max.paths.push_back(path);
//	//}
//}
//
//bool FeatureGeodesicsController::get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
//{
//	if(!is_initialzed)
//		return false;
//
//	if (!data_model.target.uncovered_V.rows())
//		return false;
//
//	
//	write_log(3) << endl << "get next geodesic as target constraints:" << endl;
//	update_geodesics_merit(data_model.current_candidates_index);
//
//	//if the current feature is covered, search forward
//	while (is_feature_covered)
//	{
//		write_log(5) << "geodesic_features: " << endl << data_model.feature_geodesics_candidates[data_model.current_candidates_index].geodesic_features << endl;
//		write_log(5) << "geodesics_merit: " << endl << data_model.feature_geodesics_candidates[data_model.current_candidates_index].geodesics_merit << endl;
//
//		data_model.current_candidates_index++;
//		if (data_model.current_candidates_index >= data_model.feature_geodesics_candidates.size())
//			return false;
//
//		//can't find geodesic in next feature neighborhood because it has no neighbors. skip to next feature
//		if (data_model.feature_geodesics_candidates[data_model.current_candidates_index].number_geodesics < 1)
//			continue;
//
//		is_feature_covered = false;
//		update_geodesics_merit(data_model.current_candidates_index);
//	}
//
//	select_geodesic(out_target_curve, out_wrapper_curve);
//
//	/*
//	int curve_length = 0;
//	while (curve_length < 4)
//	{
//		bool success = select_geodesic(out_target_curve, out_wrapper_curve);
//		if (!success)
//		{
//			is_feature_covered = true;
//			return get_next_constraints(out_target_curve, out_wrapper_curve);
//		}
//
//		curve_length = out_target_curve.rows();
//	}
//	*/
//
//	return true;
//}
//
//bool FeatureGeodesicsController::select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
//{
//	auto& current_candidates = data_model.feature_geodesics_candidates[data_model.current_candidates_index];
//	if (current_candidates.number_geodesics < 1)
//		return false;
//
//	int best_index;
//	double max_merit = current_candidates.geodesics_merit.maxCoeff(&best_index);
//	write_log(4) << endl << "selected geodesic: " << best_index << " (merit = " << max_merit << ")" << endl;
//
//	if (max_merit == -1)
//		return false;
//
//	current_candidates.current_geodesic_index = best_index;
//	current_candidates.current_merit = max_merit;
//
//	//out_target_curve = geodesics_candidates.geodesic_paths[best_index];
//	//double resampled_target_length = geodesics_candidates.geodesic_lengths(best_index);
//
//	//resample geodesic on target to match wrapper vertices
//	double resampled_target_length;
//	resample_target_geodesic(best_index, out_target_curve, resampled_target_length);
//	
//	if (current_candidates.lengths(best_index) < 4.0)
//	{
//		write_log(0) << "select geodesic < 4.0::lengths: " << endl << current_candidates.lengths << endl;
//		return false;
//	}
//
//	//write_log(0) << "out_target_curve: " << endl << out_target_curve << endl;
//
//	//write_log(0) << "out_target_curve -- check segment lengths: " << endl;
//	//for(int i = 1; i < out_target_curve.rows(); i++)
//	//	write_log(0) << i << ": " << (out_target_curve.row(i)-out_target_curve.row(i-1)).norm() << " --> (" << out_target_curve.row(i) << ") - (" << out_target_curve.row(i-1) << ")" << endl;
//
//	double target_off = abs(current_candidates.lengths(best_index) - resampled_target_length);
//	//write_log(0) << "original length: " << geodesics_candidates.geodesic_lengths(best_index) << ", resampled length: " << resampled_target_length << endl;
//	//assert(target_off > 1e-2);
//	if (target_off > 1e-2)
//		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;
//
//
//	double resampled_wrapper_length;
//	resample_wrapper_geodesic(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);
//
//	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
//	//assert(wrapper_off > 1e-2);
//	if(wrapper_off > 1e-2)
//		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;
//
//	//invalidate used geodesic
//	current_candidates.geodesics_merit(best_index) = -1;
//	//geodesics_candidates.used_geodesics.push_back(best_index);
//
//	return true;
//}
//
//void FeatureGeodesicsController::update_geodesics_merit(int index)
//{
//	if (!is_initialzed)
//		return;
//	if (data_model.feature_geodesics_candidates.size() < 1)
//		return;
//
//
//	compute_dynamic_features(index);
//	write_log(3) << endl << "  update geodesics merit..." << endl;
//
//	//auto& current_candidates = data_model.feature_geodesics_candidates[data_model.current_candidates_index];
//	auto& current_candidates = data_model.feature_geodesics_candidates[index];
//	const int number_geodesics = current_candidates.paths.size();
//	double max_merit = 0.0;
//
//	//update merit of geodesics
//	for (int i = 0; i < number_geodesics; i++)
//	{
//		if (current_candidates.geodesics_merit(i) == -1)
//			continue;
//
//		double merit = 0;
//
//		for (int j = 0; j < data_model.feature_settings.number_geodesic_features; j++)
//		{
//			double weight = data_model.feature_settings.geodesic_feature_weights[j];
//			double value = current_candidates.geodesic_features(i, j);
//
//			if (weight < 0)
//			{
//				value = 1 - value;
//				weight *= -1;
//			}
//
//			merit += weight * value;
//		}
//
//		current_candidates.geodesics_merit(i) = merit;
//		
//		if (merit > max_merit)
//			max_merit = merit;
//	}
//
//	//log_list(5, geodesics_candidates.used_geodesics, "used_geodesics: ");
//	write_log(5) << "geodesic_features: " << endl << current_candidates.geodesic_features << endl;
//	write_log(5) << "geodesics_merit: " << endl << current_candidates.geodesics_merit << endl;
//	write_log(4) << "max_merit: " << max_merit << endl;
//
//	write_log(3) << "  ...done updating merit." << endl;
//}
//
//void FeatureGeodesicsController::compute_dynamic_features(int index)
//{
//	//if (data_model.current_candidates_index >= data_model.feature_geodesics_candidates.size())
//	if (index >= data_model.feature_geodesics_candidates.size())
//		return;
//
//	write_log(3) << endl << "  compute dynamic features of all geodesics..." << endl;
//	const int loglevel = 5;
//
//	//auto& current_candidates = data_model.feature_geodesics_candidates[data_model.current_candidates_index];
//	auto& current_candidates = data_model.feature_geodesics_candidates[index];
//
//	if (!data_model.target.uncovered_V.rows())
//		return;
//
//	if (!data_model.concatenated_patches_V.rows())
//	{
//		current_candidates.geodesic_features.col(data_model.feature_settings.dynamic_geodesic_features_index).setOnes();
//		write_log(loglevel) << "geodesic_dynamic_features : " << endl << current_candidates.geodesic_features << endl;
//		return;
//	}
//
//	int current_feature = data_model.feature_settings.dynamic_geodesic_features_index;
//
//	//TODO performance: don't re-evaluate covered geodesics! only search on closest patch instead of concatenated?
//	double coverage_threshold = data_model.optimization_settings.geodesics_coverage_threshold;
//	double max_distance_absolute = 0.0;
//	double max_distance_realtive = 0.0;
//
//	for (int i = 0; i < current_candidates.paths.size(); i++)
//	{
//		if (current_candidates.geodesics_merit(i) == -1)
//		{
//			current_candidates.geodesic_features(i, current_feature) = 0.0;
//			continue;
//		}
//
//		if(is_equal(current_candidates.geodesic_features(i, current_feature), -1.0))
//			continue;
//
//		const auto& path = current_candidates.paths[i];
//
//		//TODO refactor
//		Eigen::VectorXi covered_indices, uncovered_indices, face_indices;
//		Eigen::MatrixXd closest_points;
//		Eigen::VectorXd distances_squared;
//		igl::point_mesh_squared_distance(path, data_model.concatenated_patches_V, data_model.concatenated_patches_F, distances_squared, face_indices, closest_points);
//
//		double distance_absolute = distances_squared.sum();
//		double distance_relative = distance_absolute / current_candidates.lengths(i);
//
//		if (distance_relative > max_distance_realtive)
//			max_distance_realtive = distance_relative;
//		if (distance_absolute > max_distance_absolute)
//			max_distance_absolute = distance_absolute;
//
//		//current_candidates.geodesic_features(i, current_feature) = distance_absolute <= coverage_threshold ? 0.0 : distance_relative;
//		current_candidates.geodesic_features(i, current_feature) = distance_relative;
//		if (distance_absolute <= coverage_threshold)
//			current_candidates.geodesic_features(i, current_feature) = -1.0;
//
//		write_log(loglevel) << "geodesic_features[" << i << "][" << current_feature << "]: " << distance_relative << "  ("<< distance_absolute << " / " << current_candidates.lengths(i) << ")"  << endl;
//	}
//	write_log(4) << "max_distance_absolute: " << max_distance_absolute << endl;
//	write_log(4) << "max_distance_realtive: " << max_distance_realtive << endl;
//
//	if (max_distance_realtive < coverage_threshold)
//		is_feature_covered = true;
//
//	//if (abs(max_distance_realtive) > 1e-6) //avoid division by zero
//	if (abs(max_distance_realtive) > 1)
//		current_candidates.geodesic_features.col(current_feature) /= max_distance_realtive;
//
//	//write_log(loglevel) << "geodesic_dynamic_features " << current_feature << ": " << endl << current_candidates.geodesic_features << endl;
//
//	write_log(3) << "  ...done computing dynamic features." << endl;
//}
//
//void FeatureGeodesicsController::compute_constant_features(const GeodesicFeatureConfig& feature_settings, GeodesicCandidates& geodesic_candidates)
//{
//	const int loglevel = 5;
//	int current_feature_index = 0;
//
//	const int number_geodesics = geodesic_candidates.paths.size();
//	geodesic_candidates.geodesic_features.resize(number_geodesics, feature_settings.number_geodesic_features);
//
//	add_normalized_constant_feature(current_feature_index++, geodesic_candidates.gauss_curvatures, geodesic_candidates.geodesic_features);
//	add_normalized_constant_feature(current_feature_index++, geodesic_candidates.lengths, geodesic_candidates.geodesic_features);
//	add_normalized_constant_feature(current_feature_index++, geodesic_candidates.curvatures, geodesic_candidates.geodesic_features);
//
//	write_log(loglevel) << "geodesic_constant_features " << endl << geodesic_candidates.geodesic_features << endl;
//}
//
//void FeatureGeodesicsController::add_normalized_constant_feature(const int current_feature_index, const Eigen::VectorXd& feature, Eigen::MatrixXd& out_features)
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
//void FeatureGeodesicsController::find_geodesics_on_target(const TargetMesh& target, std::vector<int>& sources, std::vector<int>& targets, std::vector<Eigen::MatrixXd>& out_geodesic_paths, Eigen::VectorXd& out_geodesic_lengths)
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
//GeodesicCandidates FeatureGeodesicsController::get_all_geodesics(DataModel& data_model, std::vector<int> sources, std::vector<std::vector<int>> targets)
//{
//	sources.pop_back(); //remove the last entry, which is -1 and signifies the end of the cluster
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
//	GeodesicCandidates candidates;
//	find_geodesics_on_target(data_model.target, sources, concatenated_targets, candidates.paths, candidates.lengths);
//
//	candidates.number_geodesics = candidates.paths.size();
//	candidates.geodesics_merit.resize(candidates.number_geodesics);
//	candidates.geodesic_features.resize(candidates.number_geodesics, data_model.feature_settings.number_geodesic_features);
//
//	if (candidates.number_geodesics < 1)
//		return candidates;
//
//	preprocess_target_geodesics(data_model.target, candidates);
//	compute_constant_features(data_model.feature_settings, candidates);
//
//	return candidates;
//}
//
//void FeatureGeodesicsController::preprocess_target_geodesics(const TargetMesh& target, GeodesicCandidates& geodesic_candidates)
//{
//	const int loglevel = 5;
//	write_log(loglevel - 1) << endl << "preprocessing geodesics..." << endl;
//
//	//coverage.resize(number_geodesics);
//	//coverage.setConstant(100000);
//
//	int number_paths = geodesic_candidates.paths.size();
//
//	geodesic_candidates.curvatures.resize(number_paths);
//	geodesic_candidates.gauss_curvatures.resize(number_paths);
//	
//	for (int i = 0; i < number_paths; i++)
//	{
//		geodesic_candidates.paths[i] = curve::remove_straight_segments(geodesic_candidates.paths[i]);
//		geodesic_candidates.curvatures(i) = curve::compute_curvature_sum(geodesic_candidates.paths[i]);
//		geodesic_candidates.gauss_curvatures_lookup(i) = curve::compute_gauss_curvature_sum(target, geodesic_candidates.paths[i]);
//	}
//
//	create_lookup(geodesic_candidates.gauss_curvatures, geodesic_candidates.gauss_curvatures_lookup);
//	create_lookup(geodesic_candidates.curvatures, geodesic_candidates.curvatures_lookup);
//	create_lookup(geodesic_candidates.lengths, geodesic_candidates.lengths_lookup);
//
//	write_log(loglevel - 1) << "...done preprocessing geodesics" << endl << endl;
//}
//
//
//void FeatureGeodesicsController::resample_wrapper_geodesic(const Eigen::MatrixXd& to_sample, const int center_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
//{
//	Eigen::MatrixXd sampled_first_half;
//	double sampled_length_first_half;
//	sample_backward(to_sample, center_index, direction*(-1), sampled_first_half, sampled_length_first_half);
//	write_log(6) << "sampled_first_half: " << endl << sampled_first_half << endl;
//
//	Eigen::MatrixXd sampled_second_half;
//	double sampled_length_second_half;
//	sample_forward(to_sample, center_index, direction, sampled_second_half, sampled_length_second_half);
//	write_log(6) << "sampled_second_half: " << endl << sampled_second_half << endl;
//
//	out_sampled_length = sampled_length_first_half + sampled_length_second_half;
//
//	out_sampled.resize(to_sample.rows(), 3);
//	out_sampled <<
//		sampled_first_half,
//		Eigen::RowVector3d::Zero(),
//		sampled_second_half;
//
//
//	if (LOG_LEVEL == 6) //CHECK
//	{
//		write_log(0) << endl << endl;
//
//		double to_sample_length = 0;
//		double sampled_length = 0;
//		for (int i = 1; i < to_sample.rows(); i++)
//		{
//			double length_original = (to_sample.row(i) - to_sample.row(i - 1)).norm();
//			double length_resampled = (out_sampled.row(i) - out_sampled.row(i - 1)).norm();
//
//			to_sample_length += length_original;
//			sampled_length += length_resampled;
//
//			write_log(0) << "length_original:  " << std::to_string(length_original) << endl;
//			write_log(0) << "length_resampled: " << std::to_string(length_resampled) << endl << endl;
//
//			if (abs(length_original - length_resampled) > 1e-3)
//				write_log(0) << "   it's off at [" << i << "] by " << std::to_string(length_original - length_resampled) << endl << endl;
//		}
//		write_log(0) << "TOTAL length_original:  " << std::to_string(to_sample_length) << endl;
//		write_log(0) << "TOTAL length_resampled: " << std::to_string(sampled_length) << endl << endl;
//		write_log(0) << endl << endl;
//	}
//
//	write_log(6) << endl;
//	write_log(4) << "wrapper sampled_length: " << out_sampled_length << endl;
//	write_log(6) << "sampled_length_first_half: " << sampled_length_first_half << endl;
//	write_log(6) << "sampled_length_second_half: " << sampled_length_second_half << endl;
//	write_log(6) << "out_sampled: " << endl << out_sampled << endl;
//	write_log(6) << endl;
//}
//
//void FeatureGeodesicsController::sample_forward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
//{
//	out_sampled.resize(to_sample.rows() - index - 1, 3);
//	out_sampled_length = 0;
//
//	int out_i = 0;
//	for (int i = index + 1; i < to_sample.rows(); i++)
//	{
//		double segment_length = (to_sample.row(i - 1) - to_sample.row(i)).norm();
//		out_sampled_length += segment_length;
//		out_sampled.row(out_i) = direction * out_sampled_length;
//
//		write_log(6) << "   geodesic: wrapper resampled[" << out_i << "] = " << out_sampled.row(out_i) << "   (segment_l=" << segment_length << ", summed_l=" << out_sampled_length << ")" << endl;
//
//		out_i++;
//	}
//}
//
//void FeatureGeodesicsController::sample_backward(const Eigen::MatrixXd& to_sample, const int index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
//{
//	out_sampled.resize(index, 3);
//	out_sampled_length = 0;
//
//	for (int i = index - 1; i >= 0; i--)
//	{
//		double segment_length = (to_sample.row(i) - to_sample.row(i + 1)).norm();
//		out_sampled_length += segment_length;
//		out_sampled.row(i) = direction * out_sampled_length;
//
//		write_log(6) << "   geodesic: wrapper resampled[" << i << "] = " << out_sampled.row(i) << "   (segment_l=" << segment_length << ", summed_l=" << out_sampled_length << ")" << endl;
//	}
//}
//
//void FeatureGeodesicsController::resample_target_geodesic(int index, Eigen::MatrixXd& out_sampled, double& out_sampled_length)
//{
//	const int loglevel = 6;
//
//	Eigen::MatrixXd original_geodesic = data_model.feature_geodesics_candidates[data_model.current_candidates_index].paths[index];
//	const int number_points = original_geodesic.rows();
//	write_log(loglevel) << "--> original_geodesic: " << endl << original_geodesic << endl;
//
//	if (!original_geodesic.allFinite())
//	{
//		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesic is not finite!" << endl << endl;
//		write_log(1) << "--> original_geodesic: " << endl << original_geodesic << endl << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	out_sampled.resize(number_points, 3);
//	out_sampled_length = 0;
//
//	Eigen::RowVectorXd current_point = original_geodesic.row(0);
//	out_sampled.row(0) = current_point;
//	write_log(loglevel) << "  adding point at " << 0 << ": " << current_point << endl << endl;
//
//	const int edge_length = 1.0;
//	int segment_start_index = 0;
//	int segment_end_index = 0;
//	double current_length = 0.0;
//	int number_points_added = 1;
//
//	//while (segment_end_index < number_points - 1)
//	while (segment_end_index < number_points)
//	//while (true)
//	{
//		if (number_points_added >= out_sampled.rows() - 1)
//			out_sampled.conservativeResize(number_points_added + number_points, Eigen::NoChange);
//
//		Eigen::RowVectorXd next_point;
//		current_length = 0.0;
//		segment_end_index = segment_start_index;
//
//		while (current_length < edge_length && segment_end_index < number_points - 1)
//		{
//			segment_end_index++;
//			next_point = original_geodesic.row(segment_end_index);
//			current_length = (next_point - current_point).norm(); //could be squared norm
//		}
//
//		if (segment_end_index - segment_start_index == 1 && current_length > edge_length)
//		{
//			write_log(loglevel) << "in ONE segment: " << current_length << endl;
//			write_log(loglevel) << "  segment (" << segment_start_index << " - " << segment_end_index << ")" << endl;
//			Eigen::RowVectorXd segment_direction = (original_geodesic.row(segment_end_index) - original_geodesic.row(segment_start_index)).normalized();
//			Eigen::RowVectorXd point = current_point + segment_direction;
//
//			double check_length = (point - current_point).norm();
//			//assert(check_length < 1e-3);
//			write_log(loglevel) << "  -> check_length: " << check_length << " (" << (point - out_sampled.row(number_points_added - 1)).norm() << ")" << endl;
//			write_log(loglevel) << "  adding point at " << number_points_added << ": " << point << endl << endl;
//
//			out_sampled_length += check_length;
//			out_sampled.row(number_points_added) = point;
//			number_points_added++;
//
//			current_point = point;
//			continue;
//		}
//
//		//if is last point to add
//		if (segment_end_index >= number_points - 1 && current_length < edge_length)
//		{
//			write_log(loglevel) << "in LAST segment: " << current_length << endl;
//			write_log(loglevel) << "  remaining current_length: " << current_length << " (end_index: " << segment_end_index << ")" << endl;
//			write_log(loglevel) << "  adding point at " << number_points_added << ": " << original_geodesic.row(segment_end_index) << endl << endl;
//
//			out_sampled_length += current_length;
//			out_sampled.row(number_points_added) = original_geodesic.row(segment_end_index);
//			number_points_added++;
//
//			break;
//		}
//
//		write_log(loglevel) << "in MULTIPLE segments: " << current_length << endl;
//		write_log(loglevel) << "  segment (" << segment_start_index << " - " << segment_end_index << ")" << endl;
//
//		segment_start_index = segment_end_index - 1;
//
//		double a = (original_geodesic.row(segment_start_index) - current_point).norm();
//		double b = (original_geodesic.row(segment_end_index) - original_geodesic.row(segment_start_index)).norm();
//		double c = (original_geodesic.row(segment_end_index) - current_point).norm();
//		write_log(loglevel) << "  length span: a = " << a << ", b = " << b << ", c = " << c << endl;
//
//		double proportional_b = (1 - a) / (c - a);
//		Eigen::RowVectorXd segment_direction = (original_geodesic.row(segment_end_index) - original_geodesic.row(segment_start_index)).normalized();
//		Eigen::RowVectorXd point = original_geodesic.row(segment_start_index) + segment_direction * (b * proportional_b);
//		//auto point = current_point + segment_direction * (b * proportional_b);
//		write_log(loglevel) << "  proportional_b: " << proportional_b << endl;
//
//		double check_length = (point - current_point).norm();
//		//assert(check_length < 1e-3);
//		write_log(loglevel) << "  -> check_length: " << check_length << " (" << (point - out_sampled.row(number_points_added - 1)).norm() << ")" << endl;
//		write_log(loglevel) << "  adding point at " << number_points_added << ": " << point << endl << endl;
//
//		out_sampled_length += check_length;
//		out_sampled.row(number_points_added) = point;
//		number_points_added++;
//
//		current_point = point;
//	}
//
//	out_sampled.conservativeResize(number_points_added, Eigen::NoChange);
//	write_log(loglevel) << "--> out_sampled: " << endl << out_sampled << endl;
//	write_log(loglevel) << "--> target resampled_length: " << out_sampled_length << ", number points: " << number_points_added << endl << endl;
//
//
//	//double geodesic_length = geodesics_candidates.geodesic_lengths(index);
//
//	////double-check resampling
//	//double off = geodesic_length - out_sampled_length;
//	//assert(abs(off) > 5e-2);
//	//write_log(2) << endl << endl << "ERROR: resampled geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(off) << endl << endl;
//
//	//geodesic_length = resampled_length;
//	//target_geodesic = resampled_geodesic;
//	//data_model.target_curve = target_geodesic;
//
//	//debug_show_curve(original_geodesic, data_model, true);
//}
//
//void scale_target(DataModel& data_model)
//{
//	//double max_length = 0.0;
//	//for (auto& geodesics_candidates : data_model.feature_geodesics_candidates)
//	//{
//	//	if (geodesics_candidates.lengths.rows() < 1)
//	//		continue;
//
//	//	double current_max_length = geodesics_candidates.lengths(geodesics_candidates.lengths_lookup(0));
//	//	if (current_max_length > max_length)
//	//		max_length = current_max_length;
//	//	
//	//	write_log(5) << "geodesic: cluster max length: " << current_max_length << endl;
//	//}
//
//	//double wrapper_length = (data_model.wrapper_scale_length - 0.0001);
//	//double scale = wrapper_length / max_length;
//	//write_log(4) << "geodesic: target scale: " << scale << " (max_length: " << max_length << ", defined wrapper_length: " << wrapper_length << ")" << endl << endl;
//
//	//if (abs(scale - 1.0) < 1e-3)
//	//	return;
//
//	//for (int i = 0; i < data_model.feature_geodesics_candidates.size(); i++)
//	//{
//	//	auto& geodesic_candidates = data_model.feature_geodesics_candidates[i];
//	//	geodesic_candidates.lengths *= scale;
//
//	//	for (int i = 0; i < geodesic_candidates.paths.size(); i++)
//	//		geodesic_candidates.paths[i] *= scale;
//	//}
//
//	//data_model.target.V *= scale;
//	//data_model.target.scale = scale;
//	////recompute normals if stored in model
//}
//
//
//bool is_probe_intersecting_features(const DataModel& data_model, const vector<int>& all_close_vertices, const int feature_index, const int candidate_feature_index)
//{
//	//each vertex index can only be part of one feature!!
//
//	//set-intersect with all other clusters, if set-intersection check with igl::segment_segment_intersect
//
//	auto feature = data_model.clustered_feature_vertices[feature_index];
//	feature.pop_back();
//	std::sort(feature.begin(), feature.end());
//
//	auto candidate_feature = data_model.clustered_feature_vertices[candidate_feature_index];
//	candidate_feature.pop_back();
//	std::sort(candidate_feature.begin(), candidate_feature.end());
//
//	for (int c = 0; c < data_model.clustered_feature_vertices.size(); c++)
//	{
//		if (c == feature_index || c == candidate_feature_index)
//			continue;
//
//		auto tested_cluster = data_model.clustered_feature_vertices[c];
//
//		std::vector<int> cluster_copy(tested_cluster.size() - 1);
//		std::copy(tested_cluster.begin(), tested_cluster.end() - 1, cluster_copy.begin());
//		std::sort(cluster_copy.begin(), cluster_copy.end());
//
//		std::vector<int> intersection(all_close_vertices.size() + cluster_copy.size());
//		auto iterator = std::set_intersection(all_close_vertices.begin(), all_close_vertices.end(), cluster_copy.begin(), cluster_copy.end(), intersection.begin());
//		intersection.resize(iterator - intersection.begin());
//
//		bool is_begin = iterator == intersection.begin();
//
//		if (intersection.size() > 0)
//		{
//			//if an intersection was found, check that it is not with the 2 features we are actually testing (happens because vertices are not unique per feature)
//
//			std::vector<int> intersection_feature(feature.size() + intersection.size());
//			auto iterator_feature = std::set_intersection(feature.begin(), feature.end(), intersection.begin(), intersection.end(), intersection_feature.begin());
//			if (iterator_feature != intersection_feature.begin())
//				continue;
//
//			std::vector<int> intersection_candidate_feature(candidate_feature.size() + intersection.size());
//			auto iterator_candidate_feature = std::set_intersection(candidate_feature.begin(), candidate_feature.end(), intersection.begin(), intersection.end(), intersection_candidate_feature.begin());
//			if (iterator_candidate_feature != intersection_candidate_feature.begin())
//				continue;
//
//			//intersection was found that is not with the original features, so return
//			return true;
//		}
//	}
//
//	return false;
//}
//
//void get_closest_vertices(const Eigen::MatrixXd& V, const Eigen::MatrixXd& path, std::vector<int>& out_closest_vertices)
//{
//	for (int i = 0; i < path.rows(); i++)
//	{
//		int index;
//		(V.rowwise() - path.row(i)).rowwise().squaredNorm().minCoeff(&index);
//		out_closest_vertices.push_back(index);
//	}
//}
