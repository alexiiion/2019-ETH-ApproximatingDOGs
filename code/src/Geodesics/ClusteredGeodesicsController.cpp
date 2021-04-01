//#include "ClusteredGeodesicsController.h"
//
//#include <igl/point_mesh_squared_distance.h>
//#include <igl/Timer.h>
//
//#include "PatchModel.h"
//#include "MeshController.h"
//#include "CurveHelper.h"
//
//#include "GeodesicWalker.h"
//#include "Rulings/RuledDevelopableSurface.h"
////#include "Spline.h"
//
//#include "Utils.h"
//#include "Logger.h"
//
//using namespace std;
//
//const Eigen::RowVector3d sample_direction = Eigen::RowVector3d::UnitX();
//
////std::vector<Mesh> ruled_developable_surfaces;
//
//
//double compute_standard_deviation(const std::vector<double>& values, const double& mean)
//{
//	double variance = 0.0;
//	for (int i = 0; i < values.size(); i++)
//		variance += (values[i] - mean) * (values[i] - mean);
//
//	variance *= 1.0 / values.size();
//	return sqrt(variance);
//}
//
//Eigen::VectorXd ortho_project2(const Eigen::VectorXd projection_target, const Eigen::VectorXd to_project)
//{
//	//return to_project - (projection_target.dot(to_project) / projection_target.dot(projection_target) * projection_target);
//	return (to_project.dot(projection_target) * projection_target) / projection_target.squaredNorm();
//}
//void print_mathematica_vector(const Eigen::VectorXd vector)
//{
//	//cout << endl << endl;
//	std::cout << "{";
//
//	for (int i = 0; i < vector.rows(); i++)
//	{
//		std::cout << vector(i);
//		if (i < vector.rows() - 1)
//			std::cout << ", ";
//	}
//
//	std::cout << "}" << std::endl;
//	//cout << endl << endl;
//}
//
//void print_mathematica_vector2d(const std::vector<std::vector<double>>& list)
//{
//	//cout << endl << endl;
//	std::cout << "{";
//	for (int r = 0; r < list.size(); r++)
//	{
//		std::cout << "{";
//		for (int c = 0; c < list[r].size(); c++)
//		{
//			std::cout << list[r][c];
//			if (c < list[r].size() - 1)
//				std::cout << ", ";
//		}
//		if (r < list.size() - 1)
//			std::cout << "}," << std::endl;
//	}
//	std::cout << "}}" << std::endl;
//	//cout << endl << endl;
//}
//void print_mathematica_vector2d(const std::vector<std::vector<int>>& list)
//{
//	//cout << endl << endl;
//	std::cout << "{";
//	for (int r = 0; r < list.size(); r++)
//	{
//		std::cout << "{";
//		for (int c = 0; c < list[r].size(); c++)
//		{
//			std::cout << list[r][c];
//			if (c < list[r].size() - 1)
//				std::cout << ", ";
//		}
//		if (r < list.size() - 1)
//			std::cout << "}," << std::endl;
//	}
//	std::cout << "}}" << std::endl;
//	//cout << endl << endl;
//}
//
//void print_measure(const std::vector<std::vector<double>>& list, const DataModel& data_model)
//{
//	const int loglevel = 4;
//	for (int i = 0; i < list.size(); i++)
//	{
//		auto neighbors = list[i];
//		write_log(loglevel) << "[" << i << "]" << endl;
//
//		for (int j = 0; j < neighbors.size(); j++)
//			write_log(loglevel) << "  [" << data_model.target.adjacency_VV[i][j] << "]: " << list[i][j] << endl;
//	}
//}
//
//std::vector<int> ClusteredGeodesicsController::get_random_geodesics(int n)
//{
//	std::vector<int> random_vertices(n);
//
//	//get n random vertices (later: farthest point sampling)
//	std::srand(unsigned(std::time(0)));
//
//	const int n_vertices = data_model.target.V.rows();
//	int i = 0;
//
//	while (i < n)
//	{
//		int random_index = std::rand() % n_vertices;
//
//		Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[random_index];
//		if (geodesic.rows() < 10)
//			continue;
//
//		auto iterator = std::find(random_vertices.begin(), random_vertices.end(), random_index);
//		bool is_found = iterator != random_vertices.end();
//
//		if (is_found)
//			continue;
//
//		random_vertices[i] = random_index;
//		i++;
//	}
//
//	log_list(4, random_vertices, "random indices for ruled developables: ", false);
//	return random_vertices;
//}
//
//void ClusteredGeodesicsController::get_random_ruled_developables(int n)
//{
//	random_geodesics = get_random_geodesics(n);
//	
//	for (int i = 0; i < n; i++)
//	{
//		Eigen::MatrixXd geodesic = data_model.geodesics_candidates.paths[random_geodesics[i]];
//
//		if (geodesic.rows() < 1)
//			continue;
//
//		//RuledDevelopableSurface surface;
//		//Mesh surface_mesh = surface.create(geodesic, 15);
//		
//		/*
//		* used the following code to compute rulings from splines
//		* removed now to work with discrete curves
//		*/
//		/*
//		//create spline
//		Spline spline(geodesic, smoothness, 5);
//		//const double sample_step = 0.5;
//		Eigen::MatrixXd sampled_curve = spline.sample_curve(sample_step);
//
//
//		//compute rulings (later: check for d''~0 and interpolate, check for singularities)
//		Eigen::MatrixXd derivative_2 = spline.sample_derivative(2, sample_step);
//		Eigen::MatrixXd derivative_3 = spline.sample_derivative(3, sample_step);
//
//		Eigen::MatrixXd rulings(derivative_2.rows(), 3);
//		for (int i = 0; i < derivative_2.rows(); i++)
//		{
//			Eigen::Vector3d order2 = derivative_2.row(i);
//			Eigen::Vector3d order3 = derivative_3.row(i);
//
//			Eigen::Vector3d ruling = order2.cross(order3);
//			rulings.row(i) = ruling.normalized();
//
//			//write_log(4) << ruling << endl;
//		}
//
//		if (LOG_LEVEL < 4)
//		{
//			std::cout << "geodesic="; print_mathematica_matrix(geodesic);
//			std::cout << "spline="; print_mathematica_matrix(sampled_curve);
//			std::cout << "spline2Derivative="; print_mathematica_matrix(derivative_2);
//			std::cout << "spline3Derivative="; print_mathematica_matrix(derivative_3);
//			std::cout << "rulings="; print_mathematica_matrix(rulings);
//			std::cout << "n=" << sampled_curve.rows() << std::endl;
//		}
//		*/
//
//		//meshhelper::add_mesh(data_model.viewer, surface_mesh.V, surface_mesh.F, false, false);
//	}
//
//
//	//select using graph cut
//}
//
//bool ClusteredGeodesicsController::initialize()
//{
//	if (is_initialzed)
//		return true;
//
//	if (data_model.geodesics_candidates.is_empty())
//	{
//		//build_geodesics();
//		//compute_geodesics_descriptors();
//		//cluster_geodesics();
//	}
//	//else
//	//	update_geodesics_merit(data_model.current_candidates_index);
//
//	//if (data_model.feature_geodesics_candidates.size() == 0)
//	//	return false;
//
//	////if (is_equal(data_model.target.scale, 1.0))
//	////	scale_target(data_model);
//
//
//	//random_geodesics = get_random_geodesics(100);
//
//
//
//	//select using graph cut
//
//
//	is_initialzed = true;
//	return is_initialzed;
//}
//
//void ClusteredGeodesicsController::build_geodesics()
//{
//	if (!data_model.geodesics_candidates.is_empty())
//		data_model.geodesics_candidates.clear();
//
//	igl::Timer timer;
//	double init_time = timer.getElapsedTime();
//
//	write_log(4) << endl << "build all walking geodesics..." << endl;
//
//	const int number_vertices = data_model.target.V.rows();
//	GeodesicCandidates candidates(number_vertices);
//	GeodesicWalker walker(data_model.target, &data_model.max_frame_error, &data_model.use_weighted_frame_error);
//
//	for (int i = 0; i < number_vertices; i++)
//	{
//		Eigen::MatrixXd geodesic;
//		double length;
//		
//		if(data_model.geodesics_direction == GeodesicsDirection::K_MIN)
//			walker.get_geodesic_kmin_at(i, geodesic, length);
//		else if(data_model.geodesics_direction == GeodesicsDirection::K_MAX)
//			walker.get_geodesic_kmax_at(i, geodesic, length);
//
//		//walker.get_geodesic_at(i, geodesic, length);
//		
//		//Eigen::RowVector3d k_max = data_model.target.surface_features.principal_k_max_direction.row(i);
//		//get_geodesic_at(i, geodesic, length);
//		
//		//Eigen::MatrixXd geodesic = get_geodesic_at(i);
//		//if (geodesic.rows() < 5)
//		//	continue;
//
//		////this is really slow! (39s instead of 1s on the bunny's 3301 vertices)
//		//double length = curve::compute_length(geodesic);
//		//double curvature = curve::compute_curvature(geodesic);
//		//double gauss_curvature = curve::compute_gauss_curvature(data_model.target, geodesic);
//
//		candidates.paths.push_back(geodesic);
//		candidates.lengths(i) = length;
//		//candidates.curvatures(i) = curvature;
//		//candidates.gauss_curvatures(i) = gauss_curvature;
//
//		candidates.number_geodesics++;
//	}
//
//	//candidates.lengths.conservativeResize(candidates.number_geodesics);
//	//candidates.curvatures.conservativeResize(candidates.number_geodesics);
//	//candidates.gauss_curvatures.conservativeResize(candidates.number_geodesics);
//
//	//create_lookup(candidates.lengths, candidates.lengths_lookup);
//	//create_lookup(candidates.curvatures, candidates.curvatures_lookup);
//	//create_lookup(candidates.gauss_curvatures, candidates.gauss_curvatures_lookup);
//
//	data_model.geodesics_candidates = candidates;
//	data_model.do_rebuild_geodesics = false;
//	//get_rulings(path, target.normals_vertices.row(source_vertex));
//
//	double t = timer.getElapsedTime();
//	write_log(4) << "...done building geodesics. found " << candidates.paths.size() << " geodesics -- ELAPSED TIME: " << t - init_time << endl << endl;
//
//	//write_log(4) << "geodesic lengths: " << endl << candidates.lengths << endl << endl;
//}
//
//int get_similar_neighbor(const GeodesicCandidates& geodesics_candidates, const std::vector<bool>& visited)
//{
//	int index = -1;
//	double min_similarity = 1e6;
//
//	for (int i = 0; i < geodesics_candidates.pairwise_similarity.size(); i++)
//	{
//		if (visited[i])
//			continue;
//
//		auto element = *std::min_element(geodesics_candidates.pairwise_similarity[i].begin(), geodesics_candidates.pairwise_similarity[i].end());
//		if (element > min_similarity)
//			continue;
//
//		min_similarity = element;
//		index = i;
//
//		write_log(5) << i << ": min_similarity = " << min_similarity << endl;
//	}
//
//	return index;
//}
//
//void ClusteredGeodesicsController::get_cluster_at(const int& start_index, std::vector<bool>& visited, std::vector<int>& out_cluster_vertices)
//{
//	const int log_level = 5;
//
//	out_cluster_vertices.clear();
//	GeodesicCandidates geodesics_candidates = data_model.geodesics_candidates;
//
//	//get median & iqr for first bounds
//	auto start_neighbors = data_model.geodesics_candidates.pairwise_similarity[start_index];
//	int num = start_neighbors.size();
//	std::sort(start_neighbors.begin(), start_neighbors.end());
//
//	int n = num % 2 == 0 ? num / 2 : (num - 1) / 2;
//	double median = num % 2 == 0 ? (start_neighbors[n-1] + start_neighbors[n]) / 2 : start_neighbors[n];
//
//	int nq = n % 2 == 0 ? n / 2 : (n - 1) / 2;
//	double q1 = n % 2 == 0 ? (start_neighbors[nq-1] + start_neighbors[nq]) / 2 : start_neighbors[nq];
//	double q3 = n % 2 == 0 ? (start_neighbors[n+nq-1] + start_neighbors[n+nq]) / 2 : start_neighbors[num-1-nq];
//	
//	// keep track of descriptive statistics for cluster adaption
//	double similarites_sum = 0.0;
//	double similarites_sum_squared = 0.0;
//
//	double min_similarity = *std::min_element(geodesics_candidates.pairwise_similarity[start_index].begin(), geodesics_candidates.pairwise_similarity[start_index].end());
//	double mean = min_similarity;
//	double bounds = median/2;
//	double stdev = 0.0;
//	//double mean = median;
//	//double bounds = (q3-q1)/2;
//	bounds = bounds < data_model.bounds_min_abs ? data_model.bounds_min_abs : bounds;
//
//	std::stack<int> open;
//	open.push(start_index);
//
//	out_cluster_vertices.push_back(start_index);
//	vector<double> cluster_similarities;
//
//
//	//// Get first neighbor and push into cluster
//	////double min_similarity = *std::min_element(geodesics_candidates.pairwise_similarity[start_index].begin(), geodesics_candidates.pairwise_similarity[start_index].end());
//	//auto min_iterator = std::min_element(geodesics_candidates.pairwise_similarity[start_index].begin(), geodesics_candidates.pairwise_similarity[start_index].end());
//	//int min_index = std::distance(geodesics_candidates.pairwise_similarity[start_index].begin(), min_iterator);
//	//double min_similarity = *min_iterator;
//
//	//int first_clustered_vertex = data_model.target.adjacency_VV[start_index][min_index];
//	//visited[first_clustered_vertex] = true;
//
//	//std::stack<int> open;
//	//open.push(start_index);
//	//open.push(first_clustered_vertex);
//	//out_cluster_vertices.push_back(start_index);
//	//out_cluster_vertices.push_back(first_clustered_vertex);
//
//	//vector<double> cluster_similarities;
//	//cluster_similarities.push_back(min_similarity);
//	//
//	//// keep track of descriptive statistics
//	//double similarites_sum = min_similarity;
//	//double similarites_sum_squared = min_similarity*min_similarity;
//	//
//	//double mean = min_similarity;
//	//double stdev = 0.0;
//	//double bounds = mean * data_model.bounds_factor_initial;
//	//bounds = bounds < data_model.bounds_min_abs ? data_model.bounds_min_abs : bounds;
//
//
//	int vertex_index;
//	while (!open.empty())
//	{
//		// Pop a vertex from stack 
//		vertex_index = open.top();
//		open.pop();
//
//		if (visited[vertex_index])
//			continue;
//
//		//TODO use n-ring-neighbors
//		//auto neighbors = data_model.target.adjacency_VV[vertex_index];
//		auto neighbors = get_n_ring_neighborhood(vertex_index, 1, data_model.target.adjacency_VV);
//
//		write_log(log_level) << endl << endl << "current index: " << vertex_index << endl;
//		log_list(log_level, neighbors, "neighbors: ", false);
//		log_list(log_level, data_model.geodesics_candidates.pairwise_similarity[vertex_index], "similarities: ", false);
//
//		for (int i = 0; i < neighbors.size(); i++)
//		{
//			int neighbor_index = neighbors[i];
//			write_log(log_level) << endl << "  neighbor_index: " << neighbor_index;
//
//			if (visited[neighbor_index])
//				continue;
//
//			double similarity = geodesics_candidates.pairwise_similarity[vertex_index][i];
//			write_log(log_level) << ", similarity: " << similarity;
//			if (similarity > 1e6)
//			{
//				visited[neighbor_index] = true;
//				continue;
//			}
//
//			write_log(log_level) << ", bounds: " << mean - bounds << " - " << mean + bounds << " (mu=" << mean << ", sd=" << bounds << ")";
//
//			//if similarity is within margin, add to cluster
//			if (similarity >= mean - bounds && similarity <= mean + bounds)
//			{
//				std::vector<int>::iterator iterator = std::find(out_cluster_vertices.begin(), out_cluster_vertices.end(), neighbor_index);
//				bool is_contained = iterator != out_cluster_vertices.end();
//				//if(is_contained && !visited[neighbor_index])
//				//	write_log(6) << "strange";
//
//				//if (is_contained)
//				//	continue;
//
//				//write_log(log_level) << endl << "    ADD index! ";
//				//out_cluster_vertices.push_back(neighbor_index);
//				//cluster_similarities.push_back(similarity);
//
//				if (!is_contained)
//				{
//					write_log(log_level) << endl << "    ADD index! ";
//					out_cluster_vertices.push_back(neighbor_index);
//					cluster_similarities.push_back(similarity);
//
//					similarites_sum += similarity;
//					similarites_sum_squared += similarity * similarity;
//				}
//
//				if (cluster_similarities.size() > 2)
//				{
//					mean = similarites_sum / cluster_similarities.size();
//
//					//bounds = mean * data_model.bounds_factor;
//					//stdev = compute_standard_deviation(cluster_similarities, mean);
//					stdev = sqrt(similarites_sum_squared/cluster_similarities.size() - mean*mean);
//					if (isnan(stdev))
//						write_log(4) << "NaN";
//
//					bounds = (stdev * data_model.bounds_factor)/2;
//					bounds = bounds < data_model.bounds_min_abs ? data_model.bounds_min_abs : bounds;
//
//					write_log(log_level) << "    updating mean: " << mean << ", bounds: " << bounds << " (stdev=" << stdev << ")";
//				}
//
//				write_log(log_level) << endl;
//				//log_list(log_level, out_cluster_vertices, "    cluster_vertices: ", false);
//				//log_list(log_level, cluster_similarities, "    cluster_similarities: ", false);
//
//				open.push(neighbor_index);
//			}
//		}
//
//		visited[vertex_index] = true;
//	}
//
//	if (cluster_similarities.size() < 2)
//		return;
//
//	//write_log(log_level) << "original   mean: " << mean << ", stdev: " << stdev << endl;
//	
//	//recompute mean
//	if (cluster_similarities.size() <= 2)
//	{
//		double sum = 0.0;
//		double squared_sum = 0.0;
//		for (int j = 0; j < cluster_similarities.size(); j++)
//		{
//			double similarity = cluster_similarities[j];
//
//			sum += similarity;
//			squared_sum += similarity * similarity;
//		}
//
//		mean = sum / cluster_similarities.size();
//		double var = (squared_sum / cluster_similarities.size()) - (mean*mean); 
//		stdev = sqrt(var); 
//
//		/*
//		write_log(0) << "recomputed mean: " << mean_re << ", stdev: " << stdev_re << ", var: " << var_re;
//		write_log(0) << " (size: " << cluster_similarities.size() << ", sum:  " << sum << ", squared_sum: " << squared_sum << ")" << endl;
//
//		if (isnan(mean) || isnan(mean_re) || isnan(stdev) || isnan(stdev_re) || abs(mean_re - mean) > 1e-6 || abs(stdev_re - stdev) > 1e-6)
//			log_list(0, cluster_similarities, "  similarities: ", false);
//
//		write_log(0) << endl;
//		*/
//	}
//
//	data_model.geodesics_clusters.push_back(out_cluster_vertices);
//	data_model.cluster_values.push_back(cluster_similarities);
//	data_model.cluster_mean.push_back(mean);
//	data_model.cluster_stdev.push_back(stdev);
//}
//
//void ClusteredGeodesicsController::cluster_geodesics()
//{
//	if (data_model.geodesics_candidates.is_empty())
//	{
//		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesics candidates not built yet!" << endl << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (data_model.geodesics_candidates.pairwise_similarity.size() < 1)
//	{
//		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesics candidates similarity not computed yet!" << endl << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	igl::Timer timer;
//	double init_time = timer.getElapsedTime();
//
//	write_log(4) << endl << "clustering geodesics..." << endl;
//
//
//	if (data_model.geodesics_clusters.size() > 0)
//		data_model.geodesics_clusters.clear();
//
//	data_model.geodesics_clusters.clear();
//	//vector<int> vertex_cluster_assingments(number_vertices, -1); //indices as #V
//	//int current_cluster_index = 0; //increment in while loop, when new cluster is started
//
//	vector<bool> visited(data_model.target.V.rows(), false);
//	int start_index = get_similar_neighbor(data_model.geodesics_candidates, visited);
//	bool has_unvisited_vertices = true;
//	
//	////for debug only!
//	//start_index = data_model.pre_selected_vertex;
//	//if (start_index < 0)
//	//	return;
//
//
//	while (has_unvisited_vertices)
//	{
//		//write_log(4) << endl << "  starting cluster at v" << start_index << endl;
//		//cluster by dfs, only mark clustered and inf vertices as visited, dont mark vertices that are not similar
//
//		std::vector<int> cluster_vertices;
//		get_cluster_at(start_index, visited, cluster_vertices);
//
//		////log_list(4, cluster_vertices, "cluster at " + to_string(start_index) + ": ", false);
//		//if(cluster_vertices.size() > 1)
//		//	data_model.geodesics_clusters.push_back(cluster_vertices);
//
//		//start next cluster at unvisited vertex with most similar neighbor geodesic
//		start_index = get_similar_neighbor(data_model.geodesics_candidates, visited);
//		if (start_index == -1)
//			has_unvisited_vertices = false;
//	}
//
//	data_model.do_recluster_geodesics = false;
//
//	double t = timer.getElapsedTime();
//	write_log(4) << "...done clustering geodesics. found " << data_model.geodesics_clusters.size() << " cluster. -- ELAPSED TIME: " << t - init_time << endl << endl;
//
//	//write_log(4) << endl << endl << "found clusters: " << endl;
//	//for (int i = 0; i < data_model.geodesics_clusters.size(); i++)
//	//	log_list(4, data_model.geodesics_clusters[i], to_string(i) + ": ", false);
//}
//
//void ClusteredGeodesicsController::compute_geodesics_descriptors()
//{
//	if (data_model.geodesics_candidates.is_empty())
//	{
//		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> geodesics candidates not built yet!" << endl << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	igl::Timer timer;
//	double init_time = timer.getElapsedTime();
//
//	write_log(4) << endl << "computing all geodesics descriptors..." << endl;
//
//	const int number_vertices = data_model.target.V.rows();
//	GeodesicCandidates& geodesics_candidates = data_model.geodesics_candidates;
//
//	std::vector<std::vector<double>> pairwise_variance(number_vertices);
//	//TODO store this somewhere?
//	std::vector<std::vector<double>> pairwise_mean(number_vertices);
//	std::vector<std::vector<double>> pairwise_distance_sum(number_vertices);
//
//	//TODO recomputes each pair right now, make some lookup?
//	for (int vertex_index = 0; vertex_index < number_vertices; vertex_index++)
//	{
//		auto neighbors = data_model.target.adjacency_VV[vertex_index];
//		int number_neighbors = neighbors.size();
//		
//		pairwise_mean[vertex_index] = std::vector<double>(number_neighbors);
//		pairwise_variance[vertex_index] = std::vector<double>(number_neighbors);
//		pairwise_distance_sum[vertex_index] = std::vector<double>(number_neighbors);
//
//		for (int j = 0; j < number_neighbors; j++)
//		{
//			int neighbor_index = neighbors[j];
//			
//			Eigen::RowVector3d edge_offset(0,0,0);
//			if(data_model.use_vertex_offset_similarity)
//				edge_offset = data_model.target.V.row(vertex_index) - data_model.target.V.row(neighbor_index);
//
//
//			//TODO compute on resampled geodesic! otherwise length influences the variance
//			//TODO cache already computed similarities (if neighbor_index < vertex_index)
//			double distance_mean, distance_variance, distance_sum;
//			distances_between_geodesics(geodesics_candidates.paths[vertex_index], geodesics_candidates.paths[neighbor_index], edge_offset, distance_mean, distance_variance, distance_sum);
//
//			//Eigen::VectorXd shortest_distances;
//			//Eigen::MatrixXd distance_vectors;
//			//Eigen::VectorXd distances_1order;
//			//distances_between_geodesics(geodesics_candidates.paths[vertex_index], geodesics_candidates.paths[neighbor_index], edge_offset, distance_vectors, shortest_distances, distances_1order, distance_mean, distance_variance, distance_sum);
//
//			pairwise_variance[vertex_index][j] = distance_variance;
//			pairwise_mean[vertex_index][j] = distance_mean;
//			pairwise_distance_sum[vertex_index][j] = distance_sum;
//
//			if (distance_variance < geodesics_candidates.min_similarity)
//			{
//				geodesics_candidates.min_similarity = distance_variance;
//				geodesics_candidates.min_similarity_vertex = vertex_index;
//			}
//		}
//	}
//
//	double t = timer.getElapsedTime();
//	write_log(4) << "...done computing geodesics descriptors -- ELAPSED TIME: " << t - init_time << endl << endl;
//
//	geodesics_candidates.pairwise_similarity = pairwise_variance;
//
//	//write_log(4) << endl << "mean: " << endl << endl;
//	//print_measure(pairwise_mean, data_model);
//	//write_log(4) << endl << "distance sum: " << endl << endl;
//	//print_measure(pairwise_distance_sum, data_model);
//	//write_log(4) << endl << "variance: " << endl << endl;
//	//print_measure(pairwise_variance, data_model);
//
//	//write_log(4) << endl << "variance (mathematica): " << endl << endl;
//	//print_mathematica_vector2d(pairwise_variance);
//	//write_log(4) << endl << "adjacency VV (mathematica): " << endl << endl;
//	//print_mathematica_vector2d(data_model.target.adjacency_VV);
//
//}
//
//void ClusteredGeodesicsController::distances_between_geodesics(const Eigen::MatrixXd& path, const Eigen::MatrixXd& path_neighbor, const Eigen::RowVector3d& edge_offset, double& out_distance_mean, double& out_distance_variance, double& out_distance_sum)
//{
//	Eigen::VectorXd shortest_distances;
//	Eigen::MatrixXd distance_vectors;
//	Eigen::VectorXd distances_1order;
//	distances_between_geodesics(path, path_neighbor, edge_offset, distance_vectors, shortest_distances, distances_1order, out_distance_mean, out_distance_variance, out_distance_sum);
//}
//
//void ClusteredGeodesicsController::distances_between_geodesics(const Eigen::MatrixXd& path, const Eigen::MatrixXd& path_neighbor, const Eigen::RowVector3d& edge_offset,
//	Eigen::MatrixXd& out_distance_vectors, Eigen::VectorXd& out_shortest_distances, Eigen::VectorXd& out_distances_1order, 
//	double& out_distance_mean, double& out_distance_variance, double& out_distance_sum)
//{
//	if (path.rows() < 1 || path_neighbor.rows() < 1)
//	{
//		out_distance_mean = std::numeric_limits<double>::infinity();
//		out_distance_variance = std::numeric_limits<double>::infinity();
//		out_distance_sum = std::numeric_limits<double>::infinity();
//		return;
//	}
//
//	out_distance_vectors.resize(path.rows(), 3);
//	out_shortest_distances.resize(path.rows());
//	out_distances_1order.resize(path.rows() - 1);
//	//Eigen::VectorXd distances_2order(path.rows() - 2);
//
//	out_distance_sum = 0.0;
//	double distance_squared_sum = 0.0;
//	double distance_sum_1order = 0.0;
//	//double error_sum_2order = 0.0;
//
//	//Eigen::MatrixXd display_distance_vectors(path.rows(), 3);
//	//int last_added_i = -1;
//
//	int number_projected = 0;
//	int last_j = 0;
//
//	for (int i = 0; i < path.rows(); i++)
//	{
//		if (path.rows() < 5 || path_neighbor.rows() < 5) //TODO do better
//			break;
//
//		Eigen::RowVector3d point = path.row(i);
//
//		//for (int j = last_j; j < path_neighbor.rows()-1; j++)
//		for (int j = 0; j < path_neighbor.rows()-1; j++)
//		{
//			Eigen::RowVector3d segment_neighbor = path_neighbor.row(j + 1) - path_neighbor.row(j);
//			
//			Eigen::RowVector3d to_project = point - path_neighbor.row(j);
//			Eigen::RowVector3d projected = ortho_project2(segment_neighbor, to_project);
//
//			double dot = projected.normalized().dot(segment_neighbor.normalized());
//			//negative dot product means the projected vector points opposite of the segement and is thus not within this segment
//			if (dot < 0)
//				continue;
//
//			if (projected.squaredNorm() > segment_neighbor.squaredNorm())
//				continue;
//
//			Eigen::RowVector3d distance_vector = to_project - projected;
//			distance_vector = distance_vector - edge_offset;
//			double distance = distance_vector.norm();
//
//			out_distance_vectors.row(number_projected) = distance_vector;
//			out_shortest_distances(number_projected) = distance;
//			out_distance_sum += distance;
//
//			distance_squared_sum += distance * distance;
//
//
//			//write_log(0) << "dot = " << dot << ", distance = " << distance << endl;
//
//			//write_log(0) << i << ", " << j << ": distance = " << distance;
//			if (number_projected > 0)
//			{
//				out_distances_1order(number_projected - 1) = distance - out_shortest_distances(number_projected - 1);
//				distance_sum_1order += abs(out_distances_1order(number_projected - 1));
//				//write_log(0) << ", 1. order = " << distances_1order(number_projected - 1);
//			}
//			//if (number_projected > 1)
//			//{
//			//	distances_2order(number_projected - 2) = distances_1order(number_projected - 1) - distances_1order(number_projected - 2);
//			//	error_sum_2order += abs(distances_2order(number_projected - 2));
//			//	//write_log(0) << ", 2. order = " << distances_2order(number_projected - 2);
//			//}
//			////write_log(0) << endl;
//
//
//			//display_distance_vectors.row(i) = distance_vector;
//			//last_added_i = i;
//
//			number_projected++;
//			last_j = j;
//			break;
//		}
//
//		//if (last_added_i < i)
//		//{
//		//	display_distance_vectors.row(i) = Eigen::RowVector3d::Zero();
//		//	last_added_i = i;
//		//}
//	}
//
//	//if (last_added_i < path.rows()-1)
//	//	display_distance_vectors.row(path.rows() - 1) = Eigen::RowVector3d::Zero();
//
//	if (number_projected < 1)
//	{
//		out_distance_mean = std::numeric_limits<double>::infinity();
//		out_distance_variance = std::numeric_limits<double>::infinity();
//		out_distance_sum = std::numeric_limits<double>::infinity();
//		
//		out_distance_vectors.conservativeResize(0, Eigen::NoChange);
//		out_shortest_distances.conservativeResize(0);
//		out_distances_1order.conservativeResize(0);
//
//		return;
//	}
//
//	out_distance_vectors.conservativeResize(number_projected, Eigen::NoChange);
//	out_shortest_distances.conservativeResize(number_projected);
//	out_distances_1order.conservativeResize(number_projected - 1);
//	//distances_2order.conservativeResize(number_projected - 2);
//
//	out_distance_mean = out_distance_sum / number_projected;
//	out_distance_variance = (distance_squared_sum / (number_projected-1)) - (out_distance_mean * out_distance_mean);
//
//
//	//write_log(0) << endl << "path = "; print_mathematica_matrix(path);
//	//write_log(0) << endl << "pathNeighbor = "; print_mathematica_matrix(path_neighbor);
//	//write_log(0) << endl << "distanceVectors = "; print_mathematica_matrix(out_distance_vectors);
//	//write_log(0) << endl << "distances = "; print_mathematica_vector(out_shortest_distances);
//	//write_log(0) << endl << "distances1order = "; print_mathematica_vector(distances_1order);
//	//write_log(0) << endl << "displayDistanceVectors = "; print_mathematica_matrix(display_distance_vectors);
//	//write_log(0) << endl << "distances2order = "; print_mathematica_vector(distances_2order);
//}
//
////old: finds shortest distance from all points, even if ends are offset
///*
//void ClusteredGeodesicsController::distances_between_geodesics(const Eigen::MatrixXd& path, const Eigen::MatrixXd& path_neighbor, Eigen::MatrixXd& out_distance_vectors, Eigen::VectorXd& out_shortest_distances)
//{
//	out_distance_vectors.resize(path.rows(), 3);
//	out_shortest_distances.resize(path.rows());
//
//	//double error_sum = 0.0;
//	//double error_1order = 0.0; //should be vector, now it's only a sum
//	//double error_2order = 0.0;
//
//	int last_assigned_i = 0;
//	int last_j = 0;
//
//	for (int i = 0; i < path.rows(); i++)
//	{
//		Eigen::RowVector3d point = path.row(i);
//
//		for (int j = last_j; j < path_neighbor.rows()-1; j++)
//		{
//			Eigen::RowVector3d ef = path_neighbor.row(j + 1) - path_neighbor.row(j);
//			double ef_length = ef.norm();
//
//			Eigen::RowVector3d to_project = point - path_neighbor.row(j);
//			Eigen::RowVector3d projected = ortho_project2(ef, to_project);
//
//			double dot = projected.normalized().dot(ef.normalized());
//			Eigen::RowVector3d distance_vector = dot >= 0 ? (to_project - projected).eval() : to_project;
//			double distance = distance_vector.norm();
//
//			//if dot < 0: shortetst distance is to first point
//			//if projected.norm() < ef_length: then projection is found on this segment
//			if (dot < 0 || projected.norm() < ef_length)
//			{
//				out_distance_vectors.row(i) = distance_vector;
//				out_shortest_distances(i) = distance;
//
//				//write_log(0) << i << ", " << j << ": distance = " << distance;
//				//error_sum += distance;
//				//if (i > 0)
//				//{
//				//	error_1order += abs(distance - out_shortest_distances(i - 1));
//				//	write_log(0) << ", 1. order = " << abs(distance - out_shortest_distances(i - 1));
//				//}
//				//if (i > 1) 
//				//{
//				//	error_2order += abs(abs(out_shortest_distances(i-2) - out_shortest_distances(i-1)) - abs(distance - out_shortest_distances(i - 1)));
//				//	write_log(0) << ", 2. order = " << abs(abs(out_shortest_distances(i - 2) - out_shortest_distances(i - 1)) - abs(distance - out_shortest_distances(i - 1)));
//				//}
//				//write_log(0) << endl;
//
//				last_j = j;
//				last_assigned_i = i;
//				break;
//			}
//		}
//	}
//
//	//if not all points where projected, then path is longer than path_neighbor and all remaining points have shortest distance to last point on path_neighbor
//	if (last_assigned_i < path.rows())
//	{
//		Eigen::RowVector3d point_neighbor = path_neighbor.row(path_neighbor.rows() - 1);
//
//		for (int i = last_assigned_i + 1; i < path.rows(); i++)
//		{
//			Eigen::RowVector3d distance_vector = path.row(i) - point_neighbor;
//			out_distance_vectors.row(i) = distance_vector;
//			out_shortest_distances(i) = distance_vector.norm();
//			//write_log(0) << i << ", " << path_neighbor.rows() - 1 << endl;
//		}
//	}
//
//	//write_log(0) << endl << "path = "; print_mathematica_matrix(path);
//	//write_log(0) << endl << "pathNeighbor = "; print_mathematica_matrix(path_neighbor);
//	//write_log(0) << endl << "distanceVectors = "; print_mathematica_matrix(out_distance_vectors);
//	//write_log(0) << endl << "distances = "; print_mathematica_vector(out_shortest_distances);
//}
//*/
//
//bool ClusteredGeodesicsController::get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
//{
//	if(!is_initialzed)
//		return false;
//
//	if (!data_model.target.uncovered_V.rows())
//		return false;
//
//
//
//
//	if (data_model.patches.size() >= random_geodesics.size())
//		return false;
//
//	select_geodesic(out_target_curve, out_wrapper_curve);
//
//
//	//
//	//write_log(3) << endl << "get next geodesic as target constraints:" << endl;
//	//update_geodesics_merit(data_model.current_candidates_index);
//
//	////if the current feature is covered, search forward
//	//while (is_feature_covered)
//	//{
//	//	write_log(5) << "geodesic_features: " << endl << data_model.feature_geodesics_candidates[data_model.current_candidates_index].geodesic_features << endl;
//	//	write_log(5) << "geodesics_merit: " << endl << data_model.feature_geodesics_candidates[data_model.current_candidates_index].geodesics_merit << endl;
//
//	//	data_model.current_candidates_index++;
//	//	if (data_model.current_candidates_index >= data_model.feature_geodesics_candidates.size())
//	//		return false;
//
//	//	//can't find geodesic in next feature neighborhood because it has no neighbors. skip to next feature
//	//	if (data_model.feature_geodesics_candidates[data_model.current_candidates_index].number_geodesics < 1)
//	//		continue;
//
//	//	is_feature_covered = false;
//	//	update_geodesics_merit(data_model.current_candidates_index);
//	//}
//
//	//select_geodesic(out_target_curve, out_wrapper_curve);
//
//	return true;
//}
//
//bool ClusteredGeodesicsController::select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve)
//{
//	int geodesic_index = random_geodesics[data_model.patches.size()];
//	Eigen::MatrixXd original_geodesic = data_model.geodesics_candidates.paths[geodesic_index];
//	
//	//resample geodesic on target to match wrapper vertices
//	double resampled_target_length;
//	curve::resample_uniformly(original_geodesic, out_target_curve, resampled_target_length);
//
//	double target_off = abs(data_model.geodesics_candidates.lengths(geodesic_index) - resampled_target_length);
//
//	if (target_off > 1e-2)
//		write_log(2) << endl << endl << "WARN: resampled *target* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(target_off) << endl << endl;
//
//
//	double resampled_wrapper_length;
//	curve::resample_matching_curve(out_target_curve, 0, sample_direction, out_wrapper_curve, resampled_wrapper_length);
//
//	double wrapper_off = abs(resampled_target_length - resampled_wrapper_length);
//	if(wrapper_off > 1e-2)
//		write_log(2) << endl << endl << "WARN: resampled *wrapper* geodesic has not the same length! Likely a bug, go fix! Length difference is " << std::to_string(wrapper_off) << endl << endl;
//
//	return true;
//}
//
