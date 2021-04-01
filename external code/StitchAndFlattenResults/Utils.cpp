#include "Utils.h"

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/barycentric_coordinates.h>
#include <igl/exact_geodesic.h>

#include <stack>

int get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point)
{
	Eigen::MatrixXd::Index index;
	(vertices.rowwise() - point).rowwise().squaredNorm().minCoeff(&index);

	return index;
}

void get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point, int& out_index)
{
	Eigen::RowVectorXd vertex;
	get_closest_vertex(vertices, point, out_index, vertex);
}

void get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point, int& out_index, Eigen::RowVectorXd& out_vertex)
{
	Eigen::MatrixXd::Index index;
	(vertices.rowwise() - point).rowwise().squaredNorm().minCoeff(&index);
	
	out_index = index;
	out_vertex = vertices.row(out_index);
}

int get_central_vertex(const std::vector<int>& indices, const Eigen::MatrixXd& lookup)	
{
	Eigen::MatrixXd values;
	index_to_value(indices, lookup, values);

	auto centroid = get_centroid(values);
	int center_index = get_closest_vertex(lookup, centroid);

	return center_index;
}

double angle(const Eigen::VectorXd& vector1, const Eigen::VectorXd& vector2)
{
	double cos_angle = clip(vector1.normalized().dot(vector2.normalized()), -1, 1);
	double angle = acos(cos_angle);
	return angle;
}

Eigen::VectorXd get_surface_normal(const Eigen::VectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	std::vector<std::vector<int>> VFi_unused;
	std::vector<std::vector<int>> adjacency_VF;
	igl::vertex_triangle_adjacency(V, F, adjacency_VF, VFi_unused);
	
	return get_surface_normal(point, V, F, N, adjacency_VF);
}

Eigen::VectorXd get_surface_normal(const Eigen::VectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& N, const std::vector<std::vector<int>>& adjacency_VF)
{
	int index = get_closest_vertex(V, point);
	//write_log(5) << "get_surface_normal:: point = " << point.transpose() << " which is close to v" << index << " (" << target.V.row(index) << ")" << endl;

	Eigen::RowVector3d barycentric_weights;
	Eigen::RowVector3i face;

	for (int fi : adjacency_VF[index])
	{
		face = F.row(fi);
		igl::barycentric_coordinates(point.transpose(), V.row(face(0)), V.row(face(1)), V.row(face(2)), barycentric_weights);
		//write_log(5) << "barycentric_weights: " << barycentric_weights << " at face: " << face << endl;

		if (barycentric_weights(0) >= 0.0 && barycentric_weights(0) <= 1.0 &&
			barycentric_weights(1) >= 0.0 && barycentric_weights(1) <= 1.0 &&
			barycentric_weights(2) >= 0.0 && barycentric_weights(2) <= 1.0)
		{
			break;
		}
	}


	Eigen::Vector3d normal(0, 0, 0);

	for (int bary_i = 0; bary_i < 3; bary_i++)
	{
		double weight = barycentric_weights(bary_i);
		int vertex_index = face(bary_i);

		normal += weight * N.row(vertex_index);
	}

	return normal;
}

double normals_deviation(const int element_index, const int element_index_compare, const Eigen::MatrixXd& N)
{
	if (element_index < 0 || element_index_compare < 0)
		return 0.0;

	Eigen::Vector3d n = N.row(element_index);
	Eigen::Vector3d n_compare = N.row(element_index_compare);

	return angle(n, n_compare);
}

void get_face_normal(const Eigen::RowVectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::RowVectorXd& out_normal)
{
	Eigen::MatrixXd query_points(1, point.cols());
	query_points << point;

	Eigen::VectorXd sq_d;
	Eigen::VectorXi closest_faces;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(query_points, V, F, sq_d, closest_faces, closest_points);

	Eigen::MatrixXd N;
	igl::per_face_normals(V, F, N);

	out_normal = N.row(closest_faces(0)); // +point;
}

Eigen::MatrixXd get_surface_curve_frame_at(const int& curve_index, const Eigen::MatrixXd& curve, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::RowVectorXd normal = get_surface_normal(curve.row(curve_index), V, F);

	Eigen::RowVector3d T = (curve.row(curve_index) - curve.row(curve_index - 1)).normalized();
	Eigen::RowVector3d N = normal;
	Eigen::RowVector3d B = T.cross(N);

	Eigen::Matrix3d frame; 
	frame.col(0) = T; 
	frame.col(1) = N; 
	frame.col(2) = B;

	return frame;
}

Eigen::RowVector3d get_centroid(Eigen::MatrixXd& vertices)
{
	const auto min_point = vertices.colwise().minCoeff();
	const auto max_point = vertices.colwise().maxCoeff();
	auto centroid = (0.5*(min_point + max_point)).eval();

	return centroid;
}

std::vector<int> get_n_ring_neighborhood(const int start_vertex, const int number_rings, const std::vector<std::vector<int>>& adjacency_VV)
{
	const int loglevel = 6;

	std::vector<int> neighbors = adjacency_VV[start_vertex];
	//log_list(loglevel, neighbors, "neighbors: ", false);

	if (number_rings == 1)
		return neighbors;

	int ring = 1;
	int last_neighbor = 0;
	int next_neighbor = 0;

	while (ring < number_rings)
	{
		int number_neighbors = neighbors.size();
		for (int i = last_neighbor; i < number_neighbors; i++)
		{
			int neighbor = neighbors[i];
			neighbors.insert(neighbors.end(), adjacency_VV[neighbor].begin(), adjacency_VV[neighbor].end());

			//log_list(loglevel, neighbors, "neighbors [" + to_string(neighbor) + "]: ", false);
		}

		//write_log(loglevel) << endl;

		last_neighbor = number_neighbors;
		ring++;
	}

	std::sort(neighbors.begin(), neighbors.end());
	neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
	//log_list(loglevel, neighbors, "unique neighbors: ", false);

	//still contains start_vertex, check if it is a problem
	return neighbors;
}

void append_vector(const Eigen::VectorXd & to_append, Eigen::VectorXd & out_target_vector)
{
	if (!out_target_vector.rows())
	{
		out_target_vector = to_append;
		return;
	}

	const auto temp = out_target_vector;
	out_target_vector.resize(temp.rows() + to_append.rows());
	out_target_vector << temp, to_append;

}

void append_matrix(const Eigen::MatrixXd& to_append, Eigen::MatrixXd& out_longer_matrix)
{
	if (!out_longer_matrix.rows())
	{
		out_longer_matrix = to_append;
		return;
	}

	const auto temp = out_longer_matrix;
	out_longer_matrix.resize(temp.rows() + to_append.rows(), 3);
	out_longer_matrix << temp, to_append;
}

void append_matrix(const Eigen::MatrixXi& to_append, Eigen::MatrixXi& out_longer_matrix)
{
	if (!out_longer_matrix.rows())
	{
		out_longer_matrix = to_append;
		return;
	}

	const auto temp = out_longer_matrix;
	out_longer_matrix.resize(temp.rows() + to_append.rows(), 3);
	out_longer_matrix << temp, to_append;
}

std::vector<double> to_std_vector(const Eigen::VectorXd& eigen_vector)
{
	std::vector<double> std_vector(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
	return std_vector;
}

Eigen::VectorXd to_eigen_vector(std::vector<double>& std_vector)
{
	Eigen::VectorXd eigen_vector = Eigen::Map<Eigen::VectorXd>(std_vector.data(), std_vector.size());
	return eigen_vector;
}

bool is_in_range(double value, double bound1, double bound2)
{
	double upper = bound1 > bound2 ? bound1 : bound2;
	double lower = bound1 <= bound2 ? bound1 : bound2;

	if (value <= upper && value >= lower)
		return true;

	return false;
}

bool is_in_range(const Eigen::RowVector3d& value, const Eigen::RowVector3d& bound1, const Eigen::RowVector3d& bound2)
{
	bool is_x_in_range = is_in_range(value(0), bound1(0), bound2(0));
	bool is_y_in_range = is_in_range(value(1), bound1(1), bound2(1));
	bool is_z_in_range = is_in_range(value(2), bound1(2), bound2(2));

	return is_x_in_range && is_y_in_range && is_z_in_range;
}

/*
   Calculate the line segment PaPb that is the shortest route between
   two lines (not segments!) P1P2 and P3P4. 
	  Pa = P1 + mua (P2 - P1)
	  Pb = P3 + mub (P4 - P3)
   Return FALSE if no solution exists.
   from http://paulbourke.net/geometry/pointlineplane/
*/
bool line_line_intersection(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3, const Eigen::RowVector3d& p4, Eigen::RowVector3d& pa, Eigen::RowVector3d& pb)
{
	// Algorithm is ported from the C algorithm of 
	// Paul Bourke at http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline3d/

	const double epsilon = 1e-9;
	const int x = 0;
	const int y = 1;
	const int z = 2;

	Eigen::RowVector3d p43 = p4 - p3;
	if (p43.squaredNorm() < epsilon)
		return false;

	Eigen::RowVector3d p21 = p2 - p1;
	if (p21.squaredNorm() < epsilon)
		return false;

	Eigen::RowVector3d p13 = p1 - p3;

	double d1343 = p13(x) * p43(x) + p13(y) * p43(y) + p13(z) * p43(z);
	double d4321 = p43(x) * p21(x) + p43(y) * p21(y) + p43(z) * p21(z);
	double d1321 = p13(x) * p21(x) + p13(y) * p21(y) + p13(z) * p21(z);
	double d4343 = p43(x) * p43(x) + p43(y) * p43(y) + p43(z) * p43(z);
	double d2121 = p21(x) * p21(x) + p21(y) * p21(y) + p21(z) * p21(z);

	double denom = d2121 * d4343 - d4321 * d4321;
	if (abs(denom) < epsilon)
		return false;

	double numer = d1343 * d4321 - d1321 * d4343;

	double mua = numer / denom;
	double mub = (d1343 + d4321 * (mua)) / d4343;

	pa(x) = p1(x) + mua * p21(x);
	pa(y) = p1(y) + mua * p21(y);
	pa(z) = p1(z) + mua * p21(z);

	pb(x) = p3(x) + mub * p43(x);
	pb(y) = p3(y) + mub * p43(y);
	pb(z) = p3(z) + mub * p43(z);

	if(is_in_range(pa, p1, p2) && is_in_range(pb, p3, p4))
		return true;
	//if(is_in_range(pa, p1, p2))
	//	return true;
	//if(is_in_range(pb, p3, p4))
	//	return true;

	return false;
	//return true;
}

void dfs(const std::vector<std::vector<int>>& search_list, const int& start_index, std::vector<int>& obstacles, std::vector<int>& out_component)
{
	// Initially mark all verices as not visited 
	std::vector<bool> visited(search_list.size(), false);

	// Create a stack for DFS 
	std::stack<int> stack;

	// Push the current source node. 
	stack.push(start_index);
	
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
			out_component.push_back(index);
		}

		std::vector<int>::iterator iterator = std::find(obstacles.begin(), obstacles.end(), index);
		bool is_contained = iterator != obstacles.end();
		if (is_contained)
			continue;


		// Get all adjacent vertices of the popped vertex s 
		// If a adjacent has not been visited, then push it to the stack. 
		for (auto i = search_list[index].begin(); i != search_list[index].end(); ++i)
			if (!visited[*i])
				stack.push(*i);
	}
}

void dfs_within_list(std::vector<int>& search_list, const std::vector<std::vector<int>>& connectivity_list, const int& start_index, std::vector<int>& out_component)
{
	// Initially mark all verices as not visited 
	std::vector<bool> visited(connectivity_list.size(), false);

	// Create a stack for DFS 
	std::stack<int> stack;

	// Push the current source node. 
	stack.push(start_index);

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

			if (is_contained(index, search_list))
				out_component.push_back(index);
		}

		if (!is_contained(index, search_list))
			continue;

		// Get all adjacent vertices of the popped vertex s 
		// If a adjacent has not been visited, then push it to the stack. 
		for (auto i = connectivity_list[index].begin(); i != connectivity_list[index].end(); ++i)
			if (!visited[*i])
				stack.push(*i);
	}
}

//double acos(double cos_angle)
//{
//	return std::acos(clip(cos_angle, -1, 1));
//}

double clip(double n, double lower, double upper)
{
	return std::max(lower, std::min(n, upper));
}

bool is_equal(double value1, double value2)
{
	return abs(value1 - value2) <= 1e-6;
}

void index_to_value(const std::vector<int>& indices, const Eigen::MatrixXi& lookup, Eigen::MatrixXi& out_values)
{
	out_values.resize(indices.size(), lookup.cols());
	for (int i = 0; i < indices.size(); i++)
		out_values.row(i) = lookup.row(indices[i]);
}

void index_to_value(const Eigen::VectorXi& indices, const Eigen::MatrixXi& lookup, Eigen::MatrixXi& out_values)
{
	out_values.resize(indices.size(), lookup.cols());
	for (int i = 0; i < indices.size(); i++)
		out_values.row(i) = lookup.row(indices[i]);
}

void index_to_value(const std::vector<int>& indices, const Eigen::MatrixXd& lookup, Eigen::MatrixXd& out_values)
{
	out_values.resize(indices.size(), lookup.cols());
	for (int i = 0; i < indices.size(); i++)
		out_values.row(i) = lookup.row(indices[i]);
}

void index_to_value(const Eigen::VectorXi& indices, const Eigen::MatrixXd& lookup, Eigen::MatrixXd& out_values)
{
	out_values.resize(indices.rows(), lookup.cols());
	for (int i = 0; i < indices.rows(); i++)
		out_values.row(i) = lookup.row(indices(i));
}

void value_to_index(const Eigen::MatrixXd& values, const Eigen::MatrixXd& lookup, std::vector<int>& out_indices)
{
	if (values.cols() != lookup.cols())
		return;

	out_indices.clear();
	for (int i = 0; i < values.rows(); i++)
	{
		int index;
		get_closest_vertex(lookup, values.row(i), index);
		out_indices.push_back(index);
	}
}

bool is_contained(const int& element, const std::vector<int>& list)
{
	auto iterator = std::find(list.begin(), list.end(), element);
	bool is_in_list = iterator != list.end();
	
	return is_in_list;
}

int index_of(const int& element, const std::vector<int>& list)
{
	auto iterator = std::find(list.begin(), list.end(), element);
	if (iterator == list.end()) //not found
		return -1;

	int index = std::distance(list.begin(), iterator);
	return index;
}

int index_of(const double& element, const std::vector<double>& list)
{
	auto iterator = std::find(list.begin(), list.end(), element);
	if (iterator == list.end()) //not found
		return -1;

	int index = std::distance(list.begin(), iterator);
	return index;
}

std::vector<int> remove_duplicates(std::vector<int>& list)
{
	std::sort(list.begin(), list.end());
	list.erase(std::unique(list.begin(), list.end()), list.end());
	return list;
}

std::vector<int> intersect(const std::vector<int>& list1, const std::vector<int>& list2)
{
	std::vector<int> intersection_list(list1.size() + list2.size());
	auto iterator = std::set_intersection(list1.begin(), list1.end(), list2.begin(), list2.end(), intersection_list.begin());
	intersection_list.resize(iterator - intersection_list.begin());

	return intersection_list;
}

void create_lookup(const Eigen::VectorXd& vector, Eigen::VectorXi& out_lookup)
{
	std::vector<int> indices_sorted = sort_indices_decending(vector);
	out_lookup = Eigen::Map<Eigen::VectorXi>(indices_sorted.data(), indices_sorted.size());
}

std::vector<int> sort_indices_decending(const Eigen::VectorXd &v)
{
	// initialize original index locations
	std::vector<int> indices(v.rows());
	std::iota(indices.begin(), indices.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(indices.begin(), indices.end(),
		[&v](int i1, int i2)
	{
		return v(i1) > v(i2);
	});

	return indices;
}

void convert_distribution_value_to_percent(const Eigen::VectorXd& distribution, const Eigen::VectorXi& distribution_lookup, const double& value, double& out_percentage) 
{
	const double max_value = distribution(distribution_lookup(0));

	if (value > max_value)
		return;

	const int number_values = distribution.rows();
	const Eigen::VectorXd ones = Eigen::VectorXd::Ones(number_values);

	int index;
	(distribution - ones * value).cwiseAbs().minCoeff(&index);

	int lookup_index;
	(distribution_lookup - ones.cast<int>()*index).cwiseAbs().minCoeff(&lookup_index);

	out_percentage = lookup_index / (double)number_values;
}

void convert_distribution_percent_to_value(const Eigen::VectorXd& distribution, const Eigen::VectorXi& distribution_lookup, const double& percentage, double& out_value)
{
	const double max_value = distribution(distribution_lookup(0));

	if (percentage < 1e-4)
	{
		out_value = max_value;
	}
	else
	{
		int lookup_index = percentage * (distribution.rows() - 1);
		out_value = distribution(distribution_lookup(lookup_index));
	}
}

bool try_get_geodesic_between(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int& source_index, const int& target_index, Eigen::MatrixXd& out_path, double& out_length)
{
	Eigen::VectorXi source_faces, target_faces;
	Eigen::VectorXi source_indices(1);
	Eigen::VectorXi target_indices(1);
	source_indices << source_index;
	target_indices << target_index;

	std::vector<Eigen::MatrixXd> paths;
	Eigen::VectorXd lengths;
	igl::exact_geodesic(V, F, source_indices, source_faces, target_indices, target_faces, lengths, paths);

	if (paths.size() < 1)
		return false;

	out_path = paths[0];
	out_length = lengths(0);

	return true;
}

void compute_mean(const std::vector<double>& data, double& out_mean, double& out_stdev)
{
	double sum = 0.0;
	double sum_squared = 0.0;

	for (double value : data)
	{
		sum += value;
		sum_squared += value * value;
	}

	out_mean = sum / data.size();
	double variance = (sum_squared / data.size()) - (out_mean * out_mean);
	out_stdev = std::sqrt(variance);
}

void compute_median(const std::vector<double>& data, double& out_median, double& out_q1, double& out_q3)
{
	std::vector<double> sorted(data);
	std::sort(sorted.begin(), sorted.end());

	int num = sorted.size();
	int n = num % 2 == 0 ? num / 2 : (num - 1) / 2;
	out_median = num % 2 == 0 ? (sorted[n - 1] + sorted[n]) / 2 : sorted[n];

	int nq = n % 2 == 0 ? n / 2 : (n - 1) / 2;
	out_q1 = n % 2 == 0 ? (sorted[nq - 1] + sorted[nq]) / 2 : sorted[nq];
	out_q3 = n % 2 == 0 ? (sorted[n + nq - 1] + sorted[n + nq]) / 2 : sorted[num - 1 - nq];
}


/* TODO write to file */
//ofstream file("RESULT_input_output.txt");
//file << "midpoints=" << to_mathematica(F_midpoints) << ";" << endl;
//file << "areas=" << to_mathematica(F_areas) << ";" << endl;

//file << "distances={";
//for (int i = 0; i < data_model.patches.size(); i++)
//	file << "{" << to_mathematica(patch_face_distances[i]) << "}," << endl;
//file << "};" << endl;

//file << "resultDistances=" << to_mathematica(resulting_distances) << ";" << endl;
//file << "resultLabels=" << to_mathematica(data_model.vertex_patch_assignment) << ";" << endl;
//file.close();
