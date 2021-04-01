#pragma once

#include <Eigen/Dense>
#include <vector>

//#include "DataModel.h"
class DataModel;

bool try_get_geodesic_between(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int& source_index, const int& target_index, Eigen::MatrixXd& out_path, double& out_length);

int get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point);
void get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point, int& out_index);
void get_closest_vertex(const Eigen::MatrixXd& vertices, const Eigen::RowVectorXd& point, int& out_index, Eigen::RowVectorXd& out_vertex);
int get_central_vertex(const std::vector<int>& indices, const Eigen::MatrixXd& lookup);

void get_face_normal(const Eigen::RowVectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::RowVectorXd& out_normal);
double normals_deviation(const int element_index, const int element_index_compare, const Eigen::MatrixXd& N);
Eigen::VectorXd get_surface_normal(const Eigen::VectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
Eigen::VectorXd get_surface_normal(const Eigen::VectorXd& point, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& N, const std::vector<std::vector<int>>& adjacency_VF);

Eigen::MatrixXd get_surface_curve_frame_at(const int& curve_index, const Eigen::MatrixXd& wrapper_curve, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

double angle(const Eigen::VectorXd& vector1, const Eigen::VectorXd& vector2);
Eigen::RowVector3d get_centroid(Eigen::MatrixXd& vertices);
std::vector<int> get_n_ring_neighborhood(const int start_vertex, const int number_rings, const std::vector<std::vector<int>>& adjacency_VV);

void append_vector(const Eigen::VectorXd& to_append, Eigen::VectorXd& out_target_vector);
void append_matrix(const Eigen::MatrixXd& to_append, Eigen::MatrixXd& out_longer_matrix);
void append_matrix(const Eigen::MatrixXi& to_append, Eigen::MatrixXi& out_longer_matrix);

std::vector<double> to_std_vector(const Eigen::VectorXd& eigen_vector);
Eigen::VectorXd to_eigen_vector(std::vector<double>& std_vector);

void index_to_value(const std::vector<int>& indices, const Eigen::MatrixXd& lookup, Eigen::MatrixXd& out_values);
void index_to_value(const Eigen::VectorXi& indices, const Eigen::MatrixXi& lookup, Eigen::MatrixXi& out_values);
void index_to_value(const std::vector<int>& indices, const Eigen::MatrixXi& lookup, Eigen::MatrixXi& out_values);
void index_to_value(const Eigen::VectorXi& indices, const Eigen::MatrixXd& lookup, Eigen::MatrixXd& out_values);
void value_to_index(const Eigen::MatrixXd& values, const Eigen::MatrixXd& lookup, std::vector<int>& out_indices);



bool is_contained(const int& element, const std::vector<int>& list);
int index_of(const int& element, const std::vector<int>& list);
int index_of(const double& element, const std::vector<double>& list);
std::vector<int> remove_duplicates(std::vector<int>& list);
std::vector<int> intersect(const std::vector<int>& list1, const std::vector<int>& list2);

void create_lookup(const Eigen::VectorXd& vector, Eigen::VectorXi& out_lookup);
std::vector<int> sort_indices_decending(const Eigen::VectorXd &v);

void convert_distribution_value_to_percent(const Eigen::VectorXd& distribution, const Eigen::VectorXi& distribution_lookup, const double& value, double& out_percentage);
void convert_distribution_percent_to_value(const Eigen::VectorXd& distribution, const Eigen::VectorXi& distribution_lookup, const double& percentage, double& out_value);


// TODO use igl::segment_segment_intersect
// TODO use igl::segments_intersect
bool line_line_intersection(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3, const Eigen::RowVector3d& p4, Eigen::RowVector3d& pa, Eigen::RowVector3d& pb);

void dfs(const std::vector<std::vector<int>>& search_list, const int& start_index, std::vector<int>& obstacles, std::vector<int>& out_component);
void dfs_within_list(std::vector<int>& search_list, const std::vector<std::vector<int>>& connectivity_list, const int& start_index, std::vector<int>& out_component);


//double acos(double cos_angle);
double clip(double n, double lower = -1, double upper = 1);
bool is_equal(double value1, double value2);

void compute_mean(const std::vector<double>& data, double& out_mean, double& out_stdev);
void compute_median(const std::vector<double>& data, double& out_median, double& out_q1, double& out_q3);