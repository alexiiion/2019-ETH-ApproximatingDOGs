#pragma once

#include <string>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

class DataModel;
struct ViewModel;
class Patch;
class RuledDevelopableSurface;

class Mesh;
class WrapperMesh;

namespace meshhelper 
{
	void add_new_target(const std::string& mesh_path, DataModel& data_model, int& out_view_index_target, int& out_view_index_result);

	void add_target(const std::string& mesh_path, DataModel& data_model, int& out_view_index_target, int& out_view_index_result);
	void update_target_scale(DataModel& data_model, float& scale_to_dimension, bool do_center = true);

	void add_wrapper(const int width, const int height, DataModel& data_model, WrapperMesh& out_wrapper);
	void create_wrapper(const int width, const int height, DataModel& data_model, WrapperMesh& out_wrapper);
	void add_wrapper(DataModel& data_model, WrapperMesh& wrapper);

	void add_mesh(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool do_center_mesh = false, bool align_camera = true);
	void add_gauss_map(DataModel& data_model, ViewModel& view_model);
	
	void implicit_smooth_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& out_smooth_V, Eigen::MatrixXi& out_smooth_F, Eigen::MatrixXd& out_smooth_N);
	void compute_mesh_properties(Mesh& mesh, bool do_upsample = true);
	void calculate_dimensions(const Eigen::MatrixXd& vertices, Eigen::RowVectorXd& out_dimensions);
	std::vector<int> boundary_loop_flat(const Eigen::MatrixXi& F);

	Eigen::MatrixXi face_edges(const Eigen::MatrixXi& F);
	std::vector<int> get_faces_from_vertices(const std::vector<int>& component_vertices, const std::vector<std::vector<int>>& adjacency_VF);
	Eigen::MatrixXd compute_face_midpoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	Eigen::VectorXd compute_face_areas(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	double compute_face_area(const std::vector<int>& faces, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	double compute_face_area(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);


	std::vector<Mesh> cut_mesh(const std::vector<std::vector<int>>& vertex_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original, const std::vector<std::vector<int>>& adjacency_VF);
	Mesh sub_mesh(std::vector<int>& vertex_indices, std::vector<int>& face_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original);
	Mesh sub_mesh(std::vector<int>& vertex_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original, const std::vector<std::vector<int>>& adjacency_VF);


	std::vector<Mesh> get_wrapper_meshes(const std::vector<Patch*>& wrappers);
	std::vector<Mesh> get_ruled_developable_meshes(const std::vector<RuledDevelopableSurface*>& ruled_developables);
	std::vector<Mesh> get_selected_ruled_developable_meshes(const std::vector<int>& selected_ruled_vertex_indices, const std::vector<int>& ruled_vertex_indices, const std::vector<RuledDevelopableSurface*> &ruled_developables);

	void concatenate_meshes(const std::vector<Mesh>& meshes, Mesh& out_concatenated, std::vector<int>& out_cumulative_offset_vertices, std::vector<int>& out_cumulative_offset_faces);
	Mesh concatenate_meshes(const std::vector<Mesh>& meshes);
	Mesh concatenate_meshes(const std::vector<Eigen::MatrixXd>& list_V, const std::vector<Eigen::MatrixXi>& list_F);
	Mesh concatenate_ruled_developables(const std::vector<RuledDevelopableSurface*>& ruled_developables);
	Mesh concatenate_wrappers(const std::vector<Patch*>& wrappers);

	void concatenate_wrappers(const DataModel& data_model, Eigen::MatrixXd& out_concatenated_wrappers_V, Eigen::MatrixXi& out_concatenated_wrappers_F);
	
	
	void project_to_meshes(const Eigen::MatrixXd& source_points, const std::vector<Mesh>& meshes, std::vector<Eigen::VectorXd>& out_patch_point_distances, std::vector<Eigen::MatrixXd>& out_patch_closest_points);
	void project_to_meshes(const Eigen::MatrixXd& source_points, const Eigen::VectorXd& weights, const std::vector<Mesh>& meshes, std::vector<Eigen::VectorXd>& out_patch_point_distances, std::vector<Eigen::MatrixXd>& out_patch_closest_points);

	void project_to_mesh(const Eigen::MatrixXd& source_points, const Mesh& mesh, Eigen::VectorXd& out_point_distances, Eigen::MatrixXd& out_closest_points);
	void project_to_mesh(Eigen::MatrixXd& points, const Mesh& projection_target);
	
	void simple_project_to_DOGs(const std::vector<Patch*>& surfaces, Mesh& to_deform);
	void simple_project_to_ruled_surfaces(const std::vector<RuledDevelopableSurface*>& surfaces, Mesh& to_deform);

	////project after graph cut label selection
	//void update_target_developable_optimized(DataModel& data_model);
	//void update_target_developable_face_optimized(DataModel& data_model);
	//void update_target_developable_vertex_optimized(DataModel& data_model);

	//void update_covering_wrappers(DataModel& data_model, Eigen::MatrixXi& out_concatenated_covering_F);
	//void update_target_coverage(DataModel& data_model, const Eigen::MatrixXi& concatenated_covering_F);
	//void update_coverage_per_vertices(DataModel& data_model); //if above works, remove this!
	//void get_target_coverage(const DataModel& data_model, const Eigen::MatrixXd& source_points, Eigen::VectorXi& out_covered, Eigen::VectorXi& out_uncovered, Eigen::VectorXd& out_distances_squared, Eigen::MatrixXd& out_closest_points, Eigen::VectorXd& out_covering_wrapper_faces_distances);
}

/*
namespace surface
{
	void compute_weighted_features(Mesh& target, double weight_K, double weight_crease, double weight_H);
	void compute_features(Mesh& target);

	void compute_principal_k(Mesh& target);
	void compute_K(Mesh& target);
	void compute_H(Mesh& target);
	void compute_creases(Mesh& target);

	void get_gaussian_curvature(Mesh& target, Eigen::VectorXd& out_K);
	void get_mean_curvature(Mesh& target, Eigen::VectorXd& out_H);
	void normalize_vector_entries(const Eigen::VectorXd& curvature, Eigen::VectorXd& out_normalized_curvature);
}
*/