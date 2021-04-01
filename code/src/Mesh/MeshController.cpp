#include "MeshController.h"

#include <igl/avg_edge_length.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/boundary_loop.h>
#include <igl/hsv_to_rgb.h>

#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

#include <igl/point_mesh_squared_distance.h>
#include <igl/upsample.h>
#include <igl/loop.h>

#include "DataModel.h"
#include "ViewModel.h"
#include "MeshModel.h"
#include "Coverage.h"
#include "DevelopableResultController.h"

#include "SurfaceFeatures.h"
#include "CreaseFinder.h"
#include "GraphCut.h"
#include "PatchModel.h"
#include "Quad.h"
#include "Logger.h"
#include "Utils.h"
#include "ViewUtils.h"



using namespace std;

void center_mesh(Eigen::MatrixXd& V);
//void find_result_creases(DataModel& data_model);
//void find_labeling_seams(DataModel& data_model, const std::vector<int>& face_patch_assignment, Eigen::MatrixXi& face_edges, const Eigen::MatrixXi& adjacency_FF);
//void find_labeling_seams(DataModel& data_model, const std::vector<int>& face_patch_assignment);

double upsample_to_resolution(const double max_edge_length, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	double current_edge_length = igl::avg_edge_length(V, F);

	write_log(0) << "before upsampling:: current_edge_length = " << current_edge_length << ", # vertices = " << V.rows() << endl;
	while (current_edge_length > max_edge_length)
	{
		igl::loop(V, F, V, F, 1);
		current_edge_length = igl::avg_edge_length(V, F);
		write_log(0) << "upsampled:: current_edge_length = " << current_edge_length << ", # vertices = " << V.rows() << endl;
	}
	return current_edge_length;
}

void meshhelper::implicit_smooth_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& out_smooth_V, Eigen::MatrixXi& out_smooth_F, Eigen::MatrixXd& out_smooth_N)
{
	const double vertexSmooth = 0.001; //0.0001;
	const double normalSmooth = 0.01; //0.001;
	//const double max_edge_length = 1.20; //0.75;

	out_smooth_V = V;
	out_smooth_F = F;

	//upsample_to_resolution(max_edge_length, out_smooth_V, out_smooth_F);
	igl::per_vertex_normals(out_smooth_V, out_smooth_F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, out_smooth_N);

	Eigen::SparseMatrix<double> L, M;
	igl::cotmatrix(out_smooth_V, out_smooth_F, L);
	igl::massmatrix(out_smooth_V, out_smooth_F, igl::MASSMATRIX_TYPE_VORONOI, M);

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholv(M - vertexSmooth * L);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> choln(M - normalSmooth * L);

	out_smooth_V = cholv.solve(M * out_smooth_V);
	out_smooth_N = choln.solve(M * out_smooth_N);
}

void meshhelper::compute_mesh_properties(Mesh& mesh, bool do_upsample)
{
	if (do_upsample)
	{
		double edge_length = upsample_to_resolution(GlobalSettings::max_average_edge, mesh.V, mesh.F);
		mesh.average_edge_length = edge_length;

		write_log(0) << "after upsampling:: average_edge_length = " << mesh.average_edge_length << ", # vertices = " << mesh.V.rows() << endl;
	}
	else
	{
		mesh.average_edge_length = igl::avg_edge_length(mesh.V, mesh.F);
	}
	//mesh.average_edge_length = igl::avg_edge_length(mesh.V, mesh.F);

	igl::per_face_normals(mesh.V, mesh.F, mesh.normals_faces);

	//store topology, normals, curvatures
	std::vector<std::vector<int>> VFi_unused;
	igl::vertex_triangle_adjacency(mesh.V, mesh.F, mesh.adjacency_VF, VFi_unused);
	igl::triangle_triangle_adjacency(mesh.F, mesh.adjacency_FF);
	igl::adjacency_list(mesh.F, mesh.adjacency_VV);

	mesh.surface_features();
	mesh.surface_area = meshhelper::compute_face_area(mesh.V, mesh.F);
}

void meshhelper::add_gauss_map(DataModel& data_model, ViewModel& view_model)
{
	//get config file path
	string file = __FILE__;
	auto index = file.find_last_of("/\\");
	auto path = file.substr(0, index + 1);
	auto sphere_filepath = path + "/../../data/sphere.obj";

	igl::readOBJ(sphere_filepath, view_model.gauss_map_sphere.V, view_model.gauss_map_sphere.F);
	add_mesh(data_model.viewer, view_model.gauss_map_sphere.V, view_model.gauss_map_sphere.F, false, false);

	data_model.viewer.data().clear();
	//data_model.viewer.data().set_mesh(view_model.gauss_map_sphere.V, view_model.gauss_map_sphere.F);
	view_model.gauss_map_view_index = data_model.viewer.selected_data_index;

	Eigen::Vector3d diffuse; diffuse << 0.98, 0.98, 0.98;
	Eigen::Vector3d ambient; ambient << 0, 0, 0;//0.05*diffuse;
	Eigen::Vector3d specular; specular << 0, 0, 0;
	data_model.viewer.data().uniform_colors(ambient, diffuse, specular);
	data_model.viewer.data().show_lines = false;
	//viewer.core().shininess = 0;

	data_model.viewer.data().point_size = 4.0;
	data_model.viewer.data().line_width = 1.0;

	Eigen::RowVectorXd dimensions;
	calculate_dimensions(view_model.gauss_map_sphere.V, dimensions);
	double current_radius = dimensions(0) * 0.5;

	double scale = view_model.gauss_map_sphere_radius / current_radius;
	view_model.gauss_map_sphere.V *= scale;

	//if (switched_mode) viewer.core().align_camera_center(dog->getVrendering(), dog->getFrendering());
	////viewer.core().align_camera_center(sphereV, sphereF);
	////viewer.core().show_lines = false;

	////Michael's code for scaling the sphere, not sure if neccessary
	//Eigen::RowVectorXd colmean = GV.colwise().mean();
	//for (int i = 0; i < GV.rows(); i++) {
	//	GV.row(i) = GV.row(i) - colmean; // TODO can't this be done in 1 line?
	//}
	//Eigen::VectorXd area_v;
	//igl::doublearea(GV, GF, area_v);
	//double area = area_v.sum() / 2.;
	//double eps = 2e-1;
	//double scale = sqrt((4 - eps)*M_PI / area); // make it a little bit smaller so we could see the lines
	//GV = GV * scale;
}

void meshhelper::add_new_target(const std::string& mesh_path, DataModel& data_model, int& out_view_index_target, int& out_view_index_result)
{
	data_model.target_filename = mesh_path;

	//Eigen::MatrixXd corner_normals; //??
	Eigen::MatrixXi fNormIndices;
	Eigen::MatrixXd UV_V;
	Eigen::MatrixXi UV_F;


	NewMesh test;
	bool success = igl::readOBJ(data_model.models_folder + mesh_path, test.V(), test.F());


	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd NV;

	bool success2 = igl::readOBJ(data_model.models_folder + mesh_path, V, UV_V, NV, F, UV_F, fNormIndices);
	NewMesh test2(V, F, NV);


	write_log(0) << "test direct method passing:" << linebreak << "V: " << test.V() << linebreak << "F: " << test.F() << endl;
	write_log(0) << "test2:" << linebreak << "V: " << test2.V() << linebreak << "F: " << test2.F() << endl;

	Eigen::MatrixXd NV2 = test.NV();
	write_log(0) << "NV2:" << linebreak << NV2 << endl;

		//<<
}

void meshhelper::add_target(const std::string& mesh_path, DataModel& data_model, int& out_view_index_target, int& out_view_index_result)
{
	data_model.target_filename = mesh_path;

	//Eigen::MatrixXd corner_normals; //??
	Eigen::MatrixXi fNormIndices;
	Eigen::MatrixXd UV_V;
	Eigen::MatrixXi UV_F;

	//igl::readOBJ(data_model.models_folder + mesh_path, data_model.target.V, data_model.target.F);
	bool success = igl::readOBJ(data_model.models_folder + mesh_path, data_model.target_original.V, UV_V, data_model.target_original.normals_vertices, data_model.target_original.F, UV_F, fNormIndices);
	//if (!success)
	//	return;

	//read or create normals
	//if(data_model.target_original.normals_vertices.rows() < 1)
		igl::per_vertex_normals(data_model.target_original.V, data_model.target_original.F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, data_model.target_original.normals_vertices);
	//igl::per_face_normals(data_model.target_original.V, data_model.target_original.F, data_model.target_original.normals_faces);


	/*
	//translate (center) and scale mesh
	if (abs(data_model.target_max_dimension) < 1e-6)
		data_model.target_max_dimension = 20;

	Eigen::RowVectorXd dimensions;
	calculate_dimensions(data_model.target_original.V, dimensions);
	write_log(4) << "dimensions: " << dimensions(0) << " x " << dimensions(1) << " x " << dimensions(2) << " // diagonal: " << dimensions.norm() << endl;

	center_mesh(data_model.target_original.V);

	double scale = data_model.target_max_dimension / dimensions.maxCoeff();
	//double scale = 30.0 / dimensions(0);
	data_model.target_original.V *= scale;
	write_log(4) << "scale = " << scale << ";" << endl;
	write_log(4) << "surfaceWidth = " << data_model.ruled_width << ";" << endl;

	calculate_dimensions(data_model.target_original.V, dimensions);
	write_log(4) << "scaled diagonal: " << dimensions.norm() << endl;



	compute_mesh_properties(data_model.target_original);

	////add to view
	//add_mesh(data_model.viewer, data_model.target_original.V, data_model.target_original.F, false, true);


	data_model.target = data_model.target_original;

	upsample_to_resolution(GlobalSettings::max_average_edge, data_model.target_original.V, data_model.target_original.F);
	//implicit_smooth_mesh(data_model.target_original.V, data_model.target_original.F, data_model.target.V, data_model.target.F, data_model.target.normals_vertices);
	compute_mesh_properties(data_model.target);
	*/
	
	meshhelper::update_target_scale(data_model, data_model.target_max_dimension, true);

	//add to view
	add_mesh(data_model.viewer, data_model.target_original.V, data_model.target_original.F, false, true);



	//update_coverage_per_vertices(data_model);
	data_model.initialization_coverage = new Coverage(data_model.target, &(data_model.optimization_settings.coverage_threshold));
	data_model.patch_coverage = new Coverage(data_model.target, &(data_model.optimization_settings.coverage_threshold));
	data_model.result_coverage = new Coverage(data_model.target, &(data_model.optimization_settings.coverage_threshold));

	//data_model.optimization_settings.outlier_threshold = 0.01 * data_model.target_max_dimension;
	//data_model.optimization_settings.coverage_threshold = 1.2 * data_model.target.average_edge_length;
	data_model.optimization_settings.coverage_threshold = 0.8;
	data_model.optimization_settings.outlier_threshold = 0.33 * data_model.target.average_edge_length;
	write_log(4) << "edge length: " << data_model.target.average_edge_length << ", outliers: " << data_model.optimization_settings.outlier_threshold << ", coverage: " << data_model.optimization_settings.coverage_threshold << endl;



	//DEBUG this currently overrides the original mesh in the view
	//TODO for the final version, show orignial and only compute with this target
	data_model.viewer.data().clear();
	data_model.viewer.data().set_mesh(data_model.target.V, data_model.target.F);
	data_model.viewer.data().set_normals(data_model.target.normals_vertices);
	out_view_index_target = data_model.viewer.selected_data_index;
	
	data_model.viewer.data().point_size = 5.0;
	data_model.viewer.data().line_width = 0.5;

	Eigen::Vector4f line_color(Colors::GRAY_MID(0), Colors::GRAY_MID(1), Colors::GRAY_MID(2),  1.0f);
	//data_model.viewer.data().line_width = 1.5;
	//Eigen::Vector4f line_color(Colors::GRAY_DARK(0), Colors::GRAY_DARK(1), Colors::GRAY_DARK(2),  1.0f);
	data_model.viewer.data().line_color = line_color;




	//copy target mesh to target_developable, which will be deformed
	//data_model.target_developable.V = data_model.target.V;
	//data_model.target_developable.F = data_model.target.F;
	data_model.developable_model.non_developable = data_model.target;

	//add_mesh(data_model.viewer, data_model.target_developable.V, data_model.target_developable.F, do_center_mesh, true);
	add_mesh(data_model.viewer, data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F, false, true);

	//REFACTOR view.update_mesh(index, V, F);		
	data_model.viewer.data().clear();
	//data_model.viewer.data().set_mesh(data_model.target_developable.V, data_model.target_developable.F);
	data_model.viewer.data().set_mesh(data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F);
	out_view_index_result = data_model.viewer.selected_data_index;

	data_model.viewer.data().point_size = 5.0;
	data_model.viewer.data().line_width = 0.5;
	data_model.viewer.data().line_color = line_color;

	data_model.vertex_labels.resize(data_model.target.V.rows());
	data_model.vertex_labels.setZero();
	data_model.face_labels.resize(data_model.target.F.rows());
	data_model.face_labels.setZero();
}

void meshhelper::update_target_scale(DataModel& data_model, float& scale_to_dimension, bool do_center)
{
	//translate (center) and scale mesh

	//if (abs(data_model.target_max_dimension) < 1e-6)
	//	data_model.target_max_dimension = 20;

	if(do_center)
		center_mesh(data_model.target_original.V);

	Eigen::RowVectorXd dimensions;
	calculate_dimensions(data_model.target_original.V, dimensions);
	write_log(4) << "(original) dimensions: " << dimensions(0) << " x " << dimensions(1) << " x " << dimensions(2) << " // diagonal: " << dimensions.norm() << endl;

	if (abs(data_model.target_max_dimension) >= 1e-6)
	{
		double scale = data_model.target_max_dimension / dimensions.maxCoeff();
		data_model.target_original.V *= scale;
		write_log(4) << "  --> scale = " << scale << ";" << endl;
		//write_log(4) << "  surfaceWidth = " << data_model.ruled_width << ";" << endl;
	}

	calculate_dimensions(data_model.target_original.V, dimensions);
	write_log(4) << "(scaled) dimensions: " << dimensions(0) << " x " << dimensions(1) << " x " << dimensions(2) << " // diagonal: " << dimensions.norm() << linebreak << endl;

	data_model.target_dimensions = dimensions;
	data_model.target_diagonal = dimensions.norm();
	data_model.target_max_dimension = dimensions.maxCoeff();

	//compute all mesh properties
	compute_mesh_properties(data_model.target_original, false);


	data_model.target = data_model.target_original;

	//upsample_to_resolution(GlobalSettings::max_average_edge, data_model.target.V, data_model.target.F);
	//implicit_smooth_mesh(data_model.target_original.V, data_model.target_original.F, data_model.target.V, data_model.target.F, data_model.target.normals_vertices);
	compute_mesh_properties(data_model.target, true);
}

void meshhelper::create_wrapper(const int width, const int height, DataModel& data_model, WrapperMesh& out_wrapper)
{
	out_wrapper.quad_width = width;
	out_wrapper.quad_height = height;

	//load_wrapper_mesh(data_model.viewer, out_wrapper.quad_width, out_wrapper.quad_height, do_center_mesh);
	get_planar_square_mesh(out_wrapper.V, out_wrapper.F, height, width);
	out_wrapper.F_quad = F_to_Fsqr(out_wrapper.F);

	if (!out_wrapper.F.rows() || !out_wrapper.F_quad.rows())
		write_log(1) << "error at creating wrapper mesh!" << endl;

	// get mesh connectivity data structure (edges, stars, etc..)
	quad_topology(out_wrapper.V, out_wrapper.F_quad, out_wrapper.quad_topology);
	//out_wrapper.render_F = out_wrapper.F;
	//out_wrapper.render_F_quad = out_wrapper.F_quad;

	// Globally scale and normalize mesh size (important for optimization) such that a (random) edge on the boundary will be of length 1
	const double edge_l = (out_wrapper.V.row(out_wrapper.quad_topology.bnd_loop[1]) - out_wrapper.V.row(out_wrapper.quad_topology.bnd_loop[0])).norm();
	out_wrapper.V *= 1. / edge_l;
	write_log(5) << "wrapper edge length: " << edge_l << endl;

	mat2_to_vec(out_wrapper.V, out_wrapper.initial_x0);

	//igl::per_vertex_normals(out_wrapper.V, out_wrapper.F, out_wrapper.normals_vertices);
	//igl::per_face_normals(out_wrapper.V, out_wrapper.F, out_wrapper.normals_faces);
	//meshhelper::compute_mesh_properties(out_wrapper);
}

void meshhelper::add_wrapper(DataModel& data_model, WrapperMesh& wrapper)
{
	add_mesh(data_model.viewer, wrapper.V, wrapper.F, false, false);

	//TODO manage indices and colors
	data_model.viewer.data().set_colors(Eigen::RowVector3d(0.1, 0.9, 0.1));
	data_model.wrapper_view_indices.push_back(data_model.viewer.selected_data_index);

	data_model.viewer.data().point_size = 10.0;
	data_model.viewer.data().line_width = 0.5;
}

void meshhelper::add_wrapper(const int width, const int height, DataModel& data_model, WrapperMesh& out_wrapper)
{
	out_wrapper.quad_width = width;
	out_wrapper.quad_height = height;

	//load_wrapper_mesh(data_model.viewer, out_wrapper.quad_width, out_wrapper.quad_height, do_center_mesh);
	get_planar_square_mesh(out_wrapper.V, out_wrapper.F, height, width);
	out_wrapper.F_quad = F_to_Fsqr(out_wrapper.F);

	if (!out_wrapper.F.rows() || !out_wrapper.F_quad.rows())
		write_log(1) << "error at creating wrapper mesh!" << endl;

	// get mesh connectivity data structure (edges, stars, etc..)
	quad_topology(out_wrapper.V, out_wrapper.F_quad, out_wrapper.quad_topology);
	//out_wrapper.render_F = out_wrapper.F;
	//out_wrapper.render_F_quad = out_wrapper.F_quad;

	// Globally scale and normalize mesh size (important for optimization) such that a (random) edge on the boundary will be of length 1
	const double edge_l = (out_wrapper.V.row(out_wrapper.quad_topology.bnd_loop[1]) - out_wrapper.V.row(out_wrapper.quad_topology.bnd_loop[0])).norm();
	out_wrapper.V *= 1. / edge_l;
	write_log(5) << "wrapper edge length: " << edge_l << endl;

	mat2_to_vec(out_wrapper.V, out_wrapper.initial_x0);

	//vector<int> visible_quad_faces = vector<int>(out_wrapper.F_quad.rows());
	//std::iota(std::begin(visible_quad_faces), std::end(visible_quad_faces), 0);
	//out_wrapper.covering_quad_faces.push_back(visible_quad_faces);


	add_mesh(data_model.viewer, out_wrapper.V, out_wrapper.F, false, false);

	//TODO manage indices and colors
	data_model.viewer.data().set_colors(Eigen::RowVector3d(0.1, 0.9, 0.1));
	data_model.wrapper_view_indices.push_back(data_model.viewer.selected_data_index);

	data_model.viewer.data().point_size = 10.0;
	data_model.viewer.data().line_width = 0.5;

	//MeshViewSettings& view_settings = MeshViewSettings();
	//view_settings.view_index = data_model.viewer.selected_data_index;
	//view_model.patch_view_settings.push_back(view_settings);
}

Mesh meshhelper::sub_mesh(std::vector<int>& vertex_indices, std::vector<int>& face_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original)
{
	Eigen::MatrixXi mapped_faces(face_indices.size(), 3);
	int added_faces = 0;

	for (int i = 0; i < face_indices.size(); i++)
	{
		int fi = face_indices[i];
		Eigen::RowVectorXi face_vertices_global = F_original.row(fi);
		//write_log(4) << i << ": global face_vertices for fi=" << fi << ": " << face_vertices_global << endl;

		vector<int> face_vertices_local;
		for (int j = 0; j < face_vertices_global.size(); j++)
		{
			int vi = face_vertices_global(j);
			face_vertices_local.push_back(index_of(vi, vertex_indices));
		}

		if (is_contained(-1, face_vertices_local))
		{
			//write_log(4) << "   NOT including this face!" << endl;
			continue;
		}

		mapped_faces.row(added_faces) = Eigen::Map<Eigen::RowVectorXi>(face_vertices_local.data(), face_vertices_local.size());
		added_faces++;

		//write_log(4) << "local component_faces.row(" << i << "): " << component_faces.row(i) << endl;
	}

	mapped_faces.conservativeResize(added_faces, 3);
	//write_log(4) << "local component_faces: " << endl << component_faces << endl;

	Eigen::MatrixXd V;
	index_to_value(vertex_indices, V_original, V);

	Mesh mesh;
	mesh.V = V;
	mesh.F = mapped_faces;

	return mesh;
}

Mesh meshhelper::sub_mesh(std::vector<int>& vertex_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original, const std::vector<std::vector<int>>& adjacency_VF)
{
	vector<int> face_indices;
	for each (int vertex in vertex_indices)
	{
		vector<int> faces = adjacency_VF[vertex];
		face_indices.insert(face_indices.end(), faces.begin(), faces.end());
	}

	std::sort(face_indices.begin(), face_indices.end());
	face_indices.erase(std::unique(face_indices.begin(), face_indices.end()), face_indices.end());


	Mesh mesh = sub_mesh(vertex_indices, face_indices, V_original, F_original);
	return mesh;
}

std::vector<Mesh> meshhelper::cut_mesh(const std::vector<std::vector<int>>& vertex_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original, const std::vector<std::vector<int>>& adjacency_VF)
{
	//re-map faces to map to local uncovered area vertices
	std::vector<Eigen::MatrixXi> remapped_faces;

	for each (vector<int> component in vertex_indices)
	{
		vector<int> face_list;
		for each (int labeled_vertex in component)
		{
			vector<int> faces = adjacency_VF[labeled_vertex];
			face_list.insert(face_list.end(), faces.begin(), faces.end());
		}

		//log_list(4, face_list, "global face_list: ", false, false);

		std::sort(face_list.begin(), face_list.end());
		face_list.erase(std::unique(face_list.begin(), face_list.end()), face_list.end());

		//log_list(4, face_list, "global unique face_list: ", false, false);


		Eigen::MatrixXi component_faces(face_list.size(), 3);
		int added_faces = 0;

		for (int i = 0; i < face_list.size(); i++)
		{
			int fi = face_list[i];
			Eigen::RowVectorXi face_vertices_global = F_original.row(fi);
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
		remapped_faces.push_back(component_faces);
		//write_log(4) << "local component_faces: " << endl << component_faces << endl;
	}

	std::vector<Mesh> meshes;
	for (int i = 0; i < vertex_indices.size(); i++)
	{
		Eigen::MatrixXd V;
		index_to_value(vertex_indices[i], V_original, V);

		Mesh mesh;
		mesh.V = V;
		mesh.F = remapped_faces[i];

		meshes.push_back(mesh);
	}

	return meshes;
}

Eigen::MatrixXd meshhelper::compute_face_midpoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::MatrixXd F_midpoints(F.rows(), 3);

	for (int fi = 0; fi < F.rows(); fi++)
	{
		Eigen::RowVectorXi face = F.row(fi);
		Eigen::RowVectorXd face_sum = Eigen::RowVectorXd::Zero(3);

		for (int vi = 0; vi < 3; vi++)
		{
			if (face(vi) >= 0) //if it's -1 then there is no face because it's the boundary
				face_sum = face_sum + V.row(face(vi));
		}

		F_midpoints.row(fi) = (face_sum / 3).eval();
	}

	return F_midpoints;
}

Eigen::VectorXd meshhelper::compute_face_areas(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::VectorXd F_areas(F.rows());

	for (int fi = 0; fi < F.rows(); fi++)
	{
		Eigen::RowVectorXi face = F.row(fi);
		Eigen::RowVector3d a = V.row(face(1)) - V.row(face(0));
		Eigen::RowVector3d b = V.row(face(2)) - V.row(face(0));
		F_areas(fi) = a.cross(b).norm() / 2;
	}

	return F_areas;
}

double meshhelper::compute_face_area(const std::vector<int>& faces, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::VectorXd F_areas = meshhelper::compute_face_areas(V, F);

	double area = 0.0;
	for (int fi : faces)
		area += F_areas(fi);

	return area;
}

double meshhelper::compute_face_area(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	Eigen::VectorXd F_areas = meshhelper::compute_face_areas(V, F);

	double area = 0.0;
	for (int i  = 0; i < F.rows(); i++)
		area += F_areas(i);

	return area;
}

std::vector<int> meshhelper::get_faces_from_vertices(const std::vector<int>& component_vertices, const std::vector<std::vector<int>>& adjacency_VF)
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
void meshhelper::update_target_developable(DataModel& data_model)
{
	if (data_model.patches.size() < 1)
		return;
	
	//Eigen::MatrixXd upsampled_V;
	//Eigen::MatrixXi upsampled_F;
	//igl::upsample(data_model.target.V, data_model.target.F, data_model.target_developable.V, data_model.target_developable.F, 2);
	//igl::upsample(data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.V, data_model.target_developable.F);

	//update coverage
	concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);
	//update_coverage_per_vertices(data_model);

	Eigen::VectorXi face_indices_unused;
	Eigen::VectorXd distances_squared;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(data_model.target_developable.V, data_model.concatenated_patches_V, data_model.concatenated_patches_F, distances_squared, face_indices_unused, closest_points);

	//move vertices to wrappers
	data_model.target_developable.V = closest_points;
	

	////project vertices on each wrapper
	////assign wrapper(s) to each vertex (d < 3*coverage_threshold)
	//const double assignment_threshold = data_model.optimization_settings.coverage_threshold * 5;
	//vector<vector<int>> vertex_assigned_patches; // indices by V, then list of P
	//
	//for (int i = 0; i < data_model.target.V.rows(); i++)
	//	vertex_assigned_patches.push_back(vector<int>());

	//vector<Eigen::VectorXd> patches_vertex_distances; //indices by P then V
	//vector<Eigen::MatrixXd> patches_vertex_points;

	//for (int i = 0; i < data_model.patches.size(); i++)
	//{
	//	WrapperMesh& wrapper = data_model.patches[i]->wrapper;

	//	Eigen::VectorXi face_indices;
	//	Eigen::VectorXd distances_squared;
	//	Eigen::MatrixXd closest_points;
	//	igl::point_mesh_squared_distance(data_model.target.V, wrapper.V, wrapper.F, distances_squared, face_indices, closest_points);
	//	
	//	patches_vertex_distances.push_back(distances_squared);
	//	patches_vertex_points.push_back(closest_points);


	//	for (int j = 0; j < distances_squared.rows(); j++)
	//	{
	//		if (distances_squared(j) < assignment_threshold)
	//			vertex_assigned_patches[j].push_back(i);
	//	}
	//}


	////check vertices with more than one wrapper for discontinuities
	//for (int i = 0; i < vertex_assigned_patches.size(); i++)
	//{
	//	vector<int> patch_candidates = vertex_assigned_patches[i];

	//	if (patch_candidates.size() < 1)
	//	{
	//		write_log(0) << "error at vertex " << i << endl;
	//		continue;
	//	}

	//	if (patch_candidates.size() == 1)
	//	{
	//		data_model.target_developable.V = patches_vertex_points[patch_candidates[0]].row(i);
	//		continue;
	//	}

	//	vector<int> vertex_neighbors = data_model.target.adjacency_VV[i];
	//	for (int j = 0; j < vertex_neighbors.size(); j++)
	//	{
	//		//check neighbors for their patch assignment, somehow try to find best assignment for this vertex...
	//	}
	//}
}
*/

/*
void print_mathematica_vector(const Eigen::VectorXd vector)
{
	//cout << endl << endl;
	std::cout << "{";

	for (int i = 0; i < vector.rows(); i++)
	{
		std::cout << vector(i);
		if (i < vector.rows() - 1)
			std::cout << ", ";
	}

	std::cout << "}" << std::endl;
	//cout << endl << endl;
}

std::string to_mathematica(const Eigen::MatrixXd matrix)
{
	int rows = matrix.rows();
	int cols = matrix.cols();

	if (rows < 1 || cols < 1)
		return "";

	std::stringstream output;

	output << "{";
	for (int r = 0; r < rows; r++)
	{
		if (cols == 1)
		{
			output << std::fixed << matrix(r, 0);
		}
		else
		{
			output << "{";
			for (int c = 0; c < cols; c++)
			{
				output << std::fixed << matrix(r, c);
				if (c < cols - 1)
					output << ", ";
			}
			output << "}";
		}

		if (r < rows - 1)
			output << "," << std::endl;
	}
	output << "}";
	return output.str();
}

//std::string to_mathematica(const vector<double>& vector)
//{
//	std::stringstream output;
//
//	output << "{";
//
//	for (int i = 0; i < vector.size(); i++)
//	{
//		output << vector[i];
//		if (i < vector.size() - 1)
//			output << ", " << std::endl;
//	}
//
//	output << "};" << std::endl;
//	return output.str();
//}
*/


/*
void meshhelper::update_target_developable_optimized(DataModel& data_model)
{
	//DevelopableResultController result_controller(data_model);
	//result_controller.optimize_face_based();

	//meshhelper::update_target_developable_face_optimized(data_model);
}
*/

/*
void meshhelper::update_target_developable_face_optimized(DataModel& data_model)
{
	if (data_model.patches.size() < 1)
		return;

	//igl::upsample(data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.V, data_model.target_developable.F, 1);
	////igl::loop(data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.V, data_model.target_developable.F, 2);


	//do graph cut on faces instead of vertices
	//write_log(0) << "adjacency_FF.cols() = " << data_model.target.adjacency_FF.cols() << endl;
	//write_log(0) << data_model.target.adjacency_FF << endl;

	Eigen::MatrixXd F_midpoints(data_model.target_developable.F.rows(), 3);
	Eigen::VectorXd F_areas(data_model.target_developable.F.rows());

	for (int fi = 0; fi < data_model.target_developable.F.rows(); fi++)
	{
		Eigen::RowVectorXi face = data_model.target_developable.F.row(fi);
		Eigen::RowVectorXd face_sum = Eigen::RowVectorXd::Zero(3);

		for (int vi = 0; vi < 3; vi++)
		{
			if (face(vi) >= 0) //if it's -1 then there is no face because it's the boundary
				face_sum = face_sum + data_model.target_developable.V.row(face(vi));
		}

		F_midpoints.row(fi) = (face_sum / 3).eval();

		Eigen::RowVector3d a = data_model.target_developable.V.row(face(1)) - data_model.target_developable.V.row(face(0));
		Eigen::RowVector3d b = data_model.target_developable.V.row(face(2)) - data_model.target_developable.V.row(face(0));
		F_areas(fi) = a.cross(b).norm() / 2;
		////write_log(0) << fi << ": F_areas.row = " << F_areas.row(fi) << ", s = " << s << endl;
	}

	//write_log(0) << "F_midpoints: " << endl << F_midpoints << endl;
	//write_log(0) << "F_areas: " << endl << F_areas << endl;

	//normalize area
	double max_area = F_areas.maxCoeff();
	F_areas /= max_area;
	//write_log(0) << "nomalized F_areas: " << endl << F_areas << endl;


	//TODO store target_developable as Mesh instead of only Mesh
	Eigen::MatrixXi adjacency_FF;
	igl::triangle_triangle_adjacency(data_model.target_developable.F, adjacency_FF);

	Eigen::VectorXi face_indices;
	vector<Eigen::VectorXd> patch_face_distances;
	//vector<Eigen::MatrixXd> patch_face_points;

	for (int i = 0; i < data_model.patches.size(); i++)
	{
		WrapperMesh& wrapper = data_model.patches[i]->wrapper;

		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd closest_points;
		igl::point_mesh_squared_distance(F_midpoints, wrapper.V, wrapper.F, distances_squared, face_indices, closest_points);

		distances_squared *= F_areas; //weigh by area

		patch_face_distances.push_back(distances_squared);
		//patch_vertex_points.push_back(closest_points);
	}



	//write_log(0) << endl << "F: " << endl << data_model.target_developable.F << endl;
	//write_log(0) << endl << "FF: " << endl << adjacency_FF << endl;

	Eigen::MatrixXi face_edges(adjacency_FF.rows() * 2, 2);
	int number_edges = 0;

	for (int fi = 0; fi < adjacency_FF.rows(); fi++)
	{
		Eigen::RowVectorXi all_neighbor_faces = adjacency_FF.row(fi);
		//write_log(0) << "  neighbor faces: " << all_neighbor_faces << endl;

		for (int n = 0; n < all_neighbor_faces.size(); n++)
		{
			int neighbor_fi = all_neighbor_faces[n];
			if (neighbor_fi < fi)
				continue;

			face_edges(number_edges, 0) = fi;
			face_edges(number_edges, 1) = neighbor_fi;
			number_edges++;
		}
	}

	face_edges.conservativeResize(number_edges, Eigen::NoChange);
	//write_log(0) << endl << "face_edges: " << endl << face_edges << endl;


	//vector<vector<int>> F_adjacency;
	//igl::adjacency_list(data_model.target_developable.F, F_adjacency);
	//write_log(0) << endl << "adjacency_list of F: " << endl;
	//for (int i = 0; i < F_adjacency.size(); i++)
	//	log_list(0, F_adjacency[i], to_string(i) + ": ", false);


	//Eigen::MatrixXi adjacency_mid_FF;
	//igl::edges(adjacency_FF, adjacency_mid_FF);
	//write_log(0) << endl << "adjacency_mid_FF: " << endl << adjacency_mid_FF << endl;

	std::vector<double> resulting_distances;
	std::vector<int> face_patch_assignment = graphcut::compute_graph_cut_face_labeling(F_midpoints, face_edges, patch_face_distances, data_model.label_selection_smoothness, resulting_distances);
	data_model.vertex_patch_assignment = face_patch_assignment;
	//write_log(4) << "assigned " << *max_element(data_model.vertex_patch_assignment.begin(), data_model.vertex_patch_assignment.end())+1 << " labels" << endl;

	//for (int i = 0; i < data_model.target_developable.V.rows(); i++)
	//{
	//	int patch_index = data_model.vertex_patch_assignment[i];
	//	data_model.target_developable.V.row(i) = patch_vertex_points[patch_index].row(i);
	//}

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

	//ofstream file2("RESULT_output_face_distances.txt");
	//file2 << "";
	//for (double distance : resulting_distances)
	//	file2 << distance << endl;
	//file2.close();


	//find_labeling_seams(data_model, face_patch_assignment, face_edges, adjacency_FF);
	//find_labeling_seams(data_model, face_patch_assignment);
}
*/

/*
void meshhelper::update_target_developable_vertex_optimized(DataModel& data_model)
{
	Eigen::VectorXi face_indices;
	vector<Eigen::VectorXd> patch_vertex_distances;
	vector<Eigen::MatrixXd> patch_vertex_points;

	for (int i = 0; i < data_model.patches.size(); i++)
	{
		WrapperMesh& wrapper = data_model.patches[i]->wrapper;

		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd closest_points;
		igl::point_mesh_squared_distance(data_model.target_developable.V, wrapper.V, wrapper.F, distances_squared, face_indices, closest_points);

		patch_vertex_distances.push_back(distances_squared);
		patch_vertex_points.push_back(closest_points);
	}

	data_model.vertex_patch_assignment = graphcut::compute_graph_cut_labeling(data_model.target_developable, patch_vertex_distances, data_model.label_selection_smoothness);
	vector<int> labels;

	for (int i = 0; i < data_model.target_developable.V.rows(); i++)
	{
		int patch_index = data_model.vertex_patch_assignment[i];
		data_model.target_developable.V.row(i) = patch_vertex_points[patch_index].row(i);

		if (!is_contained(patch_index, labels))
			labels.push_back(patch_index);
	}
	write_log(4) << "(vertex) assigned " << labels.size() << " labels: " << endl; log_list(4, labels, "", false);


	//find_result_creases(data_model);


	////reset developables parts
	//data_model.developable_parts.clear();

	//for (int i = 0; i < data_model.developable_parts_view_indices.size(); i++)
	//	data_model.viewer.data_list[data_model.developable_parts_view_indices[i]].clear();

	//data_model.developable_parts_view_indices.clear();
}
*/

//Eigen::RowVector3d get_color(int index, int number_colors)
//{
//	if (number_colors < 1)
//		return Eigen::RowVector3d(0,0,0);
//
//	double hue = index * 360/number_colors;
//	Eigen::RowVector3d hsv(hue, 1, 1);
//	Eigen::RowVector3d rgb;
//	igl::hsv_to_rgb(hsv, rgb);
//
//	return rgb;
//}

/*
//void find_labeling_seams(DataModel& data_model, const std::vector<int>& face_patch_assignment, Eigen::MatrixXi& face_edges, const Eigen::MatrixXi& adjacency_FF)
void find_labeling_seams(DataModel& data_model, const std::vector<int>& face_patch_assignment)
{
	if (data_model.patches.size() < 2)
		return;
	if (face_patch_assignment.size() < 1)
		return;

	const int loglevel = 4;
	//log_list(loglevel+1, face_patch_assignment, "face_patch_assignment: ", false);


	int number_patches = data_model.patches.size();

	vector<vector<int>> labeled_faces(number_patches);
	vector<vector<int>> labeled_vertices(number_patches);

	for (int fi = 0; fi < face_patch_assignment.size(); fi++)
	{
		int label = face_patch_assignment[fi];
		labeled_faces[label].push_back(fi);

		Eigen::RowVector3i v = data_model.target_developable.F.row(fi);
		vector<int> vertices(v.data(), v.data() + v.size());
		labeled_vertices[label].insert(labeled_vertices[label].end(), vertices.begin(), vertices.end());
	}

	//reset developables parts
	data_model.developable_parts.clear();

	for (int i = 0; i < data_model.developable_parts_view_indices.size(); i++)
		data_model.viewer.data_list[data_model.developable_parts_view_indices[i]].clear();

	data_model.developable_parts_view_indices.clear();


	vector<int> labels;

	//create developable parts
	for (int i = 0; i < number_patches; i++)
	{
		vector<int> vertices = labeled_vertices[i];
		if (vertices.size() <= 0)
			continue;

		labels.push_back(i);

		std::sort(vertices.begin(), vertices.end());
		vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
		labeled_vertices[i] = vertices;

		Mesh mesh = meshhelper::sub_mesh(labeled_vertices[i], labeled_faces[i], data_model.target_developable.V, data_model.target_developable.F);
		
		data_model.developable_parts.push_back(mesh);

		meshhelper::add_mesh(data_model.viewer, mesh.V, mesh.F, false, false);
		data_model.developable_parts_view_indices.push_back(data_model.viewer.selected_data_index);
		data_model.viewer.data().set_colors(get_color(i, number_patches));


		//write_log(0) << endl << "label " << i << ": " << endl;
		//log_list(0, labeled_faces[i], "faces: ", false);
		//log_list(0, labeled_vertices[i], "vertices: ", false);
	}

	write_log(loglevel) << "(face) assigned " << labels.size() << " labels: " << endl; log_list(loglevel, labels, "", false);


	//find seams between developable parts

	data_model.number_crease_pairs = 0;
	data_model.crease_indices_global.clear(); //global on target_developable
	data_model.patch_neighbor_crease_lookup.resize(number_patches, number_patches);
	data_model.patch_neighbor_crease_lookup.setConstant(-1);

	for (int i = 0; i < number_patches; i++)
	{
		if (labeled_vertices[i].size() <= 0)
			continue;

		for (int j = 0; j < number_patches; j++)
		{
			if (labeled_vertices[j].size() <= 0)
				continue;
			if (i <= j)
				continue;

			vector<int> seam_vertices(labeled_vertices[i].size() + labeled_vertices[j].size());
			auto iterator = std::set_intersection(labeled_vertices[i].begin(), labeled_vertices[i].end(), labeled_vertices[j].begin(), labeled_vertices[j].end(), seam_vertices.begin());
			seam_vertices.resize(iterator - seam_vertices.begin());
			write_log(loglevel) << "found " << seam_vertices.size() << " seam vertices between labels " << i << " & " << j << endl;

			if (seam_vertices.size() < 1)
				continue;

			data_model.patch_neighbor_crease_lookup(j, i) = data_model.number_crease_pairs;
			data_model.crease_indices_global.push_back(seam_vertices);
			data_model.number_crease_pairs++;



			//store local vertex pairs for stiching
			for (int k = 0; k < seam_vertices.size(); k++)
			{
				int vertex_global_index = seam_vertices[k];
				Eigen::RowVector3d vertex_global = data_model.target_developable.V.row(vertex_global_index);

				Mesh part1 = data_model.developable_parts[i];
				Mesh part2 = data_model.developable_parts[j];

				int pair_index1 = get_closest_vertex(part1.V, vertex_global);
				int pair_index2 = get_closest_vertex(part2.V, vertex_global);

				data_model.crease_indices_local.push_back({ vertex_global_index, i, pair_index1, j, pair_index2 });
			}
		}
	}
	write_log(loglevel) << "data_model.patch_neighbor_crease_lookup: " << endl << data_model.patch_neighbor_crease_lookup << endl;
	for (int i = 0; i < data_model.crease_indices_global.size(); i++)
	{
		write_log(loglevel) << i; 
		log_list(loglevel, data_model.crease_indices_global[i], ": ", false);
	}


}
*/

/*
void find_result_creases(DataModel& data_model)
{
	if (data_model.vertex_patch_assignment.size() < 1)
		return;
	if (data_model.patches.size() < 2)
		return;

	const int loglevel = 5;

	//vector<vector<int>> patch_vertex_assignment(data_model.patches.size());

	//for (int i = 0; i < data_model.vertex_patch_assignment.size(); i++)
	//{
	//	int patch_index = data_model.vertex_patch_assignment[i];
	//	patch_vertex_assignment[patch_index].push_back(i);
	//}

	//vector<Mesh> meshes = meshhelper::cut_mesh(patch_vertex_assignment, data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.adjacency_VF);
	//
	//vector<int> boundary;
	//igl::boundary_loop(meshes[0].F, boundary);

	int number_patches = data_model.patches.size();
	//int number_pairs = -1;
	data_model.number_crease_pairs = 0;

	vector<vector<int>> adjacency_VV;
	igl::adjacency_list(data_model.target_developable.F, adjacency_VV);

	//Eigen::MatrixXi patch_neighbor_lookup(number_patches, number_patches);
	data_model.patch_neighbor_crease_lookup.resize(number_patches, number_patches);
	data_model.patch_neighbor_crease_lookup.setConstant(-1);

	//vector<vector<int>> indices = vector<vector<int>>(number_patches*number_patches*2);
	data_model.crease_indices_global.clear();
	data_model.crease_indices_global = vector<vector<int>>(number_patches*number_patches*2);

	for (int vertex_index = 0; vertex_index < data_model.vertex_patch_assignment.size(); vertex_index++)
	{
		int patch_index = data_model.vertex_patch_assignment[vertex_index];
		write_log(loglevel) << "CURRENT vertex_index = " << vertex_index << ", patch_index = " << patch_index << endl;


		vector<int> all_neighbor_indices = adjacency_VV[vertex_index];
		vector<int> unique_neighbor_patches;
		vector<int> unique_neighbor_vertices;
		
		for (int n = 0; n < all_neighbor_indices.size(); n++)
		{
			int neighbor_vertex_index = all_neighbor_indices[n];
			int neighbor_patch_index = data_model.vertex_patch_assignment[neighbor_vertex_index];

			if (neighbor_patch_index == patch_index)
				continue;
			//if (is_contained(neighbor_patch_index, unique_neighbor_patches) || neighbor_patch_index == patch_index)
			//	continue;

			unique_neighbor_patches.push_back(neighbor_patch_index);
			unique_neighbor_vertices.push_back(neighbor_vertex_index);
		}

		if (unique_neighbor_patches.size() < 1)
			continue;

		log_list(loglevel, unique_neighbor_patches, "unique_neighbor_patches: ", false);
		log_list(loglevel, unique_neighbor_vertices, "unique_neighbor_vertices: ", false);

		for (int j = 0; j < unique_neighbor_vertices.size(); j++)
		{
			int v = unique_neighbor_vertices[j];
			int p = unique_neighbor_patches[j];

			int index_lookup = patch_index < p ? data_model.patch_neighbor_crease_lookup(patch_index, p) : data_model.patch_neighbor_crease_lookup(p, patch_index);

			if (index_lookup == -1)
			{
				int row = min({ patch_index, p });
				int col = max({ patch_index, p });

				data_model.patch_neighbor_crease_lookup(row, col) = data_model.number_crease_pairs;
				index_lookup = data_model.number_crease_pairs;
				data_model.number_crease_pairs++;
			}

			int v1 = patch_index < p ? vertex_index : v;
			int v2 = patch_index < p ? v : vertex_index;

			data_model.crease_indices_global[index_lookup * 2].push_back(v1);
			data_model.crease_indices_global[index_lookup * 2 + 1].push_back(v2);
		}

		write_log(loglevel) << "lookup: " << endl << data_model.patch_neighbor_crease_lookup << endl << endl;
		for (int r = 0; r < number_patches; r++)
		{
			for (int c = 0; c < number_patches; c++)
			{
				int index_lookup = data_model.patch_neighbor_crease_lookup(r, c);
				if (index_lookup == -1)
					continue;

				log_list(loglevel, data_model.crease_indices_global[index_lookup * 2], "indices[" + to_string(index_lookup * 2) + "]: ", false);
				log_list(loglevel, data_model.crease_indices_global[index_lookup * 2 + 1], "indices[" + to_string(index_lookup * 2 + 1) + "]: ", false);
			}
		}

		write_log(loglevel) << endl;
	}

	
	write_log(loglevel-1) << endl << "DONE finding crease pairings. found " << data_model.number_crease_pairs << " pairs." << endl << endl;
	write_log(loglevel-1) << "lookup: " << endl << data_model.patch_neighbor_crease_lookup << endl << endl;
	for (int r = 0; r < number_patches; r++)
	{
		for (int c = 0; c < number_patches; c++)
		{
			int index_lookup = data_model.patch_neighbor_crease_lookup(r, c);
			if (index_lookup == -1)
				continue;

			log_list(loglevel, data_model.crease_indices_global[index_lookup * 2], "indices[" + to_string(index_lookup * 2) + "]: ", false);
			log_list(loglevel, data_model.crease_indices_global[index_lookup * 2 + 1], "indices[" + to_string(index_lookup * 2 + 1) + "]: ", false);
		}
	}

	//write_log(loglevel) << endl << "vertices for mathematica" << endl << endl;
	//for (int r = 0; r < data_model.number_crease_pairs; r++)
	//{
	//	Eigen::MatrixXd vertices;
	//	index_to_value(data_model.crease_indices[r], data_model.target_developable.V, vertices);
	//	write_log(loglevel) << "vertices[" << r << "]: "; print_mathematica_matrix(vertices);
	//}

	//for (int r = 0; r < number_patches; r++)
	//{
		//for (int c = 0; c < number_patches; c++)
		//{
		//	int index_lookup = patch_neighbor_lookup(r, c);
		//	if (index_lookup == -1)
		//		continue;
		//
		//	Eigen::MatrixXd vertices;
		//	index_to_value(indices[index_lookup * 2], data_model.target_developable.V, vertices);
		//	write_log(loglevel-1) << "vertices[" + to_string(index_lookup * 2) + "]: "; print_mathematica_matrix(vertices);
		//
		//	index_to_value(indices[index_lookup * 2 + 1], data_model.target_developable.V, vertices);
		//	write_log(loglevel-1) << "vertices[" + to_string(index_lookup * 2 + 1) + "]: "; print_mathematica_matrix(vertices);
		//}
	//}

	write_log(loglevel) << endl;

}
*/

/*
void meshhelper::update_target_developable_initialization(DataModel& data_model)
{
	if (data_model.selected_ruled_vertex_indices.size() < 1)
		return;

	igl::upsample(data_model.target_developable.V, data_model.target_developable.F, data_model.target_developable.V, data_model.target_developable.F);


	//TODO fix! don't use from view, but from data_model.ruled_developables
	std::vector<Mesh> selected_ruled_developables;
	for (int view_index : data_model.selected_ruled_view_indices)
	{
		Mesh developable;
		developable.V = data_model.viewer.data_list[view_index].V;
		developable.F = data_model.viewer.data_list[view_index].F;
		selected_ruled_developables.push_back(developable);
	}

	Mesh concatenated_surfaces = concatenate_meshes(selected_ruled_developables);
	//Mesh concatenated_surfaces = concatenate_ruled_developables(data_model.ruled_developables);

	Eigen::VectorXi face_indices_unused;
	Eigen::VectorXd distances_squared;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(data_model.target_developable.V, concatenated_surfaces.V, concatenated_surfaces.F, distances_squared, face_indices_unused, closest_points);

	//move vertices to wrappers
	data_model.target_developable.V = closest_points;
}
*/

void meshhelper::calculate_dimensions(const Eigen::MatrixXd & vertices, Eigen::RowVectorXd & out_dimensions)
{
	const Eigen::RowVectorXd min_point = vertices.colwise().minCoeff();
	const Eigen::RowVectorXd max_point = vertices.colwise().maxCoeff();
	out_dimensions = max_point - min_point;
}

/*
void meshhelper::update_covering_wrappers(DataModel& data_model, Eigen::MatrixXi& out_concatenated_covering_F)
{
	if (data_model.patches.size() < 1)
		return;

	//concatenate wrappers
	meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

	//update covering faces (which are set as render faces)
	
	//find closest points from target vertices to all wrapper patches
	Eigen::VectorXi face_indices;
	Eigen::VectorXd distances_squared;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(data_model.target.V, data_model.concatenated_patches_V, data_model.concatenated_patches_F, distances_squared, face_indices, closest_points);

	//keep only distances that are under coverage threshold
	int number_covering_faces = 0;
	Eigen::VectorXi covering_face_indices(distances_squared.rows());
	Eigen::VectorXd covering_distances_squared(distances_squared.rows());
	Eigen::MatrixXd covering_closest_points(distances_squared.rows(), closest_points.cols());
	
	for (int i = 0; i < distances_squared.rows(); i++)
	{
		if (distances_squared(i) > data_model.optimization_settings.coverage_threshold)
			continue;

		covering_face_indices(number_covering_faces) = face_indices(i);
		covering_distances_squared(number_covering_faces) = covering_distances_squared(i);
		covering_closest_points.row(number_covering_faces) = covering_closest_points.row(i);
		
		number_covering_faces++;
	}

	covering_face_indices.conservativeResize(number_covering_faces);
	covering_distances_squared.conservativeResize(number_covering_faces);
	covering_closest_points.conservativeResize(number_covering_faces, Eigen::NoChange);

	data_model.target.covering_wrapper_points.resize(0, Eigen::NoChange);
	data_model.target.covering_wrapper_points = closest_points;
	//data_model.target.covering_wrapper_points = covering_closest_points;

		//not sure if i need this
	data_model.target.covering_wrapper_faces_distances.resize(0);
	data_model.target.covering_wrapper_faces_distances = distances_squared;
	//data_model.target.covering_wrapper_faces_distances = covering_distances_squared;

	//split concatenated render F into patches for visualization, while I'm at it


	index_to_value(covering_face_indices, data_model.concatenated_patches_F, out_concatenated_covering_F);
}
*/

/*
void meshhelper::update_target_coverage(DataModel& data_model, const Eigen::MatrixXi& concatenated_covering_F)
{
	//find closest point on _covering_ faces (as defined by update_covering_wrappers()) from all target vertices
	Eigen::VectorXi face_indices_unused;
	Eigen::VectorXd distances_squared;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(data_model.target.V, data_model.concatenated_patches_V, concatenated_covering_F, distances_squared, face_indices_unused, closest_points);

	//set coverage
	const double max_distance = data_model.optimization_settings.coverage_threshold;
	const int number_points = data_model.target.V.rows();

	int number_uncovered = 0;
	data_model.target.uncovered_V.resize(number_points, 3);

	int number_covered = 0;
	data_model.target.covered_V.resize(number_points, 3);

	for (int i = 0; i < number_points; i++)
	{
		if (distances_squared(i) > max_distance)
		{
			data_model.target.uncovered_V(number_uncovered) = i;
			number_uncovered++;
		}
		else
		{
			data_model.target.covered_V(number_covered) = i;
			number_covered++;
		}
	}

	data_model.target.uncovered_V.conservativeResize(number_uncovered, Eigen::NoChange);
	data_model.target.covered_V.conservativeResize(number_covered, Eigen::NoChange);
	data_model.target.coverage_distances_V = distances_squared;

	//write_log(loglevel) << "check_coverage() sqrD:" << endl << out_distances_squared << endl;
	////write_log(loglevel) << "check_coverage() face_indices:" << endl << out_face_indices << endl;
	////write_log(0) << "check_coverage() all_wrapper #F:" << endl << all_wrapper_F.rows() << endl;
	//write_log(loglevel) << "check_coverage() closest_points:" << endl << out_closest_points << endl << endl;

	//write_log(loglevel) << "check_coverage() uncovered:" << endl << out_uncovered << endl << endl;
	//write_log(loglevel) << "check_coverage() covered:" << endl << out_covered << endl << endl;
}
*/

/*
void meshhelper::update_coverage_per_vertices(DataModel& data_model)
{
	//concatenate patches
	concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

	meshhelper::get_target_coverage(data_model, data_model.target.V, data_model.target.covered_V, data_model.target.uncovered_V, data_model.target.coverage_distances_V, data_model.target.covering_wrapper_points, data_model.target.covering_wrapper_faces_distances);
              //get_target_coverage(data_model, source_points,       out_covered,                 out_uncovered,                 out_distances_squared,                  out_closest_points,                ...)
}
*/

/*
void meshhelper::get_target_coverage(const DataModel& data_model, const Eigen::MatrixXd& source_points, Eigen::VectorXi& out_covered, Eigen::VectorXi& out_uncovered, Eigen::VectorXd& out_distances_squared, Eigen::MatrixXd& out_closest_points, Eigen::VectorXd& out_covering_wrapper_faces_distances)
{
	const int loglevel = 6;
	const int number_points = source_points.rows();

	if (data_model.patches.size() < 1)
	{
		out_covered.resize(0);

		out_uncovered.resize(number_points);
		out_uncovered.setLinSpaced(0, number_points - 1);

		out_distances_squared.resize(number_points);
		out_distances_squared.setConstant(1e6);

		return;
	}

	//find closest points from target vertices to all wrapper patches
	Eigen::VectorXi face_indices_unused;
	igl::point_mesh_squared_distance(source_points, data_model.concatenated_patches_V, data_model.concatenated_patches_F, out_distances_squared, face_indices_unused, out_closest_points);

	//find distances from each wrapper face (triangle) to the target
	int cols = data_model.concatenated_patches_F.cols();
	Eigen::MatrixXd mid_faces(data_model.concatenated_patches_F.rows(), cols);
	for (int i = 0; i < data_model.concatenated_patches_F.rows(); i++)
	{
		Eigen::RowVector3d point(0,0,0);

		for (int j = 0; j < cols; j++)
			point += data_model.concatenated_patches_V.row(data_model.concatenated_patches_F(i, j));

		mid_faces.row(i) = point/cols;
	}

	Eigen::MatrixXd points_unused;
	igl::point_mesh_squared_distance(mid_faces, data_model.target.V, data_model.target.F, out_covering_wrapper_faces_distances, face_indices_unused, points_unused);


	//visualize uncovered faces
	//const double max_distance = data_model.outlier_threshold * data_model.outlier_threshold;
	const double max_distance = data_model.optimization_settings.coverage_threshold;

	int number_uncovered = 0;
	out_uncovered.resize(number_points, 3);

	int number_covered = 0;
	out_covered.resize(number_points, 3);

	for (int i = 0; i < number_points; i++)
	{
		if (out_distances_squared(i) > max_distance)
		{
			out_uncovered(number_uncovered) = i;
			number_uncovered++;
		}
		else
		{
			out_covered(number_covered) = i;
			number_covered++;
		}
	}

	out_uncovered.conservativeResize(number_uncovered, Eigen::NoChange);
	out_covered.conservativeResize(number_covered, Eigen::NoChange);

	write_log(loglevel) << "check_coverage() sqrD:" << endl << out_distances_squared << endl;
	//write_log(loglevel) << "check_coverage() face_indices:" << endl << out_face_indices << endl;
	//write_log(0) << "check_coverage() all_wrapper #F:" << endl << all_wrapper_F.rows() << endl;
	write_log(loglevel) << "check_coverage() closest_points:" << endl << out_closest_points << endl << endl;

	write_log(loglevel) << "check_coverage() uncovered:" << endl << out_uncovered << endl << endl;
	write_log(loglevel) << "check_coverage() covered:" << endl << out_covered << endl << endl;


	//int furthest_vertex_index;
	//out_distances_squared.maxCoeff(&furthest_vertex_index);
	////data_model.debug_points = mid_faces.row(furthest_face_index);

	//int vertex_index;
	//Eigen::RowVectorXd vertex;
	//get_closest_vertex(data_model.target.V, data_model.target.V.row(furthest_vertex_index), vertex_index, vertex);
	//debug_show_point(vertex, data_model);


	////debug visualize
	//append_matrix(source_points, data_model.debug_line_vertices1);
	//append_matrix(out_closest_points, data_model.debug_line_vertices2);
}
*/

//void meshhelper::get_coverage_per_faces(DataModel& data_model, Eigen::VectorXi& out_covered, Eigen::VectorXi& out_uncovered, Eigen::VectorXd& out_distances_squared)
//{
//	//create list of face mid points
//	const int cols = data_model.target.F.cols();
//	Eigen::MatrixXd mid_faces(data_model.target.F.rows(), cols);
//
//	for (int i = 0; i < data_model.target.F.rows(); i++)
//	{
//		Eigen::RowVector3d point(0,0,0);
//
//		for (int j = 0; j < cols; j++)
//			point += data_model.target.V.row(data_model.target.F(i, j));
//	
//		mid_faces.row(i) = point/cols;
//	}
//	write_log(6) << "check_coverage() mid_faces:" << endl << mid_faces << endl << endl;
//
//	get_target_coverage(data_model, mid_faces, out_covered, out_uncovered, out_distances_squared);
//
//	//debug visualize
//	data_model.debug_faces_target = out_uncovered;
//}
//
//void meshhelper::get_coverage_per_vertices(DataModel& data_model, Eigen::VectorXi& out_covered, Eigen::VectorXi& out_uncovered, Eigen::VectorXd& out_distances_squared)
//{
//	get_target_coverage(data_model, data_model.target.V, out_covered, out_uncovered, out_distances_squared);
//
//	//debug visualize
//	data_model.debug_points.resize(out_uncovered.rows(), 3);
//	for (int i = 0; i < out_uncovered.rows(); i++)
//		data_model.debug_points.row(i) = data_model.target.V.row(out_uncovered(i));
//}

void meshhelper::add_mesh(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool do_center_mesh, bool align_camera)
{
	if (do_center_mesh)
		center_mesh(V);

	viewer.append_mesh();
	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	//viewer.data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
	viewer.data().compute_normals();

	viewer.data().point_size = 6.0;
	viewer.data().line_width = 1.0;

	if (align_camera)
		viewer.core().align_camera_center(V, F);
}

//void print_mathematica_m(const Eigen::MatrixXd matrix)
//{
//	int rows = matrix.rows();
//	int cols = matrix.cols();
//
//	if (rows < 1 || cols < 1)
//		return;
//
//	std::cout << "{";
//	for (int r = 0; r < rows; r++)
//	{
//		if (cols == 1)
//		{
//			std::cout << std::fixed << matrix(r, 0);
//		}
//		else
//		{
//			if (rows > 1)
//				std::cout << "{";
//
//			for (int c = 0; c < cols; c++)
//			{
//				std::cout << std::fixed << matrix(r, c);
//				if (c < cols - 1)
//					std::cout << ", ";
//			}
//			if (rows > 1)
//				std::cout << "}";
//		}
//
//		if (r < rows - 1)
//			std::cout << "," << std::endl;
//	}
//	std::cout << "};" << std::endl;
//}


void center_mesh(Eigen::MatrixXd& V)
{
	auto centroid = get_centroid(V);
	Eigen::Vector3d shift;
	shift.setConstant(0);
	shift.head(centroid.size()) = -centroid.cast<double>();

	V = V.rowwise() + shift.transpose();
	write_log(4) << "translate = " << to_mathematica(shift.transpose()) << endl;
}

//void plot_gaussian_curvature(DataModel& data_model, const Eigen::VectorXd& K)
//{
//	if (!K.rows())
//		return;
//	if (data_model.target_view_index == -1)
//		return;
//
//	Eigen::MatrixXd C;
//	igl::jet(K, 0.0, 10.0, C); // Compute pseudocolor
//
//	// Plot the mesh with pseudocolors
//	data_model.viewer.data_list[data_model.target_view_index].set_colors(C); //REFACTOR TO debug view
//}


Mesh meshhelper::concatenate_meshes(const std::vector<Eigen::MatrixXd>& list_V, const std::vector<Eigen::MatrixXi>& list_F)
{
	if (list_V.size() != list_F.size())
		exit(EXIT_FAILURE);

	Mesh concatenated;

	for (int i = 0; i < list_V.size(); i++)
	{
		auto V = list_V[i];
		auto F = list_F[i];

		Eigen::RowVector3i offset(concatenated.V.rows(), concatenated.V.rows(), concatenated.V.rows());
		Eigen::MatrixXi offset_faces = F.rowwise() + offset;

		append_matrix(offset_faces, concatenated.F);
		append_matrix(V, concatenated.V);
	}

	return concatenated;
}

Mesh meshhelper::concatenate_meshes(const std::vector<Mesh>& meshes)
{
	Mesh concatenated;
	vector<int> cumulative_offset_vertices;
	vector<int> cumulative_offset_faces;
	concatenate_meshes(meshes, concatenated, cumulative_offset_vertices, cumulative_offset_faces);

	return concatenated;
	/*
	Mesh concatenated;

	for (auto surface : meshes)
	{
		Eigen::RowVector3i offset(concatenated.V.rows(), concatenated.V.rows(), concatenated.V.rows());
		Eigen::MatrixXi offset_faces = surface.F.rowwise() + offset;

		append_matrix(offset_faces, concatenated.F);
		append_matrix(surface.V, concatenated.V);
	}

	return concatenated;
	*/
}

std::vector<Mesh> meshhelper::get_wrapper_meshes(const std::vector<Patch*>& wrappers)
{
	std::vector<Mesh> meshes;
	for (auto wrapper : wrappers)
		meshes.push_back(wrapper->wrapper);

	return meshes;
}

std::vector<Mesh> meshhelper::get_selected_ruled_developable_meshes(const std::vector<int>& selected_ruled_vertex_indices, const std::vector<int>& ruled_vertex_indices, const std::vector<RuledDevelopableSurface*> &ruled_developables)
{
	std::vector<Mesh> meshes;
	for (auto selected_vi : selected_ruled_vertex_indices)
	{
		int i = index_of(selected_vi, ruled_vertex_indices);
		meshes.push_back(ruled_developables[i]->developable);
	}

	return meshes;
}

std::vector<Mesh> meshhelper::get_ruled_developable_meshes(const std::vector<RuledDevelopableSurface*>& ruled_developables)
{
	std::vector<Mesh> meshes;
	for (auto ruled : ruled_developables)
		meshes.push_back(ruled->developable);

	return meshes;
}

void meshhelper::concatenate_meshes(const std::vector<Mesh>& meshes, Mesh& out_concatenated, std::vector<int>& out_cumulative_offset_vertices, std::vector<int>& out_cumulative_offset_faces)
{
	out_concatenated.V.resize(0, 3);
	out_concatenated.F.resize(0, 3);
	out_cumulative_offset_vertices.clear();
	out_cumulative_offset_faces.clear();

	for (Mesh mesh : meshes)
	{
		int number_vertices = out_concatenated.V.rows();
		out_cumulative_offset_vertices.push_back(number_vertices);
		out_cumulative_offset_faces.push_back(out_concatenated.F.rows());

		Eigen::RowVector3i offset(number_vertices, number_vertices, number_vertices);
		Eigen::MatrixXi offset_faces = mesh.F.rowwise() + offset;

		append_matrix(offset_faces, out_concatenated.F);
		append_matrix(mesh.V, out_concatenated.V);
	}
}

Mesh meshhelper::concatenate_ruled_developables(const std::vector<RuledDevelopableSurface*>& ruled_developables)
{
	//std::vector<Mesh> meshes;
	//for (auto ruled : ruled_developables)
	//	meshes.push_back(ruled->developable);

	//return concatenate_meshes(meshes);

	return concatenate_meshes(meshhelper::get_ruled_developable_meshes(ruled_developables));
}

Mesh meshhelper::concatenate_wrappers(const std::vector<Patch*>& wrappers)
{
	//std::vector<Mesh> meshes;
	//for (auto wrapper : wrappers)
	//	meshes.push_back(wrapper->wrapper);

	//return concatenate_meshes(meshes);

	return concatenate_meshes(meshhelper::get_wrapper_meshes(wrappers));
}

void meshhelper::concatenate_wrappers(const DataModel& data_model, Eigen::MatrixXd& out_concatenated_wrappers_V, Eigen::MatrixXi& out_concatenated_wrappers_F)
{
	if (data_model.patches.size() < 1)
		return;

	Mesh mesh = concatenate_wrappers(data_model.patches);
	out_concatenated_wrappers_V = mesh.V;
	out_concatenated_wrappers_F = mesh.F;

	/*
	out_concatenated_wrappers_V.resize(0, data_model.patches[0]->wrapper.V.cols());
	out_concatenated_wrappers_F.resize(0, data_model.patches[0]->wrapper.F.cols());
	
	for (int i = 0; i < data_model.patches.size(); i++)
	{
		Eigen::RowVector3i offset(out_concatenated_wrappers_V.rows(), out_concatenated_wrappers_V.rows(), out_concatenated_wrappers_V.rows());
		Eigen::MatrixXi offset_faces = data_model.patches[i]->wrapper.F.rowwise() + offset;
		append_matrix(offset_faces, out_concatenated_wrappers_F);

		append_matrix(data_model.patches[i]->wrapper.V, out_concatenated_wrappers_V);
	}
	*/
	write_log(6) << "concatenated_wrappers_V:" << endl << out_concatenated_wrappers_V << endl;
	write_log(6) << "concatenated_wrappers_F:" << endl << out_concatenated_wrappers_F << endl << endl;
}

std::vector<int> meshhelper::boundary_loop_flat(const Eigen::MatrixXi& F)
{
	vector<vector<int>> boundary_indices;
	igl::boundary_loop(F, boundary_indices);

	vector<int> flat_boundary_indices;
	for (vector<int> loop : boundary_indices)
		flat_boundary_indices.insert(flat_boundary_indices.end(), loop.begin(), loop.end());

	//for (vector<int> loop : global_boundary_indices)
	//	log_list(5, loop, "loop: ", false);

	return flat_boundary_indices;
}

Eigen::MatrixXi meshhelper::face_edges(const Eigen::MatrixXi& F)
{
	Eigen::MatrixXi adjacency_FF;
	igl::triangle_triangle_adjacency(F, adjacency_FF);

	Eigen::MatrixXi F_edges(adjacency_FF.rows() * 2, 2);
	int number_edges = 0;

	for (int fi = 0; fi < adjacency_FF.rows(); fi++)
	{
		Eigen::RowVectorXi all_neighbor_faces = adjacency_FF.row(fi);
		//write_log(0) << "  neighbor faces: " << all_neighbor_faces << endl;

		for (int n = 0; n < all_neighbor_faces.size(); n++)
		{
			int neighbor_fi = all_neighbor_faces[n];
			if (neighbor_fi < fi)
				continue;

			F_edges(number_edges, 0) = fi;
			F_edges(number_edges, 1) = neighbor_fi;
			number_edges++;
		}
	}
	F_edges.conservativeResize(number_edges, Eigen::NoChange);

	return F_edges;
}


void meshhelper::project_to_meshes(const Eigen::MatrixXd& source_points, const std::vector<Mesh>& meshes, std::vector<Eigen::VectorXd>& out_patch_point_distances, std::vector<Eigen::MatrixXd>& out_patch_closest_points)
{
	Eigen::VectorXd weights(source_points.rows());
	weights.setConstant(1.0);
	project_to_meshes(source_points, weights, meshes, out_patch_point_distances, out_patch_closest_points);
}

void meshhelper::project_to_meshes(const Eigen::MatrixXd& source_points, const Eigen::VectorXd& weights, const std::vector<Mesh>& meshes, std::vector<Eigen::VectorXd>& out_patch_point_distances, std::vector<Eigen::MatrixXd>& out_patch_closest_points)
{
	Eigen::VectorXi face_indices_unused;
	out_patch_point_distances.clear();
	out_patch_closest_points.clear();

	for (int i = 0; i < meshes.size(); i++)
	{
		Mesh mesh = meshes[i];

		Eigen::VectorXd distances_squared;
		Eigen::MatrixXd closest_points;
		igl::point_mesh_squared_distance(source_points, mesh.V, mesh.F, distances_squared, face_indices_unused, closest_points);

		distances_squared *= weights; //e.g., weigh by area
		out_patch_point_distances.push_back(distances_squared);
		out_patch_closest_points.push_back(closest_points);
	}
}

void meshhelper::project_to_mesh(const Eigen::MatrixXd& source_points, const Mesh& mesh, Eigen::VectorXd& out_point_distances, Eigen::MatrixXd& out_closest_points)
{
	Eigen::VectorXi face_indices_unused;
	Eigen::VectorXd distances_squared;
	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(source_points, mesh.V, mesh.F, out_point_distances, face_indices_unused, out_closest_points);

	//out_point_distances *= weights; //e.g., weigh by area
}

void meshhelper::project_to_mesh(Eigen::MatrixXd& points, const Mesh& projection_target)
{
	Eigen::VectorXi face_indices_unused;
	Eigen::VectorXd distances_squared_unused;

	Eigen::MatrixXd closest_points;
	igl::point_mesh_squared_distance(points, projection_target.V, projection_target.F, distances_squared_unused, face_indices_unused, closest_points);
	
	points = closest_points;
}

void meshhelper::simple_project_to_DOGs(const std::vector<Patch*>& surfaces, Mesh& to_deform)
{
	auto concatenated_wrappers = meshhelper::concatenate_wrappers(surfaces);
	project_to_mesh(to_deform.V, concatenated_wrappers);
}

void meshhelper::simple_project_to_ruled_surfaces(const std::vector<RuledDevelopableSurface*>& surfaces, Mesh& to_deform)
{
	auto concatenated_surfaces = meshhelper::concatenate_ruled_developables(surfaces);
	project_to_mesh(to_deform.V, concatenated_surfaces);
}


/*
void surface::compute_weighted_features(Mesh& target, double weight_K, double weight_crease, double weight_H)
{
	target.surface_features.weighted_curvature = weight_K * target.surface_features.K_abs_normalized + weight_H * target.surface_features.H_normalized + weight_crease * target.surface_features.crease_vertices_normalized;
	create_lookup(target.surface_features.weighted_curvature, target.surface_features.weighted_curvature_lookup);
}

void surface::compute_features(Mesh& target)
{
	compute_principal_k(target);
	compute_K(target);
	compute_H(target);
	compute_creases(target);
}

void surface::compute_principal_k(Mesh& target)
{
	igl::principal_curvature(target.V, target.F, target.surface_features.principal_k_min_direction, target.surface_features.principal_k_max_direction, target.surface_features.principal_k_min, target.surface_features.principal_k_max);
}

void surface::compute_K(Mesh& target)
{
	get_gaussian_curvature(target, target.surface_features.K);
	target.surface_features.K_abs = target.surface_features.K.cwiseAbs();
	create_lookup(target.surface_features.K_abs, target.surface_features.K_abs_lookup);
	normalize_vector_entries(target.surface_features.K_abs, target.surface_features.K_abs_normalized);
}

void surface::compute_H(Mesh& target)
{
	get_mean_curvature(target, target.surface_features.H);
	create_lookup(target.surface_features.H, target.surface_features.H_lookup);
	normalize_vector_entries(target.surface_features.H, target.surface_features.H_normalized);
}

void surface::compute_creases(Mesh& target)
{
	crease::creases_ring_neighbors(target.V, target.F, target.adjacency_VF, 1, target.surface_features.crease_vertices);
	create_lookup(target.surface_features.crease_vertices, target.surface_features.crease_vertices_lookup);
	normalize_vector_entries(target.surface_features.crease_vertices, target.surface_features.crease_vertices_normalized);

	//write_log(0) << "sorted creases: " << endl;
	//for(int i = 0; i < target.surface_features.crease_vertices.rows(); i++)
	//	write_log(0) << target.surface_features.crease_vertices(target.surface_features.crease_vertices_lookup(i)) << endl;
}

void surface::normalize_vector_entries(const Eigen::VectorXd& curvature, Eigen::VectorXd& out_normalized_curvature)
{
	double min = curvature.minCoeff();
	double max = curvature.maxCoeff();
	double range = max - min;
	out_normalized_curvature = (curvature - min * Eigen::VectorXd::Ones(curvature.rows())) / range;
}

void surface::get_gaussian_curvature(Mesh& target, Eigen::VectorXd& out_K)
{
	//get gaussian curvature
	igl::gaussian_curvature(target.V, target.F, out_K); // Compute integral of Gaussian curvature

	// Compute mass matrix
	Eigen::SparseMatrix<double> M, Minv;
	igl::massmatrix(target.V, target.F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);

	out_K = (Minv*out_K).eval(); // Divide by area to get integral average
}

void surface::get_mean_curvature(Mesh& target, Eigen::VectorXd& out_H)
{
	Eigen::MatrixXd HN;
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(target.V, target.F, L);
	igl::massmatrix(target.V, target.F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);

	HN = -Minv * (L * target.V);
	out_H = HN.rowwise().norm(); //up to sign
}
*/