#include "MeshViewer.h"

#include <igl/avg_edge_length.h>
#include <igl/hsv_to_rgb.h>

#include "MeshController.h"
#include "Utils.h"
#include "Logger.h"


void red_blue_colormap(const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors);
void two_color_colormap(const Eigen::RowVector3d& color_min, const Eigen::RowVector3d& color_max, const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors);


void MeshViewer::add_mesh(Mesh& mesh)
{
	meshhelper::add_mesh(viewer, mesh.V, mesh.F);

	//settings.mesh = mesh;
	mesh_settings.cached_V = mesh.V;
	mesh_settings.view_index = viewer.selected_data_index;
	//meshhelper::calculate_dimensions(mesh.V, settings.dimensions);
}

void MeshViewer::update_view()
{
	if (mesh.V.rows() < 1)
		return;

	viewer.selected_data_index = mesh_settings.view_index;
	update_mesh();
}

void MeshViewer::update_menu()
{

}

void MeshViewer::update_debug_menu()
{

}

void MeshViewer::update_mesh()
{
	if (mesh.V.rows() < 1)
		return;

	if (mesh_settings.view_index == -1)
	{
		write_log(1) << linebreak << "ERROR in <" << __FUNCTION__ << "> target model is loaded but view index is not set. Bug -- go fix!" << linebreak << std::endl;
		exit(EXIT_FAILURE);
	}

	double difference = (mesh_settings.cached_V - mesh.V).sum();
	bool have_vertices_changed = abs(difference) < 1e-6;
	if (have_vertices_changed)
	{
		viewer.data().clear();
		viewer.data().set_mesh(mesh.V, mesh.F);
		mesh_settings.cached_V = mesh.V;
	}

	if (!mesh_settings.has_changed && !have_vertices_changed)
		return;

	update_mesh_view();
	update_mesh_fill();

	mesh_settings.has_changed = false;
}

void MeshViewer::update_mesh_view()
{
	viewer.data().show_faces = mesh_settings.show_faces && mesh_settings.show_mesh;
	viewer.data().show_lines = mesh_settings.show_wireframe && mesh_settings.show_mesh;

	viewer.data().show_vertid = mesh_settings.label_vertices && mesh_settings.show_mesh;
	viewer.data().show_faceid = mesh_settings.label_faces && mesh_settings.show_mesh;
}

void MeshViewer::update_mesh_fill()
{
	if (mesh_settings.view_index == -1)
		return;
	if (!mesh_settings.has_changed 	&& !overlay_settings.has_changed)
		return;
	if (mesh.V.rows() < 1)
		return;

	viewer.data().lines.resize(0, Eigen::NoChange);

	switch (mesh_settings.visualization)
	{
		case MeshVisualization::GaussianCurvature:
		{
			visualize_curvature(mesh.surface_features().K, mesh.surface_features().K_abs_lookup);
			break;
		}
		case MeshVisualization::MeanCurvature:
		{
			visualize_curvature(mesh.surface_features().H, mesh.surface_features().H_lookup);
			break;
		}
		case MeshVisualization::PrincipalCurvature:
		{
			visualize_principal_curvature();
			break;
		}
		case MeshVisualization::LabelAssignment:
		{
			set_assignment_fill();
			break;
		}
		case MeshVisualization::Standard:
		{
			set_standard_fill();
			break;
		}
		default:
		{
			set_standard_fill();
			break;
		}
	}

	mesh_settings.has_changed = false;
	overlay_settings.has_changed = false;
}

void MeshViewer::set_standard_fill()
{
	viewer.data().uniform_colors(
		Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
		Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
		Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

	viewer.data().grid_texture();
}

void MeshViewer::visualize_curvature(const Eigen::VectorXd& curvature, const Eigen::VectorXi& curvature_lookup)
{
	if (!curvature.rows())
		return;


	int end_index = overlay_settings.filter_curvature * (curvature.rows() - 1);
	overlay_settings.visualized_max = abs(curvature(curvature_lookup(end_index)));
	overlay_settings.visualized_min = overlay_settings.visualized_max * -1;

	int visualized_range_index = 0.3 * (curvature.rows() - 1); //clip top 30% of values to get rid of extremes and show a more consistent color range
	double visualized_range = abs(curvature(curvature_lookup(visualized_range_index)));

	Eigen::MatrixXd C;
	red_blue_colormap(overlay_settings.visualized_min, overlay_settings.visualized_max, visualized_range, curvature, C);
	viewer.data().set_colors(C);
}

void MeshViewer::visualize_principal_curvature()
{
	const double avg = igl::avg_edge_length(mesh.V, mesh.F);// *0.5;

	Eigen::MatrixXd N;
	igl::per_vertex_normals(mesh.V, mesh.F, N);

	const double offset = avg * 0.2;
	N *= offset;

	const int number_vertices = mesh.V.rows();

	Eigen::MatrixXd edges_min_start(number_vertices, 3);
	Eigen::MatrixXd edges_min_end(number_vertices, 3);
	Eigen::MatrixXd edges_max_start(number_vertices, 3);
	Eigen::MatrixXd edges_max_end(number_vertices, 3);

	for (int i = 0; i < number_vertices; i++)
	{
		double abs_min_k = abs(mesh.surface_features().principal_k_min(i));
		double abs_max_k = abs(mesh.surface_features().principal_k_max(i));

		double min_scale = abs_min_k > abs_max_k ? avg : avg * (abs_min_k / abs_max_k);
		double max_scale = abs_min_k > abs_max_k ? avg * (abs_max_k / abs_min_k) : avg;

		edges_min_start.row(i) = mesh.V.row(i) + N.row(i) + mesh.surface_features().principal_k_min_direction.row(i) * min_scale;
		edges_min_end.row(i) = mesh.V.row(i) + N.row(i) - mesh.surface_features().principal_k_min_direction.row(i) * min_scale;
		edges_max_start.row(i) = mesh.V.row(i) + N.row(i) + mesh.surface_features().principal_k_max_direction.row(i) * max_scale;
		edges_max_end.row(i) = mesh.V.row(i) + N.row(i) - mesh.surface_features().principal_k_max_direction.row(i) * max_scale;
	}

	viewer.data().add_edges(edges_min_start, edges_min_end, Colors::BLUE);
	viewer.data().add_edges(edges_max_start, edges_max_end, Colors::RED);
}

void MeshViewer::set_assignment_fill()
{
	if (!label_assignment)
		return;

	int max_label = *std::max_element(label_assignment->begin(), label_assignment->end());

	Eigen::MatrixXd colors;

	if (label_assignment->size() < 0)
	{
		colors.resize(1, 3);
		colors.row(0) = Eigen::RowVector3d(1.0, 1.0, 1.0);
	}
	else
	{
		const double hue_step = 360 / max_label;
		colors.resize(mesh.V.rows(), 3);

		for (int i = 0; i < mesh.V.rows(); i++)
		{
			int patch_index = (*label_assignment)[i];

			double hue = patch_index * hue_step;
			Eigen::RowVector3d hsv(hue, 1, 1);
			Eigen::RowVector3d rgb;
			igl::hsv_to_rgb(hsv, rgb);
			colors.row(i) = Colors::ligher(Colors::ligher(rgb));
		}
	}

	viewer.data().set_colors(colors);
}

void MeshViewer::show_filtered_curvature_points(const Eigen::VectorXd& feature, const Eigen::VectorXi& feature_lookup)
{
	if (!feature.rows())
		return;

	int end_index = overlay_settings.filter_curvature * (feature.rows() - 1);

	for (int i = 0; i <= end_index; i++)
	{
		int vertex_index = feature_lookup(i);

		if (overlay_settings.show_filtered_curvature_points)
			viewer.data().add_points(mesh.V.row(vertex_index), Eigen::RowVector3d(0.5, 0.5, 0.5));
		if (overlay_settings.show_filtered_curvature_values)
			viewer.data().add_label(mesh.V.row(vertex_index), std::to_string(feature(vertex_index)));
	}
}

void red_blue_colormap(const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors)
{
	two_color_colormap(Colors::BLUE, Colors::RED, value_min, value_max, value_range, values, out_colors);
}

void two_color_colormap(const Eigen::RowVector3d& color_min, const Eigen::RowVector3d& color_max, const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors)
{
	const int loglevel = 6;

	out_colors.resize(values.rows(), 3);
	Eigen::RowVector3d color_min_mask = (color_min - Eigen::RowVector3d::Ones()).cwiseAbs();
	Eigen::RowVector3d color_max_mask = (color_max - Eigen::RowVector3d::Ones()).cwiseAbs();

	write_log(loglevel) << std::endl;
	write_log(loglevel) << "values: " << values << std::endl;
	write_log(loglevel) << "color_min: " << color_min.transpose() << ", mask: " << color_min_mask.transpose() << std::endl;
	write_log(loglevel) << "color_max: " << color_max.transpose() << ", mask: " << color_max_mask.transpose() << std::endl;
	write_log(loglevel) << "value_min: " << value_min << ", value_max: " << value_max << ", denom: " << value_range << ", values.maxCoeff(): " << values.maxCoeff() << ", values.minCoeff(): " << values.minCoeff() << std::endl << std::endl;

	for (int i = 0; i < values.rows(); i++)
	{
		Eigen::Vector3d color = values(i) > 0 ? color_max : color_min;
		Eigen::Vector3d color_mask = values(i) > 0 ? color_max_mask : color_min_mask;
		double threshold = values(i) > 0 ? value_max : abs(value_min);

		double saturation = (abs(values(i)) - threshold) / value_range; //normalize to 1
		saturation = clip(saturation, 0.0, 1.0);
		saturation = 1 - saturation;

		out_colors.row(i) = color + color_mask * saturation;

		write_log(loglevel) << "[" << i << "]: " << values(i) << ", normalized: " << (values(i) - threshold) / value_range << ", added color: " << out_colors.row(i) << " (selected color: " << color.transpose() << ", mask: " << color_mask.transpose() << ")" << std::endl;
	}
}
