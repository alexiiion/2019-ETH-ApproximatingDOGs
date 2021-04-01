#include "ViewUtils.h"

#include <string>
#include <igl/hsv_to_rgb.h>

#include "Utils/Colors.h"
#include "Utils.h"


void show_curve(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& curve, const Eigen::RowVector3d& color, bool show_points)
{
	Eigen::MatrixXd edges_start, edges_end;
	get_curve_edges(curve, edges_start, edges_end);

	viewer.data().add_edges(edges_start, edges_end, color);

	if (show_points)
		viewer.data().add_points(curve, color);
}

void render_highlighted_curves(igl::opengl::glfw::Viewer& viewer, const std::vector<Eigen::MatrixXd>& curves, const int selected_index, const Eigen::RowVector3d& color, const bool show_points)
{
	if (curves.size() < 1)
		return;

	const auto color_unselected = Colors::darker(Colors::darker(color));

	//render the selected curve
	if (selected_index >= 0 && selected_index < curves.size())
		show_curve(viewer, curves[selected_index], color, show_points);

	//render all other curves
	for (int i = 0; i < curves.size(); i++)
	{
		if (i != selected_index)
			show_curve(viewer, curves[i], color_unselected, show_points);
	}
}

void render_highlighted_curves(igl::opengl::glfw::Viewer& viewer, const std::vector<std::vector<Eigen::MatrixXd>>& curves_list, const int selected_index, const Eigen::RowVector3d& color, const bool show_points)
{
	if (curves_list.size() < 1)
		return;

	const auto color_unselected = Colors::darker(Colors::darker(Colors::darker(color)));

	//render the selected curve
	if (selected_index >= 0 && selected_index < curves_list.size())
	{
		auto& curves = curves_list[selected_index];
		for (auto& curve : curves)
			show_curve(viewer, curve, color, show_points);
	}

	//render all other curves
	for (int i = 0; i < curves_list.size(); i++)
	{
		if (i == selected_index)
			continue;

		auto& curves = curves_list[i];
		for (auto& curve : curves)
			show_curve(viewer, curve, color_unselected, show_points);
	}
}

void show_curve(igl::opengl::glfw::Viewer& viewer, const std::vector<int>& curve_indices, const Eigen::MatrixXd& lookup, const Eigen::RowVector3d& color, bool show_points)
{
	Eigen::MatrixXd curve;
	index_to_value(curve_indices, lookup, curve);

	show_curve(viewer, curve, color, show_points);
}

void label_quad_faces(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F_quad)
{
	for (int i = 0; i < F_quad.rows(); ++i)
	{
		Eigen::RowVector3d p = Eigen::RowVector3d::Zero();
		for (int j = 0; j < F_quad.cols(); ++j)
			p += V.row(F_quad(i, j));

		p /= (double)F_quad.cols();

		viewer.data().add_label(p, std::to_string(i));
	}
}

void get_curve_edges(const Eigen::MatrixXd& curve, Eigen::MatrixXd& out_start, Eigen::MatrixXd& out_end)
{
	if (!curve.rows())
		return;

	out_start.resize(curve.rows() - 1, curve.cols());
	out_end.resize(curve.rows() - 1, curve.cols());

	for (int j = 0; j < curve.rows() - 1; j++)
	{
		out_start.row(j) = curve.row(j);
		out_end.row(j) = curve.row(j + 1);
	}
}

void get_quad_wireframe_edges(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, Eigen::MatrixXd & out_start, Eigen::MatrixXd & out_end)
{
	int number_edges = E.rows();

	out_start.resize(number_edges, 3);
	out_end.resize(number_edges, 3);
	
	for (int i = 0; i < number_edges; i++)
	{
		int v1 = E(i, 0);
		int v2 = E(i, 1);
		out_start.row(i) = V.row(v1);
		out_end.row(i) = V.row(v2);
	}
}

void get_quad_wireframe_edges(const Patch& patch, Eigen::MatrixXd & out_start, Eigen::MatrixXd & out_end)
{
	get_quad_wireframe_edges(patch.wrapper.V, patch.wrapper.quad_topology.E, out_start, out_end);
}

Eigen::RowVector3d get_color(int index, int number_colors)
{
	if (number_colors < 1)
		return Eigen::RowVector3d(0, 0, 0);

	double hue = index * 360 / number_colors;
	Eigen::RowVector3d hsv(hue, 1, 1);
	Eigen::RowVector3d rgb;
	igl::hsv_to_rgb(hsv, rgb);

	return rgb;
}
