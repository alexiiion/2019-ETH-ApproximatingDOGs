#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
//#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
//#include <imgui/imgui.h>

#include "MeshModel.h"
#include "View/Model/MeshViewSettings.h"

class MeshViewer
{
public:
	MeshViewer(igl::opengl::glfw::Viewer& viewer, igl::opengl::glfw::imgui::ImGuiMenu& menu, Mesh& mesh, OverlaySettings& overlay_settings)
		: viewer(viewer), menu(menu), mesh(mesh), overlay_settings(overlay_settings)
	{
		label_assignment = NULL;
		add_mesh(mesh);
	};
	MeshViewer(igl::opengl::glfw::Viewer& viewer, igl::opengl::glfw::imgui::ImGuiMenu& menu, Mesh& mesh, OverlaySettings& overlay_settings, std::vector<int>*& label_assignment)
		: viewer(viewer), menu(menu), mesh(mesh), overlay_settings(overlay_settings), label_assignment(label_assignment)
	{
		add_mesh(mesh);
	};

	~MeshViewer() {};

	void update_view();
	void update_menu();
	void update_debug_menu();

	Mesh& mesh;
	MeshViewSettings mesh_settings;
	OverlaySettings& overlay_settings;
	std::vector<int>* label_assignment;

private:
	igl::opengl::glfw::Viewer& viewer;
	igl::opengl::glfw::imgui::ImGuiMenu& menu;

	void add_mesh(Mesh& mesh);

	void update_mesh();
	void update_mesh_view();
	void update_mesh_fill();

	void set_standard_fill();
	void visualize_curvature(const Eigen::VectorXd& curvature, const Eigen::VectorXi& curvature_lookup);
	void visualize_principal_curvature();
	void set_assignment_fill();
	void show_filtered_curvature_points(const Eigen::VectorXd& feature, const Eigen::VectorXi& feature_lookup);

};
