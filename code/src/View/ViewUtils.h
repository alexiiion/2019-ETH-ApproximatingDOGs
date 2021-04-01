#pragma once

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include "PatchModel.h"

void show_curve(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& curve, const Eigen::RowVector3d& color, bool show_points);
void render_highlighted_curves(igl::opengl::glfw::Viewer& viewer, const std::vector<Eigen::MatrixXd>& curves, const int selected_index, const Eigen::RowVector3d& color, const bool show_points);
void render_highlighted_curves(igl::opengl::glfw::Viewer& viewer, const std::vector<std::vector<Eigen::MatrixXd>>& curves_list, const int selected_index, const Eigen::RowVector3d& color, const bool show_points);
void show_curve(igl::opengl::glfw::Viewer& viewer, const std::vector<int>& curve_indices, const Eigen::MatrixXd& lookup, const Eigen::RowVector3d& color, bool show_points);
void label_quad_faces(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F_quad);

void get_curve_edges(const Eigen::MatrixXd& curve, Eigen::MatrixXd& out_start, Eigen::MatrixXd& out_end);
void get_quad_wireframe_edges(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, Eigen::MatrixXd & out_start, Eigen::MatrixXd & out_end);
void get_quad_wireframe_edges(const Patch& patch, Eigen::MatrixXd& out_start, Eigen::MatrixXd& out_end);

Eigen::RowVector3d get_color(int index, int number_colors);
