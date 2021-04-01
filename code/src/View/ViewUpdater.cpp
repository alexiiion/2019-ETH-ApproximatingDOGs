#include "ViewUpdater.h"

#include <igl/per_vertex_normals.h>
#include <igl/edges.h>
#include <igl/avg_edge_length.h>
#include <igl/exact_geodesic.h>
#include <igl/boundary_loop.h>
#include <igl/hsv_to_rgb.h>
#include <igl/rgb_to_hsv.h>
#include <igl/per_corner_normals.h>

#include "GeodesicWalker.h"
#include "PatchModel.h"
#include "MeshController.h"
#include "ViewUtils.h"
#include "Utils/Colors.h"
#include "Utils.h"
#include "Logger.h"

#include "CurveHelper.h"
#include "TempUtils.h"

using namespace std;

Eigen::RowVector3d render_color_features = Colors::GREEN;

void red_blue_colormap__(const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors);
void two_color_colormap__(const Eigen::RowVector3d& color_min, const Eigen::RowVector3d& color_max, const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors);


void ViewUpdater::update_target_view(MeshViewSettings& mesh_settings, const Mesh& mesh)
{
	if (mesh.V.rows() < 1)
		return;

	if (mesh_settings.view_index == -1)
	{
		write_log(1) << endl << "ERROR in <" << __FUNCTION__ << "> target model is loaded but view index is not set. Bug -- go fix!" << endl << endl;
		exit(EXIT_FAILURE);
	}

	data_model.viewer.selected_data_index = mesh_settings.view_index;

	Eigen::RowVectorXd current_dimensions;
	meshhelper::calculate_dimensions(mesh.V, current_dimensions);

	bool is_target_scaled = !mesh_settings.dimensions.cols() || (current_dimensions - mesh_settings.dimensions).squaredNorm() > 1e-6;
	if (is_target_scaled)
	{
		data_model.viewer.data().clear();
		data_model.viewer.data().set_mesh(mesh.V, mesh.F);
		mesh_settings.dimensions = current_dimensions;
	}

	if (!mesh_settings.has_changed && !is_target_scaled)
		return;


	update_mesh_view(mesh_settings);
	//update_target_fill();

	//mesh_settings.has_changed = false;
}

void ViewUpdater::update_developables_parts_view()
{
	/*
	if (data_model.developable_model.developable_parts.size() < 1)
		return;
	if (!view_model.developable_parts_view_settings.has_changed)
		return;


	int selected_view_index = view_model.developable_parts_view_settings.view_index;

	if (selected_view_index == -1) //apply setting to all
	{
		for (int i = 0; i < data_model.developable_model.developable_parts.size(); i++)
		{
			view_model.developable_parts_view_settings.view_index = data_model.developable_parts_view_indices[i];
			update_target_view(view_model.developable_parts_view_settings, data_model.developable_model.developable_parts[i]);
		}
		view_model.developable_parts_view_settings.view_index = -1;
	}
	else 
	{
		for (int i = 0; i < data_model.developable_model.developable_parts.size(); i++)
		{
			int view_index = data_model.developable_parts_view_indices[i];

			if (view_index != selected_view_index) //hide unselected meshes
			{
				MeshViewSettings temp_settings = view_model.developable_parts_view_settings;
				temp_settings.view_index = view_index;
				temp_settings.show_mesh = false;
				temp_settings.has_changed = true;

				update_target_view(temp_settings, data_model.developable_model.developable_parts[i]);
			}
			else //update selected mesh
			{
				view_model.developable_parts_view_settings.view_index = view_index;
				update_target_view(view_model.developable_parts_view_settings, data_model.developable_model.developable_parts[i]);
			}
		}
	}
	*/
}

void ViewUpdater::update_patches_view()
{
	//check if new patch was added and add view settings for it
	if (data_model.patches.size() > view_model.patch_view_settings.size())
	{
		for (int i = view_model.patch_view_settings.size(); i < data_model.patches.size(); i++)
		{
			PatchViewSettings& view_settings = PatchViewSettings();
			view_settings.standard_fill_color = view_model.patch_color;

			view_settings.view_index = data_model.wrapper_view_indices[i];
			view_model.patch_view_settings.push_back(view_settings);

			if (view_model.selected_patch_index > 0)
			{
				view_model.patch_view_settings[view_model.selected_patch_index].show_position_constraints = false;
				view_model.patch_view_settings[view_model.selected_patch_index].show_target_constraints = false;
			}

			//automatically advance selected patch to newly added patch
			view_model.selected_patch_index = data_model.patches.size() - 1;
		}
	}

	if (data_model.patches.size() < 1)
		return;

	//updating all wrappers using general settings
	if (view_model.use_general_patch_settings || view_model.do_update_all_patches)
	{
		view_model.general_patch_view_settings.show_position_constraints = false;
		view_model.general_patch_view_settings.show_target_constraints = false;

		for (int i = 0; i < data_model.patches.size(); i++)
		{
			Patch& patch = *data_model.patches[i];
			data_model.viewer.selected_data_index = view_model.patch_view_settings[i].view_index;
			
			update_single_patch_view(view_model.general_patch_view_settings, patch);
		}
		view_model.do_update_all_patches = false;
	}
	else if (view_model.show_only_current_patch)
	{
		for (int i = 0; i < data_model.patches.size(); i++)
		{
			bool is_visible = i == view_model.selected_patch_index;

			Patch& patch = *data_model.patches[i];
			PatchViewSettings& view_settings = view_model.patch_view_settings[i];
			view_settings.show_mesh = is_visible;
			view_settings.show_position_constraints = is_visible && view_settings.show_position_constraints;
			view_settings.show_target_constraints = is_visible && view_settings.show_target_constraints;

			data_model.viewer.selected_data_index = view_settings.view_index;
			update_single_patch_view(view_settings, patch);
		}
	}
	//updating selected wrapper
	else
	{
		Patch& patch = *data_model.patches[view_model.selected_patch_index];
		PatchViewSettings& view_settings = view_model.patch_view_settings[view_model.selected_patch_index];

		data_model.viewer.selected_data_index = view_settings.view_index;
		update_single_patch_view(view_settings, patch);
	}
}

void ViewUpdater::update_single_patch_view(PatchViewSettings& view_settings, Patch& patch)
{
	//data_model.viewer.data().point_size = 10.0;
	//data_model.viewer.data().line_width = 1.0;

	view_settings.has_changed = true;

	//update mesh
	data_model.viewer.data().clear();
	data_model.viewer.data().set_mesh(patch.wrapper.V, patch.wrapper.F);
	update_mesh_view(view_settings);
	update_mesh_fill(view_settings, patch.wrapper);

	data_model.viewer.data().show_lines = false;
	if (view_settings.show_mesh && view_settings.show_wireframe)
	{
		Eigen::MatrixXd edges_start, edges_end;
		get_quad_wireframe_edges(patch, edges_start, edges_end);
		data_model.viewer.data().add_edges(edges_start, edges_end, Colors::GRAY_DARK);
	}

	if (view_settings.show_mesh && view_settings.label_quad_faces)
		label_quad_faces(data_model.viewer, patch.wrapper.V, patch.wrapper.F_quad);


	//visualize constraints
	if (view_model.show_debug_items)
	{
		data_model.viewer.selected_data_index = view_model.debug_view_index;

		if (view_settings.show_position_constraints)
		{
			if (patch.constraints.paired_points_from.rows() && patch.constraints.paired_points_to.rows() && patch.constraints.paired_points_from.rows() == patch.constraints.paired_points_to.rows())
			{
				data_model.viewer.data().add_edges(patch.constraints.paired_points_from, patch.constraints.paired_points_to, Colors::RED);
				data_model.viewer.data().add_points(patch.constraints.paired_points_to, Colors::RED);
				data_model.viewer.data().add_points(patch.constraints.paired_points_from, Colors::BLUE);
			}
		}
		if (view_settings.show_target_constraints)
		{
			const auto color = Colors::CYAN;
			int num_points = patch.target_curve.rows();

			for (int i = 0; i < num_points; i++)
			{
				auto current_color = color - (color * i / num_points);
				data_model.viewer.data().add_points(patch.target_curve.row(i), current_color);
			}
		}
	}

	view_settings.has_changed = false;
}

void ViewUpdater::update_gauss_map_view()
{
	if (view_model.do_update_gauss_map)
	{
		data_model.viewer.selected_data_index = view_model.gauss_map_view_index;
		data_model.viewer.data().clear();

		if (view_model.is_gauss_map_visible && view_model.gauss_mapped_mesh != NULL)
		{
			//show gauss map
			data_model.viewer.data().set_mesh(view_model.gauss_map_sphere.V, view_model.gauss_map_sphere.F);

			Eigen::Vector3d diffuse; diffuse << 0.98, 0.98, 0.98;
			Eigen::Vector3d ambient; ambient << 0, 0, 0;//0.05*diffuse;
			Eigen::Vector3d specular; specular << 0, 0, 0;
			data_model.viewer.data().uniform_colors(ambient, diffuse, specular);
			data_model.viewer.data().show_lines = false;
			//viewer.core().shininess = 0;

			Eigen::MatrixXd N;
			igl::per_vertex_normals(view_model.gauss_mapped_mesh->V, view_model.gauss_mapped_mesh->F, N);

			auto centroid = get_centroid(view_model.gauss_map_sphere.V);
			auto radius = view_model.gauss_map_sphere_radius + 0.1; //offset for visibility
			Eigen::MatrixXd N_sphere = (N * radius).rowwise() + centroid;
			data_model.viewer.data().add_points(N_sphere, Colors::BLACK);

			Eigen::MatrixXi E;
			igl::edges(view_model.gauss_mapped_mesh->F, E);

			Eigen::MatrixXd edges_start, edges_end;
			get_quad_wireframe_edges(N_sphere, E, edges_start, edges_end);
			data_model.viewer.data().add_edges(edges_start, edges_end, Colors::BLACK);

			//hide other meshes
			view_model.target_view_settings.show_mesh = false;
			//view_model.target_developable_view_settings.show_mesh = false;
			view_model.use_general_patch_settings = true;
			view_model.general_patch_view_settings.show_mesh = false;

			view_model.target_view_settings.has_changed = true;
			//view_model.target_developable_view_settings.has_changed = true;
			view_model.general_patch_view_settings.has_changed = true;
		}
		else
		{
			//show other meshes
			view_model.target_view_settings.show_mesh = true;
			//view_model.target_developable_view_settings.show_mesh = true;
			view_model.use_general_patch_settings = false;
			view_model.general_patch_view_settings.show_mesh = true;

			view_model.target_view_settings.has_changed = true;
			//view_model.target_developable_view_settings.has_changed = true;
			view_model.general_patch_view_settings.has_changed = true;
		}

		view_model.do_update_gauss_map = false;
	}
}

void ViewUpdater::update_debug_view()
{
	if (view_model.debug_view_index < 0)
		return;

	data_model.viewer.selected_data_index = view_model.debug_view_index;

#if show_coordinate_indicator
	const auto coordinate_indicator = Eigen::MatrixXd::Identity(3, 3);
	data_model.viewer.data().add_edges(Eigen::MatrixXd::Zero(3, 3), coordinate_indicator*0.2, coordinate_indicator);
#endif

	if (!view_model.show_debug_items)
		return;

	update_coverage_view();
	//update_surface_properties_view();
	update_geodesics_view();
	//update_cluster_view();
	update_result_view();
	//update_test_view();
}

void ViewUpdater::update_surface_properties_view()
{
	////if (!view_model.do_update_surface_labels)
	////	return;


	//if (view_model.show_filtered_curvature_points || view_model.show_filtered_curvature_values)
	//{
	//	//switch (view_model.current_target_visualization)
	//	switch (view_model.overlay_visualization)
	//	{
	//		case TargetVisualization::GaussianCurvature:
	//			show_filtered_curvature_points(data_model.target.surface_features().K, data_model.target.surface_features().K_abs_lookup);
	//			break;
	//		case TargetVisualization::MeanCurvature:
	//			show_filtered_curvature_points(data_model.target.surface_features().H, data_model.target.surface_features().H_lookup);
	//			break;
	//		//case TargetVisualization::WeightedCurvature:
	//		//	show_filtered_curvature_points(data_model.target.surface_features().weighted_curvature, data_model.target.surface_features().weighted_curvature_lookup);
	//		//	break;
	//		default:
	//		case TargetVisualization::Standard:
	//			break;
	//	}
	//}

	//update_crease_view();
	//view_model.do_update_surface_labels = false;
}

void ViewUpdater::update_feature_view(InteractiveSegmentSelector& segment_selection)
{
	/*
	//render preview points
	if (view_model.current_interaction_mode == InteractionMode::Select && view_model.current_selection_mode == SelectionMode::DefineFeatures)
	{
		Eigen::MatrixXd preview_segement_points;
		segment_selection.get_preview_segement_points(preview_segement_points);
		show_curve(preview_segement_points, Colors::DARK_GRAY, true);
	}


	if (!view_model.do_render_POI)
		return;

	if(!view_model.show_geodesics_candidates)
		render_color_features = Colors::GREEN;
	else 
		render_color_features = Colors::LIGHT_GRAY;

	vector<Eigen::MatrixXd> all_features;
	int size = data_model.clustered_feature_vertices.size();
	for (int i = 0; i < data_model.clustered_feature_vertices.size(); i++)
		all_features.push_back(segment_selection.get_segement_points_at(i));

	render_highlighted_curves(all_features, segment_selection.segment_index, render_color_features, true);
	*/
}

void ViewUpdater::update_geodesics_init_view()
{
	if (!view_model.do_update_init_geodesics || data_model.geodesics_candidates.is_empty())
		return;

	if (view_model.init_geodesics_view_index == -1)
	{
		data_model.viewer.append_mesh();
		view_model.init_geodesics_view_index = data_model.viewer.selected_data_index;
		data_model.viewer.data().point_size = 5.0;
		data_model.viewer.data().line_width = 1.25;
	}
	else
	{
		data_model.viewer.selected_data_index = view_model.init_geodesics_view_index;
		data_model.viewer.data().clear();
	}

	vector<int> indices_to_show = {};

	if (view_model.show_selected_geodesics && view_model.selected_ruled_index >= 0 && view_model.selected_ruled_index < data_model.selected_ruled_vertex_indices.size())
	{
		int vertex_index = data_model.selected_ruled_vertex_indices[view_model.selected_ruled_index];
		indices_to_show.push_back(vertex_index);
	}
	else
	{
		if (view_model.show_all_geodesics)
		{
			indices_to_show = vector<int>(data_model.target.V.rows());
			std::iota(indices_to_show.begin(), indices_to_show.end(), 0);
		}
		else if (view_model.show_random_geodesics)
		{
			indices_to_show = data_model.ruled_vertex_indices;
		}
		else if (view_model.show_selected_geodesics)
		{
			indices_to_show = data_model.selected_ruled_vertex_indices;
		}
	}

	//update points and curves
	for (int i = 0; i < indices_to_show.size(); i++)
	{
		int vertex = indices_to_show[i];

		if (view_model.show_curves_geodesics)
		{
			auto curve = data_model.geodesics_candidates.paths[vertex];
			show_curve(data_model.viewer, curve, Colors::BLACK, false);
		}

		if (view_model.show_points_geodesics)
		{
			data_model.viewer.data().add_points(data_model.target.V.row(vertex), Colors::BLACK);
		}
	}

	//hide all surfaces
	bool is_visible = false;

	if (data_model.selected_ruled_view_indices.size() < 1 && data_model.selected_ruled_vertex_indices.size() > 0 && data_model.ruled_developables.size() > 0) //mesh is not added to view yet
	{
		for (int geodesic_vertex_index : data_model.selected_ruled_vertex_indices)
		{
			int patch_index = index_of(geodesic_vertex_index, data_model.ruled_vertex_indices);
			if (patch_index == -1) //not contained
			{
				write_log(0) << "patch not found" << endl;
				continue;
			}

			RuledDevelopableSurface* surface = data_model.ruled_developables[patch_index];
			meshhelper::add_mesh(data_model.viewer, surface->developable.V, surface->developable.F, false, false);
			data_model.selected_ruled_view_indices.push_back(data_model.viewer.selected_data_index);
		}
	}

	for (int i = 0; i < data_model.selected_ruled_view_indices.size(); i++)
	{
		data_model.viewer.selected_data_index = data_model.selected_ruled_view_indices[i];
		data_model.viewer.data().show_faces = is_visible;
		data_model.viewer.data().show_lines = is_visible;
	}


	if (data_model.ruled_view_index == -1 && data_model.ruled_developables.size() > 0) //mesh is not added to view yet
	{
		Mesh concatenated_developables = meshhelper::concatenate_ruled_developables(data_model.ruled_developables);
		meshhelper::add_mesh(data_model.viewer, concatenated_developables.V, concatenated_developables.F, false, false);
		data_model.ruled_view_index = data_model.viewer.selected_data_index;
	}

	data_model.viewer.selected_data_index = data_model.ruled_view_index;
	data_model.viewer.data().show_faces = is_visible;
	data_model.viewer.data().show_lines = is_visible;


	if (view_model.show_surfaces_geodesics)
	{
		if (view_model.show_random_geodesics)
		{
			is_visible = true;

			data_model.viewer.selected_data_index = data_model.ruled_view_index;
			data_model.viewer.data().show_faces = is_visible;
			data_model.viewer.data().show_lines = is_visible;
		}
		else if (view_model.show_selected_geodesics && data_model.selected_ruled_view_indices.size() > 0)
		{
			//show only selected ones
			is_visible = true;
			////const double hue_step = 360 / data_model.selected_ruled_view_indices.size();
			//const double value_step = 1.0 / data_model.selected_ruled_view_indices.size();

			for (int i = 0; i < data_model.selected_ruled_view_indices.size(); i++)
			{
				int view_index = data_model.selected_ruled_view_indices[i];
				data_model.viewer.selected_data_index = view_index;

				data_model.viewer.data().show_faces = true;
				data_model.viewer.data().show_lines = true;
				data_model.viewer.data().set_colors(Colors::GRAY_LIGHT);

				////Eigen::RowVector3d hsv(i*hue_step, 1, 0.5);
				//Eigen::RowVector3d hsv(0, 0, i*value_step*0.5 + 0.25);
				//Eigen::RowVector3d rgb;
				//igl::hsv_to_rgb(hsv, rgb);
				//data_model.viewer.data().set_colors(rgb);
			}
		}
	}

	if (view_model.show_surfaces_geodesics && view_model.selected_ruled_index >= 0 && view_model.selected_ruled_index < data_model.selected_ruled_vertex_indices.size())
	{
		//int geodesic_vertex_index = data_model.selected_ruled_vertex_indices[view_model.selected_ruled_index];
		//int index = index_of(geodesic_vertex_index, data_model.ruled_vertex_indices);
		//int selected_view_index = data_model.ruled_view_indices[index];

		//for (int i = 0; i < data_model.selected_ruled_vertex_indices.size(); i++)
		//{
		//	int index = index_of(data_model.selected_ruled_vertex_indices[i], data_model.ruled_vertex_indices);
		//	int view_index = data_model.ruled_view_indices[index];
		//	data_model.viewer.selected_data_index = view_index;

		//	if (view_index == selected_view_index)
		//		is_visible = true;
		//	else
		//		is_visible = false;

		//	data_model.viewer.selected_data_index = view_index;
		//	data_model.viewer.data().show_faces = is_visible;
		//	data_model.viewer.data().show_lines = is_visible;
		//}

		int selected_view_index = data_model.selected_ruled_view_indices[view_model.selected_ruled_index];

		for (int i = 0; i < data_model.selected_ruled_vertex_indices.size(); i++)
		{
			int view_index = data_model.selected_ruled_view_indices[i];
			data_model.viewer.selected_data_index = view_index;

			if (view_index == selected_view_index)
				is_visible = true;
			else
				is_visible = false;

			data_model.viewer.selected_data_index = view_index;
			data_model.viewer.data().show_faces = is_visible;
			data_model.viewer.data().show_lines = is_visible;
		}

	}

	view_model.do_update_init_geodesics = false;
}

void ViewUpdater::update_geodesics_view()
{
}

void try_update_debug_geodesics_view(ViewModel& view_model, DataModel& data_model)
{
	if (!view_model.do_update_geodesics)
		return;

	if (view_model.do_render_debug_walking_geodesics)
	{
		if (view_model.debug_geodesics_view_index == -1)
		{
			data_model.viewer.append_mesh();
			view_model.debug_geodesics_view_index = data_model.viewer.selected_data_index;
			data_model.viewer.data().point_size = 10.0;
			data_model.viewer.data().line_width = 1.0;
		}
		else
		{
			data_model.viewer.selected_data_index = view_model.debug_geodesics_view_index;
			data_model.viewer.data().clear();
		}

		for (Eigen::MatrixXd geodesic : view_model.debug_walking_geodesics)
			show_curve(data_model.viewer, geodesic, Colors::BLUE, false);
	}
	//else
	//{
	//	if (view_model.geodesics_view_index == -1)
	//	{
	//		data_model.viewer.append_mesh();
	//		view_model.geodesics_view_index = data_model.viewer.selected_data_index;
	//		data_model.viewer.data().point_size = 5.0;
	//		data_model.viewer.data().line_width = 1.0;
	//	}
	//	else
	//	{
	//		data_model.viewer.selected_data_index = view_model.geodesics_view_index;
	//		data_model.viewer.data().clear();
	//	}
	//}

	view_model.do_update_geodesics = false;
}

void try_update_debug_rulings_view(ViewModel& view_model, DataModel& data_model)
{
	//view_model.selected_debug_ruled_surface_index = 544;
	//view_model.do_update_debug_ruled_surfaces = true;
	//view_model.is_debugging_ruled_surfaces = true;
	//view_model.use_debug_rulings_type = 3;

	if (!view_model.do_update_debug_ruled_surfaces || !view_model.is_debugging_ruled_surfaces)
		return;

	if (view_model.ruled_debug_view_index == -1)
	{
		data_model.viewer.append_mesh();
		view_model.ruled_debug_view_index = data_model.viewer.selected_data_index;
		data_model.viewer.data().point_size = 10.0;
		data_model.viewer.data().line_width = 1.0;
	}

	data_model.viewer.selected_data_index = view_model.ruled_debug_view_index;
	data_model.viewer.data().clear();
	view_model.do_update_debug_ruled_surfaces = false;

	if (!view_model.is_debugging_ruled_surfaces || view_model.selected_debug_ruled_surface_index < 0)
		return;


	int vertex = view_model.selected_debug_ruled_surface_index;
	write_log(4) << linebreak << linebreak << " -- updating debug ruled surface at v" << vertex << linebreak << endl;

	
	//get geodesic
	GeodesicWalker walker(data_model.target, data_model.geodesics_stopping);
	Eigen::MatrixXd geodesic;
	double length = 0;

	temp::get_geodesic(vertex, data_model.target, walker, data_model.geodesics_directions, geodesic, length);

	if (geodesic.rows() < 6)
	{
		view_model.selected_debug_ruled_surface_index = -1;
		return;
	}


	//get ruled surface
	RulingsType rulings_type;
	switch (view_model.use_debug_rulings_type)
	{
		case(0):
		{
			bool is_flat = curve::is_flat(geodesic);
			rulings_type = is_flat ? RulingsType::Discrete : RulingsType::Analytic;
			break;
		}
		case(1):
		{
			rulings_type = RulingsType::Analytic;
			break;
		}
		case(2):
		{
			rulings_type = RulingsType::Discrete;
			break;
		}
		case(3):
		{
			rulings_type = RulingsType::Compound;
			break;
		}
	}

	view_model.debug_ruled_surface = new RuledDevelopableSurface(data_model.target, geodesic, view_model.selected_debug_ruled_surface_index, rulings_type);
	view_model.debug_ruled_surface->create(data_model.ruled_width);
	//view_model.debug_ruled_surface->print_mathematica_data();

	show_curve(data_model.viewer, geodesic, Colors::BLUE, true);
	show_curve(data_model.viewer, view_model.debug_ruled_surface->generator_curve(), Colors::GRAY_MID, false);
	data_model.viewer.data().add_points(data_model.target.V.row(view_model.selected_debug_ruled_surface_index), Colors::BLUE);
	data_model.viewer.data().set_mesh(view_model.debug_ruled_surface->developable.V, view_model.debug_ruled_surface->developable.F);

	write_log(4) << linebreak << "  Ruled selection:: at v" << vertex << ": rows = " << geodesic.rows() << ", length = " << length << endl;
}

void ViewUpdater::update_debug_geodesics_view()
{
	try_update_debug_geodesics_view(view_model, data_model);
	try_update_debug_rulings_view(view_model, data_model);


	/*
	if (!view_model.do_update_geodesics)
		return;

	if (view_model.do_render_debug_walking_geodesics)
	{
		if (view_model.debug_geodesics_view_index == -1)
		{
			data_model.viewer.append_mesh();
			view_model.debug_geodesics_view_index = data_model.viewer.selected_data_index;
			data_model.viewer.data().point_size = 10.0;
			data_model.viewer.data().line_width = 1.0;
		}
		else
		{
			data_model.viewer.selected_data_index = view_model.debug_geodesics_view_index;
			data_model.viewer.data().clear();
		}

		for (Eigen::MatrixXd geodesic : view_model.debug_walking_geodesics)
			show_curve(data_model.viewer, geodesic, Colors::BLUE, false);
	}
	//else
	//{
	//	if (view_model.geodesics_view_index == -1)
	//	{
	//		data_model.viewer.append_mesh();
	//		view_model.geodesics_view_index = data_model.viewer.selected_data_index;
	//		data_model.viewer.data().point_size = 5.0;
	//		data_model.viewer.data().line_width = 1.0;
	//	}
	//	else
	//	{
	//		data_model.viewer.selected_data_index = view_model.geodesics_view_index;
	//		data_model.viewer.data().clear();
	//	}
	//}

	view_model.do_update_geodesics = false;
	*/
}

void ViewUpdater::update_coverage_view()
{
	/*
	//if (!view_model.show_coverage)
	//	return;

	if(!data_model.patch_coverage || data_model.patch_coverage == NULL)
		return;

	//if (!data_model.target.covering_wrapper_points.rows())
	if (!data_model.patch_coverage->covering_points().rows())
		return;

	if(view_model.show_coverage)
		data_model.viewer.data().add_edges(data_model.target.V, data_model.patch_coverage->covering_points(), Colors::RED);

	//if (data_model.target.is_covered())
	if (data_model.patch_coverage->is_covered())
		return;

	if (view_model.show_coverage_points)
	{
		Eigen::MatrixXd uncovered_vertices;
		index_to_value(data_model.patch_coverage->uncovered_vertices(), data_model.target.V, uncovered_vertices);
		data_model.viewer.data().add_points(uncovered_vertices, Eigen::RowVector3d(0.6, 0.0, 0.0)); //dark red
		//index_to_value(data_model.target.uncovered_V, data_model.target.V, uncovered_vertices);
	}
	*/

	if (view_model.current_coverage_choice == -1)
		return;

	if (!view_model.show_coverage && !view_model.show_coverage_points)
		return;

	Coverage* selected_coverage = data_model.patch_coverage;

	if (view_model.current_coverage_choice == 0)
		selected_coverage = data_model.initialization_coverage;
	else if (view_model.current_coverage_choice == 1)
		selected_coverage = data_model.patch_coverage;
	else if (view_model.current_coverage_choice == 2)
		selected_coverage = data_model.result_coverage;

	if (!selected_coverage || selected_coverage == NULL)
		return;

	if (!selected_coverage->covering_points().rows())
		return;

	if (view_model.show_coverage)
		data_model.viewer.data().add_edges(data_model.target.V, selected_coverage->covering_points(), Colors::RED);

	if (view_model.show_coverage_points && !selected_coverage->is_covered())
	{
		Eigen::MatrixXd uncovered_vertices;
		index_to_value(selected_coverage->uncovered_vertices(), data_model.target.V, uncovered_vertices);
		data_model.viewer.data().add_points(uncovered_vertices, Eigen::RowVector3d(0.6, 0.0, 0.0)); //dark red
	}




	//HOLES

	if (!view_model.show_coverage || !view_model.show_coverage_holes)
		return;

	if (data_model.hole_meshes.size() < 1)
		return;

	const int hue_step = 360 / data_model.hole_meshes.size();

	for(int i = 0; i < data_model.hole_meshes.size(); i++)
	{
		Eigen::MatrixXd uncovered_vertices = data_model.hole_meshes[i].V;

		Eigen::RowVector3d hsv(i*hue_step, 1, 1);
		Eigen::RowVector3d rgb;
		igl::hsv_to_rgb(hsv, rgb);
		
		data_model.viewer.data().add_points(uncovered_vertices, rgb);
		data_model.viewer.data().add_label(data_model.target.V.row(data_model.hole_average_vertex[i]), "a");

		//Eigen::RowVector3d average = uncovered_vertices.colwise().sum() / uncovered_vertices.rows();
		//data_model.viewer.data().add_points(average, Colors::BLACK);
	}

	Mesh uncovered = meshhelper::concatenate_meshes(data_model.hole_meshes);
	data_model.viewer.data().set_mesh(uncovered.V, uncovered.F);
	data_model.viewer.data().compute_normals();
	data_model.viewer.data().set_colors(Colors::GRAY_LIGHT);
}

void ViewUpdater::update_result_view()
{
	/*
	//visualize labeled faces
	if (data_model.vertex_patch_assignment.size() < 1)
		return;

	Eigen::MatrixXd F_midpoints(data_model.target_developable.F.rows(), 3);
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
	}

	int number_labels = *max_element(data_model.vertex_patch_assignment.begin(), data_model.vertex_patch_assignment.end()) + 1;
	int hue_step = 360 / number_labels;

	for (int li = 0; li < number_labels; li++)
	{
		Eigen::MatrixXd points(data_model.target_developable.F.rows(), 3);
		int n = 0;

		for (int i = 0; i < data_model.target_developable.F.rows(); i++)
		{
			if (data_model.vertex_patch_assignment[i] != li)
				continue;

			points.row(n) = F_midpoints.row(i);
			n++;
		}

		Eigen::RowVector3d hsv(li*hue_step, 1, 1);
		Eigen::RowVector3d rgb;
		igl::hsv_to_rgb(hsv, rgb);

		points.conservativeResize(n, Eigen::NoChange);
		data_model.viewer.data().add_points(points, rgb);
	}
	*/

	/*
	// visualize seams
	const int number_patches = data_model.patches.size();

	if (number_patches < 2)
		return;
	if (data_model.crease_indices.size() < 1)
		return;
	
	if (data_model.number_crease_pairs < 1)
		return;
	if (data_model.patch_neighbor_crease_lookup.rows() != number_patches || data_model.patch_neighbor_crease_lookup.cols() != number_patches)
		return;

	int hue_step = 360 / data_model.number_crease_pairs;
	for (int i = 0; i < data_model.number_crease_pairs; i++)
	{
		Eigen::MatrixXd vertices1;
		index_to_value(data_model.crease_indices[i], data_model.target_developable.V, vertices1);

		double hue = i * hue_step;
		hue = fmod((hue + 30.0), 360.0);
		Eigen::RowVector3d hsv(hue, 1, 1);
		Eigen::RowVector3d rgb;
		igl::hsv_to_rgb(hsv, rgb);

		data_model.viewer.data().add_points(vertices1, rgb);
		data_model.viewer.data().add_points(vertices1, Colors::darker(Colors::darker(rgb)));
	}
	*/

	//hue_step = 360 / data_model.number_crease_pairs;
	//for (int i = 0; i < data_model.number_crease_pairs; i++)
	//{
	//	Eigen::MatrixXd vertices1;
	//	index_to_value(data_model.crease_indices[i*2], data_model.target_developable.V, vertices1);
	//
	//	Eigen::MatrixXd vertices2;
	//	index_to_value(data_model.crease_indices[i*2+1], data_model.target_developable.V, vertices2);
	//
	//	double hue = i * hue_step;
	//	hue = fmod((hue + 30.0), 360.0);
	//	Eigen::RowVector3d hsv(hue, 1, 1);
	//	Eigen::RowVector3d rgb;
	//	igl::hsv_to_rgb(hsv, rgb);
	//
	//	data_model.viewer.data().add_points(vertices1, rgb);
	//	data_model.viewer.data().add_points(vertices1, Colors::darker(Colors::darker(rgb)));
	//}
	
	//for (int i = 0; i < data_model.number_crease_pairs; i++)
	//{
	//	double hue = i * hue_step;
	//	Eigen::RowVector3d hsv(hue, 1, 1);
	//	Eigen::RowVector3d rgb;
	//	igl::hsv_to_rgb(hsv, rgb);

	//	data_model.viewer.append_mesh();
	//	data_model.viewer.data().set_mesh(data_model.developable_parts[i].V, data_model.developable_parts[i].F);
	//	data_model.viewer.data().set_colors(rgb);
	//}



	//if (data_model.crease_indices[0].size() < 1)
	//	return;

	//Eigen::MatrixXd F_midpoints(data_model.target_developable.F.rows(), 3);
	//for (int fi = 0; fi < data_model.target_developable.F.rows(); fi++)
	//{
	//	Eigen::RowVectorXi face = data_model.target_developable.F.row(fi);
	//	Eigen::RowVectorXd face_sum = Eigen::RowVectorXd::Zero(3);
	//	for (int vi = 0; vi < 3; vi++)
	//	{
	//		if (face(vi) >= 0) //if it's -1 then there is no face because it's the boundary
	//			face_sum = face_sum + data_model.target_developable.V.row(face(vi));
	//	}

	//	F_midpoints.row(fi) = (face_sum / 3).eval();
	//}

	//Eigen::MatrixXd seam_midpoints;
	//index_to_value(data_model.crease_indices[0], F_midpoints, seam_midpoints);
	//data_model.viewer.data().add_points(F_midpoints, Colors::CYAN);





	/*
	Eigen::MatrixXd F_midpoints(data_model.target_developable.F.rows(), 3);
	Eigen::VectorXd F_areas(data_model.target_developable.F.rows());

	for (int fi = 0; fi < data_model.target_developable.F.rows(); fi++)
	{
		Eigen::RowVectorXi face = data_model.target_developable.F.row(fi);

		Eigen::RowVectorXd face_sum = Eigen::RowVectorXd::Zero(3);
		Eigen::RowVectorXd face_sides = Eigen::RowVectorXd::Zero(3);

		for (int vi = 0; vi < 3; vi++)
		{
			if (face(vi) < 0) //if it's -1 then there is no face because it's the boundary
				continue;

			face_sum = face_sum + data_model.target_developable.V.row(face(vi));
			face_sides(vi) = (data_model.target_developable.V.row(face(vi)) - data_model.target_developable.V.row(face((vi + 1) % 3))).norm();
		}

		F_midpoints.row(fi) = (face_sum / 3).eval();

		//Eigen::RowVectorXd a = data_model.target_developable.V.row(face(1)) - data_model.target_developable.V.row(face(0));
		//Eigen::RowVectorXd b = data_model.target_developable.V.row(face(2)) - data_model.target_developable.V.row(face(0));
		//double proj_ab = a.norm() * a.dot(b);

		double s = face_sides.sum() / 2;
		F_areas(fi) = sqrt(s * (s - face_sides(0)) * (s - face_sides(1)) * (s - face_sides(2)));

		data_model.viewer.data().add_label(F_midpoints.row(fi), to_string(F_areas(fi)));
	}

	data_model.viewer.data().add_points(F_midpoints, Colors::MAGENTA);
	*/
}

void ViewUpdater::update_mesh_view(const MeshViewSettings& mesh_settings)
{
	data_model.viewer.data().show_faces = mesh_settings.show_faces && mesh_settings.show_mesh;
	data_model.viewer.data().show_lines = mesh_settings.show_wireframe && mesh_settings.show_mesh;

	data_model.viewer.data().show_vertid = mesh_settings.label_vertices && mesh_settings.show_mesh;
	data_model.viewer.data().show_faceid = mesh_settings.label_faces && mesh_settings.show_mesh;
}

void ViewUpdater::update_target_fill()
{
	update_mesh_fill(view_model.target_view_settings, data_model.target);
	

	//add boundary first
	if (view_model.target_view_settings.show_boundary && view_model.target_view_settings.show_mesh)
	{
		vector<vector<int>> boundaries;
		igl::boundary_loop(data_model.target.F, boundaries);

		for (vector<int>& boundary : boundaries)
			show_curve(data_model.viewer, boundary, data_model.target.V, Colors::GRAY_DARK, false);
	}


	//show error on target, if selected
	if (view_model.target_view_settings.view_index == -1)
		return;
	//if (!view_model.target_view_settings.has_changed)
	//	return;
	if (data_model.target.V.rows() < 1)
		return;

	if (view_model.overlay_visualization == MeshVisualization::ErrorDistances)
	{
		if (view_model.current_coverage_choice == -1)
		{
			view_model.overlay_visualization = MeshVisualization::Standard;
			return;
		}

		Coverage* error;
		if (view_model.current_coverage_choice == 0)
			error = data_model.initialization_coverage;
		else if (view_model.current_coverage_choice == 1)
			error = data_model.patch_coverage;
		else if (view_model.current_coverage_choice == 2)
			error = data_model.result_coverage;
		else
			return;

		if (error == NULL)
			return;

		Eigen::VectorXd vertex_error = error->distances();
		if (vertex_error.rows() < 1)
			return;

		//fixed color range
		//double value_min = vertex_error.minCoeff();
		double value_max = vertex_error.maxCoeff();
		value_max = value_max > data_model.target_diagonal * 0.05 ? data_model.target_diagonal * 0.05 : value_max;
		write_log(0) << linebreak << "showing error to max = " << value_max << endl;

		//Eigen::MatrixXd colors;
		////const Eigen::RowVector3d purple = Eigen::RowVector3d(0.999, 0.001, 0.001);
		//two_color_colormap__(Colors::WHITE, Colors::RED, value_min, value_max, value_max, vertex_error, colors);
		//write_log(0) << "error colors.rows = " << colors.rows() << ",  colors.cols = " << colors.cols() << endl;

		//data_model.viewer.data().set_colors(colors);


		Eigen::MatrixXd colors(vertex_error.rows(), 3);
		Eigen::RowVector3d colorHSV;
		igl::rgb_to_hsv(Colors::RED, colorHSV);

		for (int i = 0; i < vertex_error.rows(); i++)
		{
			double saturation = (abs(vertex_error(i))) / value_max; //normalize to 1
			saturation = clip(saturation, 0.0, 1.0);
			//saturation = 1 - saturation;

			Eigen::RowVector3d hsv = colorHSV;
			hsv(1) = saturation;
			Eigen::RowVector3d rgb;
			igl::hsv_to_rgb(hsv, rgb);

			colors.row(i) = rgb;
			write_log(5) << "[" << i << "]: " << vertex_error(i) << ", normalized: " << (abs(vertex_error(i))) / value_max << ", added color: " << colors.row(i) << endl;
		}

		data_model.viewer.data().set_colors(colors);
		data_model.viewer.data().F_material_specular = Colors::BLACK;
		data_model.viewer.data().V_material_specular = Colors::BLACK;
	}
}

void ViewUpdater::update_mesh_fill(MeshViewSettings& mesh_settings, Mesh& mesh)
{
	if (mesh_settings.view_index == -1)
		return; 
	if (!mesh_settings.has_changed)
		return;
	if (mesh.V.rows() < 1)
		return;

	//data_model.viewer.selected_data_index = mesh_settings.view_index;
	
	//Eigen::MatrixXd empty_d(1, 3);
	//Eigen::MatrixXi empty_i(1, 3);
	//data_model.viewer.data().set_edges(empty_d, empty_i, empty_d); //clear edges from principal curvature
	data_model.viewer.data().lines.resize(0, Eigen::NoChange);

	//switch (view_model.current_target_visualization)
	switch (view_model.overlay_visualization)
	{
		case MeshVisualization::GaussianCurvature:
		{
			visualize_curvature(mesh.surface_features().K, mesh.surface_features().K_abs_lookup);
			break;
		}
		case MeshVisualization::GaussianCurvaturePrincipal:
		{
			visualize_curvature(mesh.surface_features().Kp, mesh.surface_features().Kp_abs_lookup);
			break;
		}
		case MeshVisualization::MeanCurvature:
		{
			visualize_curvature(mesh.surface_features().H, mesh.surface_features().H_lookup);
			break;
		}
		case MeshVisualization::PrincipalCurvature:
		{
			const double avg = igl::avg_edge_length(mesh.V, mesh.F);// *0.5;

			Eigen::MatrixXd N;
			igl::per_vertex_normals(mesh.V, mesh.F, N);

			const double offset = avg*0.2;
			N *= offset;

			//data_model.viewer.data().add_edges(data_model.target.V + N + data_model.target.surface_features().principal_k_min_direction * avg, data_model.target.V + N - data_model.target.surface_features().principal_k_min_direction * avg, Colors::BLUE);
			//data_model.viewer.data().add_edges(data_model.target.V + N + data_model.target.surface_features().principal_k_max_direction * avg, data_model.target.V + N - data_model.target.surface_features().principal_k_max_direction * avg, Colors::RED);
			
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
				edges_min_end.row(i)   = mesh.V.row(i) + N.row(i) - mesh.surface_features().principal_k_min_direction.row(i) * min_scale;
				edges_max_start.row(i) = mesh.V.row(i) + N.row(i) + mesh.surface_features().principal_k_max_direction.row(i) * max_scale;
				edges_max_end.row(i)   = mesh.V.row(i) + N.row(i) - mesh.surface_features().principal_k_max_direction.row(i) * max_scale;
			}

			data_model.viewer.data().add_edges(edges_min_start, edges_min_end, Colors::BLUE);
			data_model.viewer.data().add_edges(edges_max_start, edges_max_end, Colors::RED);

			break;
		}
		case MeshVisualization::LabelAssignment:
		{
			update_mesh_color_assignment_fill(mesh, data_model.vertex_ruled_assignment, data_model.selected_ruled_vertex_indices.size());
			break;
		}
		case MeshVisualization::VertexLabelling:
		{
			int number_labels = data_model.vertex_labels.maxCoeff();
			if (number_labels < 1)
				break;

			const double hue_step = 360 / number_labels;

			Eigen::MatrixXd colors;
			colors.resize(mesh.V.rows(), 3);


			for (int i = 0; i < data_model.vertex_labels.rows(); i++)
			{
				int label = data_model.vertex_labels(i);

				double hue = label * hue_step;
				double saturation = label == 0 ? 0 : 1;
				Eigen::RowVector3d hsv(hue, saturation, 1);

				Eigen::RowVector3d rgb;
				igl::hsv_to_rgb(hsv, rgb);
				colors.row(i) = Colors::ligher(Colors::ligher(rgb));
			}

			data_model.viewer.data().set_colors(colors);
			
			
			//Eigen::MatrixXd points;
			//index_to_value(data_model.vertex_ruled_assignment, data_model.target.V, points);
			//data_model.viewer.data().add_points(points, Colors::RED);
			
			break;
		}
		case MeshVisualization::FaceLabelling:
		{
			break;
		}
		case MeshVisualization::Standard:
		default:
		{
			//set_standard_fill();
			data_model.viewer.data().uniform_colors(Colors::GRAY_DARK, mesh_settings.standard_fill_color, Colors::BLACK);
			data_model.viewer.data().grid_texture();

			//Eigen::MatrixXd colors(mesh.V.rows(), mesh_settings.standard_fill_color.size());
			//for (int i = 0; i < mesh.V.rows(); i++)
			//	colors.row(i) = mesh_settings.standard_fill_color;
			//data_model.viewer.data().set_colors(colors);
			break;
		}
		//default:
		//{
		//	set_standard_fill();
		//	break;
		//}
	}

	mesh_settings.has_changed = false;
}

void ViewUpdater::update_mesh_color_assignment_fill(const Mesh& mesh, const std::vector<int>& label_assignment, const int number_labels)
{
	if (number_labels < 1)
		return;

	Eigen::MatrixXd colors;

	if (label_assignment.size() < 0)
	{
		colors.resize(1, 3);
		colors.row(0) = Eigen::RowVector3d(1.0, 1.0, 1.0);
	}
	else
	{
		const double hue_step = 360 / number_labels;
		colors.resize(mesh.V.rows(), 3);

		for (int i = 0; i < mesh.V.rows(); i++)
		{
			int patch_index = label_assignment[i];

			double hue = patch_index * hue_step;
			Eigen::RowVector3d hsv(hue, 1, 1);
			Eigen::RowVector3d rgb;
			igl::hsv_to_rgb(hsv, rgb);
			colors.row(i) = Colors::ligher(Colors::ligher(rgb));
		}
	}

	data_model.viewer.data().set_colors(colors);
}


void ViewUpdater::set_standard_fill()
{
	data_model.viewer.data().uniform_colors(
		Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
		Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
		Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

	//data_model.viewer.data_list[data_model.target_view_index].grid_texture();
	data_model.viewer.data().grid_texture();
}

void ViewUpdater::visualize_curvature(const Eigen::VectorXd& curvature, const Eigen::VectorXi& curvature_lookup)
{
	if (!curvature.rows())
		return;

	//int end_index = view_model.filter_curvature * (curvature.rows() - 1);
	//view_model.visualized_max = abs(curvature(curvature_lookup(end_index)));
	//view_model.visualized_min = view_model.visualized_max * -1;

	//int visualized_range_index = 0.3 * (curvature.rows() - 1); //clip top 30% of values to get rid of extremes and show a more consistent color range
	//double visualized_range = abs(curvature(curvature_lookup(visualized_range_index)));

	double curve_max = curvature.cwiseAbs().maxCoeff();
	double visualized_range = view_model.visualized_max;

	if (view_model.visualized_max > curve_max)
	{
		view_model.visualized_max = curve_max;
		int visualized_range_index = 0.3 * (curvature.rows() - 1); //clip top 30% of values to get rid of extremes and show a more consistent color range
		visualized_range = abs(curvature(curvature_lookup(visualized_range_index)));
	}


	Eigen::MatrixXd C;
	//igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_JET, curvature, view_model.visualized_min, view_model.visualized_max, C);
	red_blue_colormap__(view_model.visualized_min, view_model.visualized_max, visualized_range, curvature, C);
	data_model.viewer.data().set_colors(C);
	//data_model.viewer.data_list[view_model.target_view_settings.view_index].set_colors(C);
}

void ViewUpdater::show_filtered_curvature_points(const Eigen::VectorXd& feature, const Eigen::VectorXi& feature_lookup)
{
	if (!feature.rows())
		return;

	int end_index = view_model.filter_curvature * (feature.rows() - 1);

	for (int i = 0; i <= end_index; i++)
	{
		int vertex_index = feature_lookup(i);

		if (view_model.show_filtered_curvature_points)
			data_model.viewer.data().add_points(data_model.target.V.row(vertex_index), Eigen::RowVector3d(0.5, 0.5, 0.5));
		if (view_model.show_filtered_curvature_values)
			data_model.viewer.data().add_label(data_model.target.V.row(vertex_index), std::to_string(feature(vertex_index)));
	}
}

void red_blue_colormap__(const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors)
{
	two_color_colormap__(Colors::BLUE, Colors::RED, value_min, value_max, value_range, values, out_colors);
}

void two_color_colormap__(const Eigen::RowVector3d& color_min, const Eigen::RowVector3d& color_max, const double& value_min, const double& value_max, const double& value_range, const Eigen::VectorXd& values, Eigen::MatrixXd& out_colors)
{
	const int loglevel = 6;

	out_colors.resize(values.rows(), 3);
	Eigen::RowVector3d color_min_mask = (color_min - Eigen::RowVector3d::Ones()).cwiseAbs();
	Eigen::RowVector3d color_max_mask = (color_max - Eigen::RowVector3d::Ones()).cwiseAbs();

	write_log(loglevel) << endl;
	write_log(loglevel) << "values: " << values << endl;
	write_log(loglevel) << "color_min: " << color_min.transpose() << ", mask: " << color_min_mask.transpose() << endl;
	write_log(loglevel) << "color_max: " << color_max.transpose() << ", mask: " << color_max_mask.transpose() << endl;
	write_log(loglevel) << "value_min: " << value_min << ", value_max: " << value_max << ", denom: " << value_range << ", values.maxCoeff(): " << values.maxCoeff() << ", values.minCoeff(): " << values.minCoeff() << endl << endl;

	for (int i = 0; i < values.rows(); i++)
	{
		Eigen::Vector3d color = values(i) > 0 ? color_max : color_min;
		Eigen::Vector3d color_mask = values(i) > 0 ? color_max_mask : color_min_mask;
		double threshold = values(i) > 0 ? value_max : abs(value_min);

		double saturation = (abs(values(i)) - threshold) / value_range; //normalize to 1
		saturation = clip(saturation, 0.0, 1.0);
		saturation = 1 - saturation;

		out_colors.row(i) = color + color_mask * saturation;
		
		write_log(loglevel) << "[" << i << "]: " << values(i) << ", normalized: " << (values(i) - threshold) / value_range << ", added color: " << out_colors.row(i) << " (selected color: " << color.transpose() << ", mask: " << color_mask.transpose() << ")" << endl;
	}
}


// colormap seismic: https://matplotlib.org/examples/color/colormaps_reference.html
/*
let seismic = [
	[0.000, [0.000, 0.000, 0.300]],
		[0.002, [0.000, 0.000, 0.300]],
		[0.004, [0.000, 0.000, 0.311]],
		[0.006, [0.000, 0.000, 0.311]],
		[0.008, [0.000, 0.000, 0.322]],
		[0.010, [0.000, 0.000, 0.322]],
		[0.012, [0.000, 0.000, 0.333]],
		[0.014, [0.000, 0.000, 0.333]],
		[0.016, [0.000, 0.000, 0.344]],
		[0.018, [0.000, 0.000, 0.344]],
		[0.020, [0.000, 0.000, 0.355]],
		[0.022, [0.000, 0.000, 0.355]],
		[0.023, [0.000, 0.000, 0.366]],
		[0.025, [0.000, 0.000, 0.366]],
		[0.027, [0.000, 0.000, 0.377]],
		[0.029, [0.000, 0.000, 0.377]],
		[0.031, [0.000, 0.000, 0.388]],
		[0.033, [0.000, 0.000, 0.388]],
		[0.035, [0.000, 0.000, 0.399]],
		[0.037, [0.000, 0.000, 0.399]],
		[0.039, [0.000, 0.000, 0.410]],
		[0.041, [0.000, 0.000, 0.410]],
		[0.043, [0.000, 0.000, 0.421]],
		[0.045, [0.000, 0.000, 0.421]],
		[0.047, [0.000, 0.000, 0.432]],
		[0.049, [0.000, 0.000, 0.432]],
		[0.051, [0.000, 0.000, 0.443]],
		[0.053, [0.000, 0.000, 0.443]],
		[0.055, [0.000, 0.000, 0.454]],
		[0.057, [0.000, 0.000, 0.454]],
		[0.059, [0.000, 0.000, 0.465]],
		[0.061, [0.000, 0.000, 0.465]],
		[0.063, [0.000, 0.000, 0.476]],
		[0.065, [0.000, 0.000, 0.476]],
		[0.067, [0.000, 0.000, 0.487]],
		[0.068, [0.000, 0.000, 0.487]],
		[0.070, [0.000, 0.000, 0.498]],
		[0.072, [0.000, 0.000, 0.498]],
		[0.074, [0.000, 0.000, 0.509]],
		[0.076, [0.000, 0.000, 0.509]],
		[0.078, [0.000, 0.000, 0.520]],
		[0.080, [0.000, 0.000, 0.520]],
		[0.082, [0.000, 0.000, 0.531]],
		[0.084, [0.000, 0.000, 0.531]],
		[0.086, [0.000, 0.000, 0.542]],
		[0.088, [0.000, 0.000, 0.542]],
		[0.090, [0.000, 0.000, 0.553]],
		[0.092, [0.000, 0.000, 0.553]],
		[0.094, [0.000, 0.000, 0.564]],
		[0.096, [0.000, 0.000, 0.564]],
		[0.098, [0.000, 0.000, 0.575]],
		[0.100, [0.000, 0.000, 0.575]],
		[0.102, [0.000, 0.000, 0.585]],
		[0.104, [0.000, 0.000, 0.585]],
		[0.106, [0.000, 0.000, 0.596]],
		[0.108, [0.000, 0.000, 0.596]],
		[0.110, [0.000, 0.000, 0.607]],
		[0.112, [0.000, 0.000, 0.607]],
		[0.114, [0.000, 0.000, 0.618]],
		[0.115, [0.000, 0.000, 0.618]],
		[0.117, [0.000, 0.000, 0.629]],
		[0.119, [0.000, 0.000, 0.629]],
		[0.121, [0.000, 0.000, 0.640]],
		[0.123, [0.000, 0.000, 0.640]],
		[0.125, [0.000, 0.000, 0.651]],
		[0.127, [0.000, 0.000, 0.651]],
		[0.129, [0.000, 0.000, 0.662]],
		[0.131, [0.000, 0.000, 0.662]],
		[0.133, [0.000, 0.000, 0.673]],
		[0.135, [0.000, 0.000, 0.673]],
		[0.137, [0.000, 0.000, 0.684]],
		[0.139, [0.000, 0.000, 0.684]],
		[0.141, [0.000, 0.000, 0.695]],
		[0.143, [0.000, 0.000, 0.695]],
		[0.145, [0.000, 0.000, 0.706]],
		[0.147, [0.000, 0.000, 0.706]],
		[0.149, [0.000, 0.000, 0.717]],
		[0.151, [0.000, 0.000, 0.717]],
		[0.153, [0.000, 0.000, 0.728]],
		[0.155, [0.000, 0.000, 0.728]],
		[0.157, [0.000, 0.000, 0.739]],
		[0.159, [0.000, 0.000, 0.739]],
		[0.160, [0.000, 0.000, 0.750]],
		[0.162, [0.000, 0.000, 0.750]],
		[0.164, [0.000, 0.000, 0.761]],
		[0.166, [0.000, 0.000, 0.761]],
		[0.168, [0.000, 0.000, 0.772]],
		[0.170, [0.000, 0.000, 0.772]],
		[0.172, [0.000, 0.000, 0.783]],
		[0.174, [0.000, 0.000, 0.783]],
		[0.176, [0.000, 0.000, 0.794]],
		[0.178, [0.000, 0.000, 0.794]],
		[0.180, [0.000, 0.000, 0.805]],
		[0.182, [0.000, 0.000, 0.805]],
		[0.184, [0.000, 0.000, 0.816]],
		[0.186, [0.000, 0.000, 0.816]],
		[0.188, [0.000, 0.000, 0.827]],
		[0.190, [0.000, 0.000, 0.827]],
		[0.192, [0.000, 0.000, 0.838]],
		[0.194, [0.000, 0.000, 0.838]],
		[0.196, [0.000, 0.000, 0.849]],
		[0.198, [0.000, 0.000, 0.849]],
		[0.200, [0.000, 0.000, 0.860]],
		[0.202, [0.000, 0.000, 0.860]],
		[0.204, [0.000, 0.000, 0.871]],
		[0.205, [0.000, 0.000, 0.871]],
		[0.207, [0.000, 0.000, 0.882]],
		[0.209, [0.000, 0.000, 0.882]],
		[0.211, [0.000, 0.000, 0.893]],
		[0.213, [0.000, 0.000, 0.893]],
		[0.215, [0.000, 0.000, 0.904]],
		[0.217, [0.000, 0.000, 0.904]],
		[0.219, [0.000, 0.000, 0.915]],
		[0.221, [0.000, 0.000, 0.915]],
		[0.223, [0.000, 0.000, 0.926]],
		[0.225, [0.000, 0.000, 0.926]],
		[0.227, [0.000, 0.000, 0.937]],
		[0.229, [0.000, 0.000, 0.937]],
		[0.231, [0.000, 0.000, 0.948]],
		[0.233, [0.000, 0.000, 0.948]],
		[0.235, [0.000, 0.000, 0.959]],
		[0.237, [0.000, 0.000, 0.959]],
		[0.239, [0.000, 0.000, 0.970]],
		[0.241, [0.000, 0.000, 0.970]],
		[0.243, [0.000, 0.000, 0.981]],
		[0.245, [0.000, 0.000, 0.981]],
		[0.247, [0.000, 0.000, 0.992]],
		[0.249, [0.000, 0.000, 0.992]],
		[0.250, [0.004, 0.004, 1.000]],
		[0.252, [0.004, 0.004, 1.000]],
		[0.254, [0.020, 0.020, 1.000]],
		[0.256, [0.020, 0.020, 1.000]],
		[0.258, [0.035, 0.035, 1.000]],
		[0.260, [0.035, 0.035, 1.000]],
		[0.262, [0.051, 0.051, 1.000]],
		[0.264, [0.051, 0.051, 1.000]],
		[0.266, [0.067, 0.067, 1.000]],
		[0.268, [0.067, 0.067, 1.000]],
		[0.270, [0.082, 0.082, 1.000]],
		[0.272, [0.082, 0.082, 1.000]],
		[0.274, [0.098, 0.098, 1.000]],
		[0.276, [0.098, 0.098, 1.000]],
		[0.278, [0.114, 0.114, 1.000]],
		[0.280, [0.114, 0.114, 1.000]],
		[0.282, [0.129, 0.129, 1.000]],
		[0.284, [0.129, 0.129, 1.000]],
		[0.286, [0.145, 0.145, 1.000]],
		[0.288, [0.145, 0.145, 1.000]],
		[0.290, [0.161, 0.161, 1.000]],
		[0.292, [0.161, 0.161, 1.000]],
		[0.294, [0.176, 0.176, 1.000]],
		[0.295, [0.176, 0.176, 1.000]],
		[0.297, [0.192, 0.192, 1.000]],
		[0.299, [0.192, 0.192, 1.000]],
		[0.301, [0.208, 0.208, 1.000]],
		[0.303, [0.208, 0.208, 1.000]],
		[0.305, [0.224, 0.224, 1.000]],
		[0.307, [0.224, 0.224, 1.000]],
		[0.309, [0.239, 0.239, 1.000]],
		[0.311, [0.239, 0.239, 1.000]],
		[0.313, [0.255, 0.255, 1.000]],
		[0.315, [0.255, 0.255, 1.000]],
		[0.317, [0.271, 0.271, 1.000]],
		[0.319, [0.271, 0.271, 1.000]],
		[0.321, [0.286, 0.286, 1.000]],
		[0.323, [0.286, 0.286, 1.000]],
		[0.325, [0.302, 0.302, 1.000]],
		[0.327, [0.302, 0.302, 1.000]],
		[0.329, [0.318, 0.318, 1.000]],
		[0.331, [0.318, 0.318, 1.000]],
		[0.333, [0.333, 0.333, 1.000]],
		[0.335, [0.333, 0.333, 1.000]],
		[0.337, [0.349, 0.349, 1.000]],
		[0.339, [0.349, 0.349, 1.000]],
		[0.341, [0.365, 0.365, 1.000]],
		[0.342, [0.365, 0.365, 1.000]],
		[0.344, [0.380, 0.380, 1.000]],
		[0.346, [0.380, 0.380, 1.000]],
		[0.348, [0.396, 0.396, 1.000]],
		[0.350, [0.396, 0.396, 1.000]],
		[0.352, [0.412, 0.412, 1.000]],
		[0.354, [0.412, 0.412, 1.000]],
		[0.356, [0.427, 0.427, 1.000]],
		[0.358, [0.427, 0.427, 1.000]],
		[0.360, [0.443, 0.443, 1.000]],
		[0.362, [0.443, 0.443, 1.000]],
		[0.364, [0.459, 0.459, 1.000]],
		[0.366, [0.459, 0.459, 1.000]],
		[0.368, [0.475, 0.475, 1.000]],
		[0.370, [0.475, 0.475, 1.000]],
		[0.372, [0.490, 0.490, 1.000]],
		[0.374, [0.490, 0.490, 1.000]],
		[0.376, [0.506, 0.506, 1.000]],
		[0.378, [0.506, 0.506, 1.000]],
		[0.380, [0.522, 0.522, 1.000]],
		[0.382, [0.522, 0.522, 1.000]],
		[0.384, [0.537, 0.537, 1.000]],
		[0.386, [0.537, 0.537, 1.000]],
		[0.387, [0.553, 0.553, 1.000]],
		[0.389, [0.553, 0.553, 1.000]],
		[0.391, [0.569, 0.569, 1.000]],
		[0.393, [0.569, 0.569, 1.000]],
		[0.395, [0.584, 0.584, 1.000]],
		[0.397, [0.584, 0.584, 1.000]],
		[0.399, [0.600, 0.600, 1.000]],
		[0.401, [0.600, 0.600, 1.000]],
		[0.403, [0.616, 0.616, 1.000]],
		[0.405, [0.616, 0.616, 1.000]],
		[0.407, [0.631, 0.631, 1.000]],
		[0.409, [0.631, 0.631, 1.000]],
		[0.411, [0.647, 0.647, 1.000]],
		[0.413, [0.647, 0.647, 1.000]],
		[0.415, [0.663, 0.663, 1.000]],
		[0.417, [0.663, 0.663, 1.000]],
		[0.419, [0.678, 0.678, 1.000]],
		[0.421, [0.678, 0.678, 1.000]],
		[0.423, [0.694, 0.694, 1.000]],
		[0.425, [0.694, 0.694, 1.000]],
		[0.427, [0.710, 0.710, 1.000]],
		[0.429, [0.710, 0.710, 1.000]],
		[0.431, [0.725, 0.725, 1.000]],
		[0.432, [0.725, 0.725, 1.000]],
		[0.434, [0.741, 0.741, 1.000]],
		[0.436, [0.741, 0.741, 1.000]],
		[0.438, [0.757, 0.757, 1.000]],
		[0.440, [0.757, 0.757, 1.000]],
		[0.442, [0.773, 0.773, 1.000]],
		[0.444, [0.773, 0.773, 1.000]],
		[0.446, [0.788, 0.788, 1.000]],
		[0.448, [0.788, 0.788, 1.000]],
		[0.450, [0.804, 0.804, 1.000]],
		[0.452, [0.804, 0.804, 1.000]],
		[0.454, [0.820, 0.820, 1.000]],
		[0.456, [0.820, 0.820, 1.000]],
		[0.458, [0.835, 0.835, 1.000]],
		[0.460, [0.835, 0.835, 1.000]],
		[0.462, [0.851, 0.851, 1.000]],
		[0.464, [0.851, 0.851, 1.000]],
		[0.466, [0.867, 0.867, 1.000]],
		[0.468, [0.867, 0.867, 1.000]],
		[0.470, [0.882, 0.882, 1.000]],
		[0.472, [0.882, 0.882, 1.000]],
		[0.474, [0.898, 0.898, 1.000]],
		[0.476, [0.898, 0.898, 1.000]],
		[0.477, [0.914, 0.914, 1.000]],
		[0.479, [0.914, 0.914, 1.000]],
		[0.481, [0.929, 0.929, 1.000]],
		[0.483, [0.929, 0.929, 1.000]],
		[0.485, [0.945, 0.945, 1.000]],
		[0.487, [0.945, 0.945, 1.000]],
		[0.489, [0.961, 0.961, 1.000]],
		[0.491, [0.961, 0.961, 1.000]],
		[0.493, [0.976, 0.976, 1.000]],
		[0.495, [0.976, 0.976, 1.000]],
		[0.497, [0.992, 0.992, 1.000]],
		[0.499, [0.992, 0.992, 1.000]],
		[0.501, [1.000, 0.992, 0.992]],
		[0.503, [1.000, 0.992, 0.992]],
		[0.505, [1.000, 0.976, 0.976]],
		[0.507, [1.000, 0.976, 0.976]],
		[0.509, [1.000, 0.961, 0.961]],
		[0.511, [1.000, 0.961, 0.961]],
		[0.513, [1.000, 0.945, 0.945]],
		[0.515, [1.000, 0.945, 0.945]],
		[0.517, [1.000, 0.929, 0.929]],
		[0.519, [1.000, 0.929, 0.929]],
		[0.521, [1.000, 0.914, 0.914]],
		[0.523, [1.000, 0.914, 0.914]],
		[0.524, [1.000, 0.898, 0.898]],
		[0.526, [1.000, 0.898, 0.898]],
		[0.528, [1.000, 0.882, 0.882]],
		[0.530, [1.000, 0.882, 0.882]],
		[0.532, [1.000, 0.867, 0.867]],
		[0.534, [1.000, 0.867, 0.867]],
		[0.536, [1.000, 0.851, 0.851]],
		[0.538, [1.000, 0.851, 0.851]],
		[0.540, [1.000, 0.835, 0.835]],
		[0.542, [1.000, 0.835, 0.835]],
		[0.544, [1.000, 0.820, 0.820]],
		[0.546, [1.000, 0.820, 0.820]],
		[0.548, [1.000, 0.804, 0.804]],
		[0.550, [1.000, 0.804, 0.804]],
		[0.552, [1.000, 0.788, 0.788]],
		[0.554, [1.000, 0.788, 0.788]],
		[0.556, [1.000, 0.773, 0.773]],
		[0.558, [1.000, 0.773, 0.773]],
		[0.560, [1.000, 0.757, 0.757]],
		[0.562, [1.000, 0.757, 0.757]],
		[0.564, [1.000, 0.741, 0.741]],
		[0.566, [1.000, 0.741, 0.741]],
		[0.568, [1.000, 0.725, 0.725]],
		[0.569, [1.000, 0.725, 0.725]],
		[0.571, [1.000, 0.710, 0.710]],
		[0.573, [1.000, 0.710, 0.710]],
		[0.575, [1.000, 0.694, 0.694]],
		[0.577, [1.000, 0.694, 0.694]],
		[0.579, [1.000, 0.678, 0.678]],
		[0.581, [1.000, 0.678, 0.678]],
		[0.583, [1.000, 0.663, 0.663]],
		[0.585, [1.000, 0.663, 0.663]],
		[0.587, [1.000, 0.647, 0.647]],
		[0.589, [1.000, 0.647, 0.647]],
		[0.591, [1.000, 0.631, 0.631]],
		[0.593, [1.000, 0.631, 0.631]],
		[0.595, [1.000, 0.616, 0.616]],
		[0.597, [1.000, 0.616, 0.616]],
		[0.599, [1.000, 0.600, 0.600]],
		[0.601, [1.000, 0.600, 0.600]],
		[0.603, [1.000, 0.584, 0.584]],
		[0.605, [1.000, 0.584, 0.584]],
		[0.607, [1.000, 0.569, 0.569]],
		[0.609, [1.000, 0.569, 0.569]],
		[0.611, [1.000, 0.553, 0.553]],
		[0.613, [1.000, 0.553, 0.553]],
		[0.614, [1.000, 0.537, 0.537]],
		[0.616, [1.000, 0.537, 0.537]],
		[0.618, [1.000, 0.522, 0.522]],
		[0.620, [1.000, 0.522, 0.522]],
		[0.622, [1.000, 0.506, 0.506]],
		[0.624, [1.000, 0.506, 0.506]],
		[0.626, [1.000, 0.490, 0.490]],
		[0.628, [1.000, 0.490, 0.490]],
		[0.630, [1.000, 0.475, 0.475]],
		[0.632, [1.000, 0.475, 0.475]],
		[0.634, [1.000, 0.459, 0.459]],
		[0.636, [1.000, 0.459, 0.459]],
		[0.638, [1.000, 0.443, 0.443]],
		[0.640, [1.000, 0.443, 0.443]],
		[0.642, [1.000, 0.427, 0.427]],
		[0.644, [1.000, 0.427, 0.427]],
		[0.646, [1.000, 0.412, 0.412]],
		[0.648, [1.000, 0.412, 0.412]],
		[0.650, [1.000, 0.396, 0.396]],
		[0.652, [1.000, 0.396, 0.396]],
		[0.654, [1.000, 0.380, 0.380]],
		[0.656, [1.000, 0.380, 0.380]],
		[0.658, [1.000, 0.365, 0.365]],
		[0.659, [1.000, 0.365, 0.365]],
		[0.661, [1.000, 0.349, 0.349]],
		[0.663, [1.000, 0.349, 0.349]],
		[0.665, [1.000, 0.333, 0.333]],
		[0.667, [1.000, 0.333, 0.333]],
		[0.669, [1.000, 0.318, 0.318]],
		[0.671, [1.000, 0.318, 0.318]],
		[0.673, [1.000, 0.302, 0.302]],
		[0.675, [1.000, 0.302, 0.302]],
		[0.677, [1.000, 0.286, 0.286]],
		[0.679, [1.000, 0.286, 0.286]],
		[0.681, [1.000, 0.271, 0.271]],
		[0.683, [1.000, 0.271, 0.271]],
		[0.685, [1.000, 0.255, 0.255]],
		[0.687, [1.000, 0.255, 0.255]],
		[0.689, [1.000, 0.239, 0.239]],
		[0.691, [1.000, 0.239, 0.239]],
		[0.693, [1.000, 0.224, 0.224]],
		[0.695, [1.000, 0.224, 0.224]],
		[0.697, [1.000, 0.208, 0.208]],
		[0.699, [1.000, 0.208, 0.208]],
		[0.701, [1.000, 0.192, 0.192]],
		[0.703, [1.000, 0.192, 0.192]],
		[0.705, [1.000, 0.176, 0.176]],
		[0.706, [1.000, 0.176, 0.176]],
		[0.708, [1.000, 0.161, 0.161]],
		[0.710, [1.000, 0.161, 0.161]],
		[0.712, [1.000, 0.145, 0.145]],
		[0.714, [1.000, 0.145, 0.145]],
		[0.716, [1.000, 0.129, 0.129]],
		[0.718, [1.000, 0.129, 0.129]],
		[0.720, [1.000, 0.114, 0.114]],
		[0.722, [1.000, 0.114, 0.114]],
		[0.724, [1.000, 0.098, 0.098]],
		[0.726, [1.000, 0.098, 0.098]],
		[0.728, [1.000, 0.082, 0.082]],
		[0.730, [1.000, 0.082, 0.082]],
		[0.732, [1.000, 0.067, 0.067]],
		[0.734, [1.000, 0.067, 0.067]],
		[0.736, [1.000, 0.051, 0.051]],
		[0.738, [1.000, 0.051, 0.051]],
		[0.740, [1.000, 0.035, 0.035]],
		[0.742, [1.000, 0.035, 0.035]],
		[0.744, [1.000, 0.020, 0.020]],
		[0.746, [1.000, 0.020, 0.020]],
		[0.748, [1.000, 0.004, 0.004]],
		[0.750, [1.000, 0.004, 0.004]],
		[0.751, [0.994, 0.000, 0.000]],
		[0.753, [0.994, 0.000, 0.000]],
		[0.755, [0.986, 0.000, 0.000]],
		[0.757, [0.986, 0.000, 0.000]],
		[0.759, [0.978, 0.000, 0.000]],
		[0.761, [0.978, 0.000, 0.000]],
		[0.763, [0.971, 0.000, 0.000]],
		[0.765, [0.971, 0.000, 0.000]],
		[0.767, [0.963, 0.000, 0.000]],
		[0.769, [0.963, 0.000, 0.000]],
		[0.771, [0.955, 0.000, 0.000]],
		[0.773, [0.955, 0.000, 0.000]],
		[0.775, [0.947, 0.000, 0.000]],
		[0.777, [0.947, 0.000, 0.000]],
		[0.779, [0.939, 0.000, 0.000]],
		[0.781, [0.939, 0.000, 0.000]],
		[0.783, [0.931, 0.000, 0.000]],
		[0.785, [0.931, 0.000, 0.000]],
		[0.787, [0.924, 0.000, 0.000]],
		[0.789, [0.924, 0.000, 0.000]],
		[0.791, [0.916, 0.000, 0.000]],
		[0.793, [0.916, 0.000, 0.000]],
		[0.795, [0.908, 0.000, 0.000]],
		[0.796, [0.908, 0.000, 0.000]],
		[0.798, [0.900, 0.000, 0.000]],
		[0.800, [0.900, 0.000, 0.000]],
		[0.802, [0.892, 0.000, 0.000]],
		[0.804, [0.892, 0.000, 0.000]],
		[0.806, [0.884, 0.000, 0.000]],
		[0.808, [0.884, 0.000, 0.000]],
		[0.810, [0.876, 0.000, 0.000]],
		[0.812, [0.876, 0.000, 0.000]],
		[0.814, [0.869, 0.000, 0.000]],
		[0.816, [0.869, 0.000, 0.000]],
		[0.818, [0.861, 0.000, 0.000]],
		[0.820, [0.861, 0.000, 0.000]],
		[0.822, [0.853, 0.000, 0.000]],
		[0.824, [0.853, 0.000, 0.000]],
		[0.826, [0.845, 0.000, 0.000]],
		[0.828, [0.845, 0.000, 0.000]],
		[0.830, [0.837, 0.000, 0.000]],
		[0.832, [0.837, 0.000, 0.000]],
		[0.834, [0.829, 0.000, 0.000]],
		[0.836, [0.829, 0.000, 0.000]],
		[0.838, [0.822, 0.000, 0.000]],
		[0.840, [0.822, 0.000, 0.000]],
		[0.841, [0.814, 0.000, 0.000]],
		[0.843, [0.814, 0.000, 0.000]],
		[0.845, [0.806, 0.000, 0.000]],
		[0.847, [0.806, 0.000, 0.000]],
		[0.849, [0.798, 0.000, 0.000]],
		[0.851, [0.798, 0.000, 0.000]],
		[0.853, [0.790, 0.000, 0.000]],
		[0.855, [0.790, 0.000, 0.000]],
		[0.857, [0.782, 0.000, 0.000]],
		[0.859, [0.782, 0.000, 0.000]],
		[0.861, [0.775, 0.000, 0.000]],
		[0.863, [0.775, 0.000, 0.000]],
		[0.865, [0.767, 0.000, 0.000]],
		[0.867, [0.767, 0.000, 0.000]],
		[0.869, [0.759, 0.000, 0.000]],
		[0.871, [0.759, 0.000, 0.000]],
		[0.873, [0.751, 0.000, 0.000]],
		[0.875, [0.751, 0.000, 0.000]],
		[0.877, [0.743, 0.000, 0.000]],
		[0.879, [0.743, 0.000, 0.000]],
		[0.881, [0.735, 0.000, 0.000]],
		[0.883, [0.735, 0.000, 0.000]],
		[0.885, [0.727, 0.000, 0.000]],
		[0.886, [0.727, 0.000, 0.000]],
		[0.888, [0.720, 0.000, 0.000]],
		[0.890, [0.720, 0.000, 0.000]],
		[0.892, [0.712, 0.000, 0.000]],
		[0.894, [0.712, 0.000, 0.000]],
		[0.896, [0.704, 0.000, 0.000]],
		[0.898, [0.704, 0.000, 0.000]],
		[0.900, [0.696, 0.000, 0.000]],
		[0.902, [0.696, 0.000, 0.000]],
		[0.904, [0.688, 0.000, 0.000]],
		[0.906, [0.688, 0.000, 0.000]],
		[0.908, [0.680, 0.000, 0.000]],
		[0.910, [0.680, 0.000, 0.000]],
		[0.912, [0.673, 0.000, 0.000]],
		[0.914, [0.673, 0.000, 0.000]],
		[0.916, [0.665, 0.000, 0.000]],
		[0.918, [0.665, 0.000, 0.000]],
		[0.920, [0.657, 0.000, 0.000]],
		[0.922, [0.657, 0.000, 0.000]],
		[0.924, [0.649, 0.000, 0.000]],
		[0.926, [0.649, 0.000, 0.000]],
		[0.928, [0.641, 0.000, 0.000]],
		[0.930, [0.641, 0.000, 0.000]],
		[0.932, [0.633, 0.000, 0.000]],
		[0.933, [0.633, 0.000, 0.000]],
		[0.935, [0.625, 0.000, 0.000]],
		[0.937, [0.625, 0.000, 0.000]],
		[0.939, [0.618, 0.000, 0.000]],
		[0.941, [0.618, 0.000, 0.000]],
		[0.943, [0.610, 0.000, 0.000]],
		[0.945, [0.610, 0.000, 0.000]],
		[0.947, [0.602, 0.000, 0.000]],
		[0.949, [0.602, 0.000, 0.000]],
		[0.951, [0.594, 0.000, 0.000]],
		[0.953, [0.594, 0.000, 0.000]],
		[0.955, [0.586, 0.000, 0.000]],
		[0.957, [0.586, 0.000, 0.000]],
		[0.959, [0.578, 0.000, 0.000]],
		[0.961, [0.578, 0.000, 0.000]],
		[0.963, [0.571, 0.000, 0.000]],
		[0.965, [0.571, 0.000, 0.000]],
		[0.967, [0.563, 0.000, 0.000]],
		[0.969, [0.563, 0.000, 0.000]],
		[0.971, [0.555, 0.000, 0.000]],
		[0.973, [0.555, 0.000, 0.000]],
		[0.975, [0.547, 0.000, 0.000]],
		[0.977, [0.547, 0.000, 0.000]],
		[0.978, [0.539, 0.000, 0.000]],
		[0.980, [0.539, 0.000, 0.000]],
		[0.982, [0.531, 0.000, 0.000]],
		[0.984, [0.531, 0.000, 0.000]],
		[0.986, [0.524, 0.000, 0.000]],
		[0.988, [0.524, 0.000, 0.000]],
		[0.990, [0.516, 0.000, 0.000]],
		[0.992, [0.516, 0.000, 0.000]],
		[0.994, [0.508, 0.000, 0.000]],
		[0.996, [0.508, 0.000, 0.000]],
		[0.998, [0.500, 0.000, 0.000]],
		[1.000, [0.500, 0.000, 0.000]]
];
*/
