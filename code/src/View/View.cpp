#include "View.h"

#include <igl/unproject_in_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/hsv_to_rgb.h>

//#include <imgui/imgui.h>

#include "GlobalSettings.h"

#include "GeodesicWalker.h"

#include "PatchModel.h"
#include "MeshController.h"

#include "CurveHelper.h"
#include "ViewUtils.h"
#include "Utils/Colors.h"
#include "Utils.h"
#include "Logger.h"

#include "TempUtils.h"

using namespace std;


void View::initialize()
{
	// Attach a menu plugin
	data_model.viewer.plugins.push_back(&data_model.menu);

	data_model.viewer.core().background_color << 1.0f, 1.0f, 1.0f, 1.0f;
	data_model.viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

	//data_model.viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
	//data_model.viewer.core().rotation_type = igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL;
	//data_model.viewer.resize(1600, 1400);
	//data_model.viewer.resize(1920*0.65, 1080*0.9); //for recording
	data_model.viewer.resize(1920*0.9, 1080); //for rotation recording

	//debug view settings
	data_model.viewer.selected_data_index = view_model.debug_view_index;
#if show_coordinate_indicator
	const auto coordinate_indicator = Eigen::MatrixXd::Identity(3, 3);
	data_model.viewer.data().add_edges(Eigen::MatrixXd::Zero(3, 3), coordinate_indicator*0.2, coordinate_indicator);
#endif
	data_model.viewer.data().point_size = 7.0;
	data_model.viewer.data().line_width = 1.5;
	//data_model.viewer.data().line_width = 0.5;

	view_model.target_developable_view_settings.show_mesh = false;
	view_model.developable_parts_view_settings.view_index = -1;

	//add gauss map 
	meshhelper::add_gauss_map(data_model, view_model);


	//if (data_model.clustered_feature_vertices.size() > 0)
	//	segment_selection.segment_index = data_model.clustered_feature_vertices.size() - 1;

	view_model.general_patch_view_settings.view_index = -8;
	
	//view_model.patch_color = Eigen::RowVector3d(176.0/255.0, 168.0/255.0, 158.0/255.0); //taupe 
	view_model.patch_color = Eigen::RowVector3d(179.0/255.0, 176.0/255.0, 174.0/255.0); //taupe V 70%
	//view_model.patch_color = Eigen::RowVector3d(154.0/255.0, 154.0/255.0, 154.0/255.0); //gray V 60%
	//view_model.patch_color = Eigen::RowVector3d(172.0 / 255.0, 187.0 / 255.0, 198.0 / 255.0); 
	//view_model.patch_color = Colors::WHITE;
	view_model.general_patch_view_settings.standard_fill_color = view_model.patch_color;

	////set initial geodesics direction
	//if (view_model.do_trace_kmin && view_model.do_trace_kmax)
	//	data_model.geodesics_direction = GeodesicsDirection::K_BOTH;
	//else if (view_model.do_trace_kmin)
	//	data_model.geodesics_direction = GeodesicsDirection::K_MIN;
	//else if (view_model.do_trace_kmax)
	//	data_model.geodesics_direction = GeodesicsDirection::K_MAX;

	view_model.timer.start();
}

bool View::update_view()
{
	if (!is_initialized)
	{
		style_viewer();
		is_initialized = true;
	}

	//if (data_model.do_rebuild_geodesics)
	//	optimization.geodesics_controller->build_geodesics();
	//if (data_model.do_recluster_geodesics)
	//	optimization.geodesics_controller->cluster_geodesics();


	//data_model.viewer.selected_data_index = data_model.debug_view_index;
	data_model.viewer.data_list[view_model.debug_view_index].clear();

	//render selection point
	if (data_model.pre_selected_vertex >= 0 && data_model.target.V.rows())
		data_model.viewer.data().add_points(data_model.target.V.row(data_model.pre_selected_vertex), Colors::GRAY_DARK);

	data_model.viewer.data().add_points(view_model.selected_points, Colors::YELLOW);
	//if (view_model.do_update_points)
	//{
	//	data_model.viewer.data().add_points(view_model.selected_points, Colors::YELLOW);
	//	view_model.do_update_points = false;
	//}


	view_updater.update_gauss_map_view();


	//TODO update overlay visualization
	//TODO make list of ViewSettings such that can loop through, store of all visible meshes the original fill, set to the visualization and then set back


	////update target mesh
	//view_updater.update_target_view(view_model.target_view_settings, data_model.target);

	if (view_model.is_rotating)
	{
		Eigen::Vector3d rot_axis = Eigen::Vector3d::UnitY();
		switch (view_model.rotation_axis)
		{
		case 1:
			rot_axis = Eigen::Vector3d::UnitZ();
			break;
		case 2:
			rot_axis = Eigen::Vector3d::UnitX();
			break;
		}

		float animation_speed = .05 + view_model.animation_speed_flag * .05;

		Eigen::Matrix3d rot = Eigen::AngleAxisd(view_model.timer.getElapsedTimeInSec() * animation_speed * M_PI, rot_axis).toRotationMatrix();

		if (view_model.mid.size() < 3)
			view_model.mid = get_centroid(data_model.target.V);

		for (int i = 0; i < data_model.target.V.rows(); ++i)
		{
			data_model.target.V.row(i) = view_model.mid + (data_model.target_original.V.row(i) - view_model.mid) * rot;
		}

		data_model.viewer.selected_data_index = view_model.target_view_settings.view_index;
		data_model.viewer.data().set_vertices(data_model.target.V);
		data_model.viewer.data().compute_normals();
	}

	//update target mesh
	view_updater.update_target_view(view_model.target_view_settings, data_model.target);

	view_model.target_view_settings.has_changed = true;
	view_updater.update_target_fill();
	   	  
	//update developable mesh
	if (view_model.target_developable_view_settings.view_index != -1 && view_model.target_developable_view_settings.show_mesh && view_model.target_developable_view_settings.has_changed)
	{
		data_model.viewer.selected_data_index = view_model.target_developable_view_settings.view_index;
		data_model.viewer.data().clear();
		data_model.viewer.data().set_mesh(data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F);
		//data_model.viewer.data().set_mesh(data_model.target_developable.V, data_model.target_developable.F);
		view_model.target_developable_view_settings.show_faces = false;
	}

	view_updater.update_target_view(view_model.target_developable_view_settings, data_model.developable_model.non_developable);
	view_updater.update_mesh_fill(view_model.target_developable_view_settings, data_model.developable_model.non_developable);
	/*
	if (view_model.target_developable_view_settings.show_mesh && view_model.target_developable_view_settings.has_changed)
	{
		if (data_model.vertex_patch_assignment.size() == data_model.target_developable.V.rows())
		{
			data_model.viewer.selected_data_index = view_model.target_developable_view_settings.view_index;
			view_updater.update_mesh_color_assignment_fill(data_model.target_developable, data_model.vertex_patch_assignment, data_model.patches.size());
			view_model.target_developable_view_settings.has_changed = false;

			//reset developables parts
			data_model.developable_parts.clear();

			for (int i = 0; i < data_model.developable_parts_view_indices.size(); i++)
				data_model.viewer.data_list[data_model.developable_parts_view_indices[i]].clear();

			data_model.developable_parts_view_indices.clear();
		}
		else if (data_model.vertex_patch_assignment.size() == data_model.target_developable.F.rows())
		{

		}
	}
	*/

	//TODO update developable parts
	//view_updater.update_developables_parts_view();


	//update patches
	view_updater.update_patches_view();


	view_updater.update_geodesics_init_view();


	//update debug items	
	view_updater.update_debug_view();
	view_updater.update_feature_view(segment_selection);
	view_updater.update_debug_geodesics_view();

	/*
	if (data_model.vertex_labels.size() > 0)
	{
		//Eigen::MatrixXd points;
		//index_to_value(data_model.vertex_ruled_assignment, data_model.target.V, points);
		//data_model.viewer.data().add_points(points, Colors::RED);

		int number_labels = data_model.vertex_labels.maxCoeff();
		if (number_labels < 1)
			return false;

		vector<vector<int>> flat_areas(number_labels);
		const double hue_step = 360 / number_labels;



		for (int i = 0; i < data_model.vertex_labels.rows(); i++)
		{
			int label = data_model.vertex_labels(i);
			if (label >= 1)
				flat_areas[label - 1].push_back(i);
		}

		for (int i = 0; i < flat_areas.size(); i++)
		{
			vector<int> area = flat_areas[i];

			double hue = (i+1) * hue_step;
			Eigen::RowVector3d hsv(hue, 1, 1);
			Eigen::RowVector3d rgb;
			igl::hsv_to_rgb(hsv, rgb);

			Eigen::MatrixXd points;
			index_to_value(area, data_model.target.V, points);
			data_model.viewer.data().add_points(points, rgb);
		}
	}
	*/

	return false;
}


void View::update_menu()
{
	//styles: https://github.com/ocornut/imgui/issues/707 (the post from 22.8.2018)
	//toggle button: https://github.com/ocornut/imgui/issues/1537


	menu_updater.update_editor_menu();

	ImGui::Spacing();
	menu_updater.update_mode_menu();

	ImGui::Spacing();
	if (ImGui::CollapsingHeader("Gauss map"))
	{
		ImGui::Indent();


		ImGui::Text("apply gauss map to:");

		if (ImGui::RadioButton("target", &view_model.choice_gauss_map, 1))
			view_model.gauss_mapped_mesh = &data_model.target;
		ImGui::SameLine();
		if (ImGui::RadioButton("result", &view_model.choice_gauss_map, 2))
			view_model.gauss_mapped_mesh = &data_model.developable_model.result;
		ImGui::SameLine();
		if (ImGui::RadioButton("current patch", &view_model.choice_gauss_map, 3))
		{
			if (view_model.selected_patch_index >= 0 && view_model.selected_patch_index < data_model.patches.size())
				view_model.gauss_mapped_mesh = &(data_model.patches[view_model.selected_patch_index]->wrapper);
		}

		if (ImGui::Checkbox("show gauss map", &view_model.is_gauss_map_visible))
			view_model.do_update_gauss_map = true;

		ImGui::Unindent();
	}


	ImGui::Spacing();
	if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Indent();

		menu_updater.update_target_mesh_menu(optimization);
		menu_updater.update_wrapper_mesh_menu();

		menu_updater.update_surface_analysis_menu();
		ImGui::Unindent();
	}

	ImGui::Spacing();
	//menu_updater.update_features_menu(segment_selection);
	//menu_updater.update_geodesics_menu();

	menu_updater.update_ruled_geodesics_menu(*(optimization.geodesics_controller));
	menu_updater.update_optimization_menu(optimization);
}

void View::update_debug_menu()
{
	menu_updater.update_debug_window_menu(optimization);
	menu_updater.update_geodesics_menu();
}

bool View::callback_key_down(unsigned int key)
{
	switch (key)
	{
	case ' ':
		data_model.is_optimizing = !data_model.is_optimizing;
		break;
	case 341: //ctrl
		view_model.current_interaction_mode = InteractionMode::Select;
		bool is_selecting_features = view_model.current_interaction_mode == InteractionMode::Select && view_model.current_selection_mode == SelectionMode::DefineFeatures;
		segment_selection.toggle_segmentation(is_selecting_features);
		break;
	//case 340: //shift
	//	break;
	}

	//write_log(4) << "key DOWN::view_model.current_interaction_mode: " << view_model.current_interaction_mode << endl << endl;
	return false;
}

bool View::callback_key_up(unsigned int key)
{
	switch (key)
	{
	case 341: //ctrl
	{
		write_log(3) << "vertex index: " << data_model.pre_selected_vertex << endl;
		view_model.current_interaction_mode = InteractionMode::Navigate;
		segment_selection.toggle_segmentation(false);
		break;
	}
	case 'I':
	{
		optimization.initialize();
		break;
	}
	case 'F':
	{	//view_model.current_interaction_mode = InteractionMode::Navigate;
		segment_selection.finalize_segment();
		break;
	}
	// keyboard shortcuts for rotating animation
	case GLFW_KEY_M:
	{
		view_model.timer.start();
		data_model.target.V = data_model.target_original.V;
		view_model.target_view_settings.has_changed = true;
		break;
	}
	case GLFW_KEY_N:
	{
		view_model.timer.start();
		data_model.target.V = data_model.target_original.V;
		view_model.rotation_axis = (view_model.rotation_axis + 1) % 3;
		view_model.target_view_settings.has_changed = true;
		break;
	}
	case GLFW_KEY_B:
	{
		view_model.is_rotating = !view_model.is_rotating;
		view_model.target_view_settings.has_changed = true;
		break;
	}
	case ',':
		view_model.timer.start();
		view_model.animation_speed_flag = (view_model.animation_speed_flag + 1) % 5;
		view_model.target_view_settings.has_changed = true; 
		write_log(0) << "animation speed:: " << view_model.animation_speed_flag << ". Can be changed from 0 to 4. " << endl;
		break;
	}

	//write_log(4) << "key UP::view_model.current_interaction_mode: " << view_model.current_interaction_mode << endl << endl;
	return false;
}

bool View::callback_mouse_move(int mouse_x, int mouse_y)
{
	if (view_model.current_interaction_mode == InteractionMode::Navigate)
		return false;

	data_model.pre_selected_vertex = get_vertex_from_screen(data_model.viewer, mouse_x, mouse_y, data_model.target.V, data_model.target.F);

	if (view_model.current_interaction_mode == InteractionMode::Select && view_model.current_selection_mode == SelectionMode::DefineFeatures)
		segment_selection.find_edge_path();


	if (data_model.target.surface_features().are_built)
	{
		int vi = data_model.pre_selected_vertex;
		write_log(3) << "K(" << vi << ")  = " << data_model.target.surface_features().K(vi) << endl;
		write_log(3) << "Kp(" << vi << ") = " << data_model.target.surface_features().Kp(vi) << " (min=" << data_model.target.surface_features().principal_k_min(vi) << " * max=" << data_model.target.surface_features().principal_k_max(vi) << ")" << endl;
	}


	return false;
}

bool View::callback_mouse_down(int button)
{
	view_model.is_mouse_down = true;
	return false;
}

bool View::callback_mouse_up(int button)
{
	view_model.is_mouse_down = false;

	if (data_model.viewer.mouse_mode == igl::opengl::glfw::Viewer::MouseMode::None) //likely that the menu was clicked
		return false;

	if (view_model.current_interaction_mode == InteractionMode::Navigate)
		return false;

	if (view_model.current_interaction_mode == InteractionMode::Select)
	{
		//if (view_model.current_selection_mode == SelectionMode::DefineFeatures)
		//{
		//	segment_selection.add_path();
		//}
		//else if (view_model.current_selection_mode == SelectionMode::GeodesicsToAll && data_model.pre_selected_vertex >= 0)
		//{
		//	view_model.selected_geodesic_source_vertex = data_model.pre_selected_vertex;
		//	view_model.has_geodesic_selection_changed = true;
		//	//view_model.target_view_settings.has_changed = true;
		//}
		//else if (view_model.current_selection_mode == SelectionMode::SpecificGeodesic && data_model.pre_selected_vertex >= 0)
		//{
		//	if (view_model.paired_geodesic_sources.size() == view_model.paired_geodesic_targets.size())
		//		view_model.paired_geodesic_sources.push_back(data_model.pre_selected_vertex);
		//	else 
		//		view_model.paired_geodesic_targets.push_back(data_model.pre_selected_vertex);

		//	view_model.selected_vertex = view_model.paired_geodesic_targets.size() - 1;
		//	view_model.has_geodesic_selection_changed = true;
		//	//view_model.target_view_settings.has_changed = true;
		//}
		//else 

		//if (view_model.current_selection_mode == SelectionMode::WalkingGeodesic && data_model.pre_selected_vertex >= 0)
		//{
		if (view_model.current_selection_mode == SelectionMode::VertexPicking && data_model.pre_selected_vertex >= 0)
		{
			if (data_model.pre_selected_vertex >= data_model.target.V.rows())
				return false;

			if (view_model.is_selecting_debug_geodesic)
			{
				if (view_model.do_render_debug_walking_geodesics)
				{
					/*
					GeodesicWalker walker(data_model.target, data_model.geodesics_stopping);

					Eigen::RowVector3d k_max = data_model.target.surface_features().principal_k_max_direction.row(data_model.pre_selected_vertex);
					Eigen::RowVector3d k_min = data_model.target.surface_features().principal_k_min_direction.row(data_model.pre_selected_vertex);
					Eigen::RowVector3d axis = k_max.cross(k_min).normalized();
					int angle_deg = data_model.geodesics_directions[0];

					Eigen::RowVector3d direction;

					if (angle_deg == 0)
					{
						direction = k_max;
					}
					else if (angle_deg == 90)
					{
						direction = k_min;
					}
					else
					{
						double angle = angle_deg / 180.0 * igl::PI;

						Eigen::Matrix3d R;
						R = Eigen::AngleAxisd(angle, axis);
						direction = k_max * R;
					}
					write_log(4) << "axis: " << axis << ", angle_deg: " << angle_deg << ", rotated direction: " << direction << endl;


					Eigen::MatrixXd current_geodesic;
					double current_length;
					walker.get_geodesic_at(data_model.pre_selected_vertex, direction, current_geodesic, current_length);
					*/

					//get geodesic
					GeodesicWalker walker(data_model.target, data_model.geodesics_stopping);
					Eigen::MatrixXd geodesic;
					double length = 0;

					temp::get_geodesic(data_model.pre_selected_vertex, data_model.target, walker, data_model.geodesics_directions, geodesic, length);
					view_model.debug_walking_geodesics.push_back(geodesic);
					write_log(4) << "Geodesic at v" << data_model.pre_selected_vertex << ": length = " << length << ", rows = " << geodesic.rows() << ", is flat ? " << boolalpha << curve::is_flat(geodesic) << endl << endl;
				}
				else
				{
					view_model.walking_geodesics.push_back(data_model.pre_selected_vertex);
				}

				view_model.do_update_geodesics = true;
				view_model.has_geodesic_selection_changed = true;
			}
			else if (view_model.is_debugging_ruled_surfaces)
			{
				view_model.selected_debug_ruled_surface_index = data_model.pre_selected_vertex;
				view_model.do_update_debug_ruled_surfaces = true;

				/*
				//get geodesic
				GeodesicWalker walker(data_model.target, data_model.geodesics_stopping);

				Eigen::RowVector3d k_max = data_model.target.surface_features().principal_k_max_direction.row(view_model.selected_debug_ruled_surface_index);
				Eigen::RowVector3d k_min = data_model.target.surface_features().principal_k_min_direction.row(view_model.selected_debug_ruled_surface_index);
				Eigen::RowVector3d axis = k_max.cross(k_min).normalized();
				int angle_deg = data_model.geodesics_directions[0];

				Eigen::RowVector3d direction;

				if (angle_deg == 0)
				{
					direction = k_max;
				}
				else if (angle_deg == 90)
				{
					direction = k_min;
				}
				else
				{
					double angle = angle_deg / 180.0 * igl::PI;

					Eigen::Matrix3d R;
					R = Eigen::AngleAxisd(angle, axis);
					direction = k_max * R;
				}
				write_log(4) << "axis: " << axis << ", angle_deg: " << angle_deg << ", rotated direction: " << direction << endl;

				Eigen::MatrixXd current_geodesic;
				double current_length;
				walker.get_geodesic_at(view_model.selected_debug_ruled_surface_index, direction, current_geodesic, current_length);


				//get ruled surface
				RulingsType rulings_type;
				switch (view_model.use_debug_rulings_type)
				{
					case(0):
					{
						bool is_flat = curve::is_flat(current_geodesic);
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
				}

				if(view_model.debug_ruled_surface != NULL)
					delete view_model.debug_ruled_surface;

				view_model.debug_ruled_surface = new RuledDevelopableSurface(data_model.target, current_geodesic, view_model.selected_debug_ruled_surface_index, rulings_type);
				view_model.debug_ruled_surface->create(data_model.ruled_width);

				view_model.do_update_debug_ruled_surfaces = true;
				//view_model.has_debug_ruled_surface_changed = true;
				*/
			}
		}
	}

	return false;
}

int View::get_vertex_from_screen(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	// Cast a ray in the view direction starting from the mouse position
	double x = mouse_x;
	double y = viewer.core().viewport(3) - mouse_y;

	Eigen::RowVector3d pt;
	int vi = -1;
	std::vector<igl::Hit> hits;

	/*
	Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;
	igl::unproject_in_mesh(Eigen::Vector2f(x,y),
						   modelview,
						   viewer.core().proj,
						   viewer.core().viewport,
						   ei,pt,hits);
	*/

	igl::unproject_in_mesh(Eigen::Vector2f(x, y), viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, pt, hits);

	if (hits.size() > 0)
	{
		int fi = hits[0].id;
		Eigen::RowVector3d bc;
		bc << 1.0 - hits[0].u - hits[0].v, hits[0].u, hits[0].v;

		auto coeff = bc.maxCoeff(&vi);
		//write_log(0) << endl << "get_vertex_from_screen: hits.size() = " << hits.size() << ", max coeff = " << coeff << endl;

		vi = F(fi, vi);
	}
	return vi;
}


void View::style_viewer()
{
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsLight();
	ImGuiStyle& style = ImGui::GetStyle();
	style.FrameRounding = 0.0f;
}
