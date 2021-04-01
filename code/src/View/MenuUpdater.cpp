#include "MenuUpdater.h"

#include <stdio.h> 
#include <time.h>
#include <experimental/filesystem>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <igl/boundary_loop.h>
#include <igl/random_points_on_mesh.h>
#include <igl/upsample.h>

#include <igl/readOBJ.h>
#include <igl/per_corner_normals.h>

#include "OptimizationController.h"
#include "Serializer.h"
#include "MeshController.h"
#include "PatchModel.h"
#include "GeodesicWalker.h"
#include "CurveHelper.h"
#include "CoordinateConverter.h"
#include "Rulings/RuledDevelopableSurface.h"

#include "Utils.h"
#include "Logger.h"

using namespace std;

std::string get_time_string()
{
	time_t rawtime;
	time(&rawtime);
	struct tm* timeinfo = localtime(&rawtime);

	char buffer[80];
	strftime(buffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

	return buffer;
}

void MenuUpdater::update_debug_window_menu(OptimizationController& optimization)
{
	int loglevel = LOG_LEVEL;
	if (ImGui::InputInt("debug output level", &loglevel)) 
		LOG_LEVEL = loglevel < 1 ? 1 : loglevel;

	if (ImGui::InputInt("optimization output level", &DataModel::log_level_optimization)) 
		DataModel::log_level_optimization = max(1, DataModel::log_level_optimization);

	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	if (ImGui::Button("implicitly smooth target"))
		meshhelper::implicit_smooth_mesh(data_model.target.V, data_model.target.F, data_model.target.V, data_model.target.F, data_model.target.normals_vertices);

	/*
	ImGui::Checkbox("show debug items", &view_model.show_debug_items);

	ImGui::Separator();
	ImGui::Spacing();

	ImGui::Text("INFO:");

	double coverage_info = data_model.target.coverage_distances_V.rows() ? data_model.target.coverage_distances_V.sum() : -1;
	ImGui::Text(("current coverage: " + to_string(coverage_info)).c_str());

	if (data_model.current_candidates_index >= 0 && data_model.current_candidates_index < data_model.feature_geodesics_candidates.size())
	{
		auto& geodesics = data_model.feature_geodesics_candidates[data_model.current_candidates_index];
		ImGui::Text(("feature neighborhood: " + to_string(data_model.current_candidates_index)).c_str());
		ImGui::Text(("geodesic: " + to_string(geodesics.current_geodesic_index)).c_str());
		ImGui::Text(("merit: " + to_string(geodesics.current_merit)).c_str());

		ImGui::Spacing();
		double max_merit_info = geodesics.geodesics_merit.rows() ? geodesics.geodesics_merit.maxCoeff() : -1;
		double max_cover_info = geodesics.geodesic_features.rows() ? geodesics.geodesic_features.col(data_model.feature_settings.dynamic_geodesic_features_index).maxCoeff() : -1;
		ImGui::Text(("max merit: " + to_string(max_merit_info)).c_str());
		ImGui::Text(("max cover: " + to_string(max_cover_info)).c_str());
		//int index;
		//double max_cover = geodesics.geodesic_features.col(data_model.feature_settings.number_geodesic_features - 1).maxCoeff(index);
		//ImGui::Text(("max cover absolute: " + to_string(max_cover)).c_str());
		//ImGui::Text(("max cover relative: " + to_string(max_cover)).c_str());
	}
	

	ImGui::Spacing();
	ImGui::Text("SETTINGS:");

	bool has_geodesics_weight_changed = false;
	ImGui::Text("weight for geodesic candidates selection:");

	if (ImGui::SliderFloat("gauss curvature", &data_model.feature_settings.geodesic_feature_weights[0], -1.0, 1.0))
		has_geodesics_weight_changed = true;
	if (ImGui::SliderFloat("length", &data_model.feature_settings.geodesic_feature_weights[1], -1.0, 1.0))
		has_geodesics_weight_changed = true;
	if (ImGui::SliderFloat("curvature", &data_model.feature_settings.geodesic_feature_weights[2], -1.0, 1.0))
		has_geodesics_weight_changed = true;
	if (ImGui::SliderFloat("coverage", &data_model.feature_settings.geodesic_feature_weights[3], -5.0, 5.0))
		has_geodesics_weight_changed = true;

	ImGui::InputFloat("geodesics coverage threshold", &data_model.optimization_settings.geodesics_coverage_threshold);
	ImGui::InputFloat("wrapper width scale", &data_model.wrapper_width_scale);


	ImGui::Separator();
	ImGui::Spacing();
	
	ImGui::Text("VISUALIZE:");

	if(ImGui::Checkbox("show all geodesic candidates", &view_model.show_all_geodesics_candidates))
		view_model.do_update_geodesics = true;

	if (ImGui::InputInt("visualized cluster", &view_model.selected_feature_cluster))
	{
		view_model.selected_feature_cluster = clip(view_model.selected_feature_cluster, -1, data_model.feature_connectivity_graph.size() - 1);
		has_geodesics_weight_changed = true;
	}

	ImGui::Checkbox("show feature clusters", &view_model.show_feature_cluster_graph);
	if(ImGui::Checkbox("show selected candidates", &view_model.show_geodesics_candidates))
		view_model.do_update_geodesics = true;

	ImGui::Indent();

		if(ImGui::Checkbox("show labels##geodesics", &view_model.show_geodesics_candidates_labels))
			view_model.do_update_geodesics = true;

	ImGui::Unindent();

	//if (has_geodesics_weight_changed)
	//{
	//	const int selected_feature_index = view_model.selected_feature_cluster;
	//	if (selected_feature_index < 0 || selected_feature_index >= data_model.feature_geodesics_candidates.size())
	//		return;

	//	const auto& geodesic_candidates = data_model.feature_geodesics_candidates[selected_feature_index];
	//	if (geodesic_candidates.number_geodesics < 1)
	//		return;

	//	view_model.merit_max_index = geodesic_candidates.geodesics_merit.maxCoeff(&view_model.merit_max_index);
	//	optimization.geodesics_controller->update_geodesics_merit(selected_feature_index);

	//	view_model.do_update_geodesics = true;
	//}

	ImGui::Text(("best index preview: " + to_string(view_model.merit_max_index)).c_str());

	if (ImGui::Button("print candidate features"))
	{
		if (view_model.selected_feature_cluster < 0 || view_model.selected_feature_cluster > data_model.feature_geodesics_candidates.size())
			return;
		write_log(0) << "geodesic_features: " << endl << data_model.feature_geodesics_candidates[view_model.selected_feature_cluster].geodesic_features << endl;
		write_log(0) << "geodesics_merit: " << endl << data_model.feature_geodesics_candidates[view_model.selected_feature_cluster].geodesics_merit << endl;
	}
	*/
}

void MenuUpdater::update_editor_menu()
{
	if (!ImGui::CollapsingHeader("Editor"))
		return;


	ImGui::Text("workspace: "); 


	if (ImGui::Button("load##Workspace"))
	{
		std::string fname = igl::file_dialog_open();
		if (fname.length() == 0)
			return;

		deserialize(fname, data_model, view_model, true);
	}
	ImGui::SameLine();
	if (ImGui::Button("save##Workspace"))
	{
		std::string fname = igl::file_dialog_save();
		if (fname.length() == 0)
			return;

		serialize(fname, data_model);
	}
	ImGui::SameLine();
	if (ImGui::Button("save (default location)##Workspace"))
	{
		if (data_model.target_filename.length() == 0)
			return;

		string filename = data_model.target_filename;
		size_t last_dot = filename.rfind('.');
		if (last_dot != std::string::npos)
			filename.erase(last_dot, filename.length() - last_dot);

		filename += "__scene";

		//filename += "__";
		//if (data_model.do_persist_features && data_model.clustered_feature_vertices.size() > 0)
		//	filename += "_features";
		//if (data_model.do_persist_graph && data_model.feature_connectivity_graph.size() > 0)
		//	filename += "_graph";
		//if (data_model.do_persist_geodesic_candidates && data_model.feature_geodesics_candidates.size() > 0)
		//	filename += "_geodesics";

		serialize(filename, data_model);
	}


	ImGui::Separator();

	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;
	auto button_size = ImVec2((w - p) / 3.f, 0);
	auto label_size = ImVec2((w - p) / 3.f, 0);

	Eigen::Vector4f wireframe_color(Colors::GRAY_MID(0), Colors::GRAY_MID(1), Colors::GRAY_MID(2), 1.0f);

	ImGui::Text("mesh (.obj): "); 
	if (ImGui::Button("load target", button_size))
	{
		std::string filename = igl::file_dialog_open();
		if (filename.length() == 0)
			return;

		//meshhelper::add_target(filename, data_model, view_model.target_view_settings.view_index, view_model.target_developable_view_settings.view_index);
		
		bool success = igl::readOBJ(filename, data_model.target.V, data_model.target.F);
		if (!success)
			write_log(1) << "error at loading target. (filename: " << filename << ")" << endl;
		
		meshhelper::add_mesh(data_model.viewer, data_model.target.V, data_model.target.F, false, true);
		data_model.target_original.V = data_model.target.V;
		data_model.target_original.F = data_model.target.F;

		Eigen::RowVectorXd dimensions;
		meshhelper::calculate_dimensions(data_model.target.V, dimensions);
		data_model.target_max_dimension = dimensions.maxCoeff();
		data_model.target_dimensions = dimensions;
		data_model.target_diagonal = dimensions.norm();

		//data_model.developable_model.non_developable = data_model.target;
		view_model.target_view_settings.view_index = data_model.viewer.selected_data_index;
		//data_model.viewer.data().clear();
		//data_model.viewer.data().set_mesh(data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F);
		data_model.viewer.data().point_size = 5.0;
		data_model.viewer.data().line_width = 0.5;
		data_model.viewer.data().line_color = wireframe_color;


		//data_model.viewer.data().face_based = true;
		//data_model.viewer.data().F_material_specular = Colors::BLACK;
		//data_model.viewer.data().V_material_specular = Colors::BLACK;
	}
	if (ImGui::Button("load result", button_size))
	{
		std::string filename_result = igl::file_dialog_open();
		if (filename_result.length() == 0)
			return;

		bool success = igl::readOBJ(filename_result, data_model.developable_model.concatenated.V, data_model.developable_model.concatenated.F);
		if (!success)
			write_log(1) << "error at loading result. (filename: " << filename_result << ")" << endl;

		meshhelper::add_mesh(data_model.viewer, data_model.developable_model.concatenated.V, data_model.developable_model.concatenated.F, false, false);

		data_model.developable_model.non_developable = data_model.target;
		view_model.target_developable_view_settings.view_index = data_model.viewer.selected_data_index;
		//data_model.viewer.data().clear();
		//data_model.viewer.data().set_mesh(data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F);
		data_model.viewer.data().point_size = 5.0;
		data_model.viewer.data().line_width = 0.5;
		data_model.viewer.data().line_color = wireframe_color;

		//data_model.viewer.data().face_based = true;

		//Eigen::MatrixXd CN;
		//igl::per_corner_normals(data_model.developable_model.concatenated.V, data_model.developable_model.concatenated.F, 20, CN);
		//data_model.viewer.data().set_normals(CN);

		data_model.viewer.data().F_material_specular = Colors::BLACK;
		data_model.viewer.data().V_material_specular = Colors::BLACK;

		data_model.result_coverage = new Coverage(data_model.target, &(data_model.optimization_settings.coverage_threshold));
		data_model.result_coverage->compute(data_model.developable_model.concatenated);

		Eigen::VectorXd vertex_error = data_model.result_coverage->distances();
		write_log(0) << std::fixed << vertex_error << endl;
		ofstream file("coverage_result.txt");
		file << std::fixed << vertex_error << endl;
		file.close();
	}

	ImGui::SameLine();
	if (ImGui::Button("save target", button_size))
	{
		std::string filename = igl::file_dialog_save();
		if (filename.length() == 0)
			return;

		size_t last_dot = filename.rfind('.');
		if (last_dot != std::string::npos)
			filename.erase(last_dot, filename.length() - last_dot);
		filename += ".obj";

		igl::writeOBJ(filename, data_model.target.V, data_model.target.F);
	}


	if (ImGui::Button("save result package"))
	{
		std::string extension = ".obj";

		//data_model.target_filename: everything after the models folder (eg "/91 models for region finding - DUPLICATES/cone_argyle_w30.obj")

		string filepath = data_model.target_filename; 
		size_t last_slash = filepath.rfind('/');
		if (last_slash != std::string::npos)
			filepath.erase(last_slash, filepath.length() - last_slash);


		string filename = data_model.target_filename;
		filename.erase(0, last_slash);
		
		size_t last_dot = filename.rfind('.');
		if (last_dot != std::string::npos)
			filename.erase(last_dot, filename.length() - last_dot);
		
		//if (filename[0] == '/')
		//	filename.erase(0, 0);
		size_t first_slash = filename.rfind('/');
		if (first_slash != std::string::npos)
			filename.erase(0, first_slash+1);


		string datetime = get_time_string();
		string path = filepath + "/" + datetime + "__" + filename + "/";
		
		bool success = false;
		if(!std::experimental::filesystem::exists(path))
			success = std::experimental::filesystem::create_directories(data_model.models_folder + path);


		serialize(path + filename + "__scene", data_model);

		path = data_model.models_folder + path + filename;

		if (data_model.patches.size() > 0)
		{
			meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);
			igl::writeOBJ(path + "__DOGs" + extension, data_model.concatenated_patches_V, data_model.concatenated_patches_F);

			for(int i = 0; i < data_model.patches.size(); i++)
				igl::writeOBJ(path + "__DOG" + to_string(i+1) + extension, data_model.patches[i]->wrapper.V, data_model.patches[i]->wrapper.F);
		}

		if (data_model.developable_model.developable_parts.size() > 0)
		{
			for (int i = 0; i < data_model.developable_model.developable_parts.size(); i++)
				igl::writeOBJ(path + "__parts" + to_string(i+1) + extension, data_model.developable_model.developable_parts[i].V, data_model.developable_model.developable_parts[i].F);
		}
		if (data_model.developable_model.developable_parts_2D.size() > 0)
		{
			for (int i = 0; i < data_model.developable_model.developable_parts_2D.size(); i++)
				igl::writeOBJ(path + "__parts_flat" + to_string(i+1) + extension, data_model.developable_model.developable_parts_2D[i].V, data_model.developable_model.developable_parts_2D[i].F);
		}

		if(data_model.target.V.rows() > 0)
			igl::writeOBJ(path + "__input" + extension, data_model.target.V, data_model.target.F);

		if(data_model.developable_model.result.V.rows() > 0)
			igl::writeOBJ(path + "__result" + extension, data_model.developable_model.result.V, data_model.developable_model.result.F);
		
		if(data_model.developable_model.concatenated.V.rows() > 0)
			igl::writeOBJ(path + "__concatenated" + extension, data_model.developable_model.concatenated.V, data_model.developable_model.concatenated.F);

		//save seam pairs
		ofstream seam_file(path + "__concatenated_seams.txt");
		for (int i = 0; i < data_model.developable_model.seam_pairs.rows(); i++)
		{
			Eigen::Vector2i pair = data_model.developable_model.seam_pairs.row(i);
			seam_file << to_string(pair(0)) << " " << to_string(pair(1)) << endl;
		}
		seam_file.close();
	}
	/*
	if (ImGui::Button("save result", button_size))
	{
		std::string filename = igl::file_dialog_save();
		if (filename.length() == 0)
			return;

		size_t last_dot = filename.rfind('.');
		if (last_dot != std::string::npos)
			filename.erase(last_dot, filename.length() - last_dot);
		filename += ".obj";

		igl::writeOBJ(filename, data_model.developable_model.result.V, data_model.developable_model.result.F);
	}
	if (ImGui::Button("save result parts", button_size))
	{
		if (data_model.developable_model.developable_parts.size() > 0)
		{
			std::string filename = igl::file_dialog_save();
			if (filename.length() == 0)
				return;

			size_t last_dot = filename.rfind('.');
			if (last_dot != std::string::npos)
				filename.erase(last_dot, filename.length() - last_dot);
			filename += ".obj";

			Mesh concatenated_developable_parts = meshhelper::concatenate_meshes(data_model.developable_model.developable_parts);
			igl::writeOBJ(filename, concatenated_developable_parts.V, concatenated_developable_parts.F);
		}
	}
	//ImGui::SameLine(0, p);
	if (ImGui::Button("save DOGs", button_size))
	{
		std::string filename = igl::file_dialog_save();
		if (filename.length() == 0)
			return;

		size_t last_dot = filename.rfind('.');
		//std::string extension = fname.substr(last_dot + 1);
		if (last_dot != std::string::npos)
			filename.erase(last_dot, filename.length() - last_dot);
		filename += ".obj";

		meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);
		igl::writeOBJ(filename, data_model.concatenated_patches_V, data_model.concatenated_patches_F);
	}
	*/


	ImGui::Separator();

	// Select rotation type
	int rotation_type = static_cast<int>(data_model.viewer.core().rotation_type);
	static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
	static bool orthographic = true;
	if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
	{
		using RT = igl::opengl::ViewerCore::RotationType;
		auto new_type = static_cast<RT>(rotation_type);
		if (new_type != data_model.viewer.core().rotation_type)
		{
			if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
			{
				trackball_angle = data_model.viewer.core().trackball_angle;
				orthographic = data_model.viewer.core().orthographic;
				data_model.viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
				data_model.viewer.core().orthographic = true;
			}
			else if (data_model.viewer.core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
			{
				data_model.viewer.core().trackball_angle = trackball_angle;
				data_model.viewer.core().orthographic = orthographic;
			}
			data_model.viewer.core().set_rotation_type(new_type);
		}
	}

	ImGui::Checkbox("Orthographic view", &(data_model.viewer.core().orthographic));
	ImGui::ColorEdit4("Background", data_model.viewer.core().background_color.data(),
		ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);

	if (ImGui::Button("center mesh", ImVec2(-1, 0)))
		data_model.viewer.core().align_camera_center(data_model.target.V, data_model.target.F);
}

void MenuUpdater::update_mode_menu()
{
	ImGui::Text("MODE:");

	int choice = view_model.current_interaction_mode;

	ImGui::RadioButton("navigate", &choice, InteractionMode::Navigate);
	ImGui::RadioButton("select", &choice, InteractionMode::Select);

	view_model.current_interaction_mode = static_cast<InteractionMode>(choice);


	//ImGui::Checkbox("show vertex", &view_model.do_show_vertex);
	ImGui::InputInt("show/select vertex", &data_model.pre_selected_vertex);

	const string V_range_info = "(" + to_string(0) + " - " + to_string(data_model.target.V.rows()) + ")";
	ImGui::Text(V_range_info.c_str());
}

void MenuUpdater::update_target_mesh_menu(OptimizationController& optimization)
{
	ImGui::Spacing();
	if (ImGui::CollapsingHeader("Target", ImGuiTreeNodeFlags_DefaultOpen))
	{
		mesh_view_options(view_model.target_view_settings);

		ImGui::Text("Show coverage:");

		ImGui::RadioButton("none##coverage", &view_model.current_coverage_choice, -1); ImGui::SameLine();
		ImGui::RadioButton("init##coverage", &view_model.current_coverage_choice, 0); ImGui::SameLine();
		ImGui::RadioButton("DOGs##coverage", &view_model.current_coverage_choice, 1); ImGui::SameLine();
		ImGui::RadioButton("result##coverage", &view_model.current_coverage_choice, 2);
			
		ImGui::Checkbox("show lines##coverage", &view_model.show_coverage); ImGui::SameLine();
		ImGui::Checkbox("show points##coverage", &view_model.show_coverage_points);
		ImGui::Checkbox("show holes##coverage", &view_model.show_coverage_holes);
		//ImGui::Checkbox("show coloring##coverage", &view_model.show_coverage_coloring);

		//Coverage* selected_coverage = data_model.patch_coverage;

		//if (view_model.current_coverage_choice == 0)
		//	selected_coverage = data_model.initialization_coverage;
		//else if (view_model.current_coverage_choice == 1)
		//	selected_coverage = data_model.patch_coverage;
		//else if (view_model.current_coverage_choice == 2)
		//	selected_coverage = data_model.result_coverage;

		//bool is_coverage_updated = false;
		//if (ImGui::Checkbox("show lines", &view_model.show_coverage))
		//{
		//	if (view_model.show_coverage)
		//	{
		//		data_model.patch_coverage->compute(meshhelper::get_wrapper_meshes(data_model.patches));
		//		is_coverage_updated = true;
		//	}
		//}
		//if (ImGui::Checkbox("show points", &view_model.show_coverage_points))
		//{
		//	if (view_model.show_coverage_points && !is_coverage_updated)
		//		data_model.patch_coverage->compute(meshhelper::get_wrapper_meshes(data_model.patches));
		//}

	}


	//ImGui::Spacing();
	if (ImGui::CollapsingHeader("Developable Mesh"))
	{
		/*
		if (data_model.developable_parts_view_indices.size() < 1)
			deactivate_elements();
		

		ImGui::Text("developable PARTS");

		mesh_view_options(view_model.developable_parts_view_settings);

		if (ImGui::InputInt("show developable part", &view_model.selected_developable_part))
		{
			view_model.selected_developable_part = clip(view_model.selected_developable_part, -1, data_model.developable_parts_view_indices.size() - 1);
			view_model.developable_parts_view_settings.has_changed = true;

			if(view_model.selected_developable_part == -1)
				view_model.developable_parts_view_settings.view_index = view_model.selected_developable_part;
			else
				view_model.developable_parts_view_settings.view_index = data_model.developable_parts_view_indices[view_model.selected_developable_part];
		}

		if (data_model.developable_parts_view_indices.size() < 1)
			activate_elements();
		*/


		ImGui::Text("global developable mesh");

		mesh_view_options(view_model.target_developable_view_settings);

		ImGui::Text("create RESULT");


		ImGui::Indent();
			ImGui::Text("RESULT optimization settings:");
			
			ImGui::Checkbox("use DEBUG mode? ", &data_model.developable_model.optimization_settings.use_debug_mode);
			ImGui::InputFloat("label_selection_smoothness", &data_model.developable_model.optimization_settings.label_selection_smoothness);
		
			ImGui::Text("weights:");
			ImGui::InputFloat("stiching", &data_model.developable_model.optimization_settings.weight_stiching_objective);
			ImGui::InputFloat("developable", &data_model.developable_model.optimization_settings.weight_developable_objective);
			ImGui::InputFloat("proximity", &data_model.developable_model.optimization_settings.weight_proximity_objective);
			ImGui::InputFloat("boundary", &data_model.developable_model.optimization_settings.weight_boundary_objective);

			ImGui::InputFloat("smoothness", &data_model.developable_model.optimization_settings.weight_smooth_objective);
			ImGui::InputFloat("proximity inner", &data_model.developable_model.optimization_settings.weight_proximity_inner_objective);
			ImGui::InputFloat("boundary seams", &data_model.developable_model.optimization_settings.weight_boundary_seam_objective);

		ImGui::Unindent();

		if (ImGui::Button("upsample result", ImVec2(-1, 0)))
		{
			igl::upsample(data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F, data_model.developable_model.non_developable.V, data_model.developable_model.non_developable.F, 1);
			view_model.target_developable_view_settings.has_changed = true;

			////reset developables parts
			//for (int i = 0; i < data_model.developable_parts_view_indices.size(); i++)
			//	data_model.viewer.data_list[data_model.developable_parts_view_indices[i]].clear();
		}
		if (ImGui::Button("reset result", ImVec2(-1, 0)))
		{
			optimization.reset_mesh_optimiztaion();
			
			data_model.developable_model.non_developable.V = data_model.target.V;
			data_model.developable_model.non_developable.F = data_model.target.F;
			data_model.developable_model.face_patch_assignment.clear();

			view_model.target_developable_view_settings.has_changed = true;
		}
		if (ImGui::Button("finish result", ImVec2(-1, 0)))
		{
			optimization.result_controller->finish_optimiztaion();
		}
		//if (ImGui::Button("re-initialize result", ImVec2(-1, 0)))
		//{
		//	optimization.reset_mesh_optimiztaion();
		//	
		//	data_model.developable_model.non_developable.V = data_model.target.V;
		//	data_model.developable_model.non_developable.F = data_model.target.F;
		//	data_model.developable_model.face_patch_assignment.clear();

		//	view_model.target_developable_view_settings.has_changed = true;
		//	optimization.run_optimiztaion();
		//}
		if (ImGui::Button("simple project to DOGs", ImVec2(-1, 0)))
		{
			meshhelper::simple_project_to_DOGs(data_model.patches, data_model.developable_model.non_developable);
			view_model.target_developable_view_settings.has_changed = true;
		}
		if (ImGui::Button("project to selected ruled developables", ImVec2(-1, 0)))
		{
			meshhelper::simple_project_to_ruled_surfaces(data_model.ruled_developables, data_model.developable_model.non_developable);
			view_model.target_developable_view_settings.has_changed = true;
		}
	}
}

void MenuUpdater::update_wrapper_mesh_menu()
{
	ImGui::Spacing();
	if (!ImGui::CollapsingHeader("DOGs##visualization"))
		return;

	if (data_model.patches.size() <= 0)
		deactivate_elements();

	ImGui::Checkbox("apply to all", &view_model.use_general_patch_settings);
	ImGui::Checkbox("show only", &view_model.show_only_current_patch);

	if (view_model.use_general_patch_settings)
		deactivate_elements();

		if (ImGui::InputInt("selected DOG", &view_model.selected_patch_index))
			view_model.selected_patch_index = clip(view_model.selected_patch_index, 0, data_model.patches.size() - 1);
	
	if (view_model.use_general_patch_settings)
		activate_elements();


	if (view_model.use_general_patch_settings)
	{
		patch_view_options(view_model.general_patch_view_settings);
	}
	else
	{
		if (view_model.selected_patch_index == -1)
			patch_view_options(PatchViewSettings());
		else
			patch_view_options(view_model.patch_view_settings[view_model.selected_patch_index]);
	}


	if (data_model.patches.size() <= 0)
		activate_elements();

	if (ImGui::Button("delete all patches"))
	{
		data_model.patches.clear();
		data_model.wrapper_view_indices.clear();

		for (auto settings : view_model.patch_view_settings)
			data_model.viewer.data_list[settings.view_index].clear();

		view_model.patch_view_settings.clear();
		view_model.selected_patch_index = -1;
	}
}

void MenuUpdater::update_optimization_menu(OptimizationController& optimization)
{
	ImGui::Spacing();
	if (!ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen))
		return;

	ImGui::Indent();

	if (ImGui::CollapsingHeader("Optimize stepwise"))
	{
		ImGui::Text("optimize stepwise:");

		if (ImGui::Button("initialize geodesics", ImVec2(-1, 0)))
		{
			optimization.geodesics_controller->initialize();
		}
		if (ImGui::Button("initialize optimization", ImVec2(-1, 0)))
		{
			optimization.initialize();
		}
		if (ImGui::Button("add patch", ImVec2(-1, 0)))
		{
			optimization.try_add_patch();
		}
		if (ImGui::Button("run 1 iteration", ImVec2(-1, 0)))
		{
			data_model.is_optimizing = true;
			optimization.run_optimiztaion();
			data_model.is_optimizing = false;
		}
		if (ImGui::Button((data_model.is_optimizing) ? "stop optimiztaion" : "run optimiztaion", ImVec2(-1, 0)))
		{
			data_model.is_optimizing = !data_model.is_optimizing;
			optimization.run_optimiztaion();
		}
	}

	if (ImGui::CollapsingHeader("Optimization settings"))
	{
		ImGui::Checkbox("pause on new patch?", &data_model.is_pausing_on_adding_patch);
		ImGui::Checkbox("pause on new constraint strategy?", &data_model.is_pausing_on_changing_constraint_strategy);
		ImGui::Checkbox("render only covering faces?", &data_model.is_rendering_patch_coverage);
		if (ImGui::Button("update patch rendering"))
		{
			view_model.do_update_all_patches = true;
			optimization.update_patch_rendering(true);
		}

		ImGui::Text("energy weights:");
		if (ImGui::InputFloat("bending", &data_model.optimization_settings.weight_bending_energy))
			write_log(3) << "set weight_bending_energy: " << data_model.optimization_settings.weight_bending_energy << endl;

		if (ImGui::InputFloat("isometry", &data_model.optimization_settings.weight_isometry_energy))
			write_log(3) << "set weight_isometry_energy: " << data_model.optimization_settings.weight_isometry_energy << endl;

		if (ImGui::InputFloat("length regularizer", &data_model.optimization_settings.weight_regularizing_energy))
			write_log(3) << "set weight_regularizing_energy: " << data_model.optimization_settings.weight_regularizing_energy << endl;

		if (ImGui::InputFloat("fitting", &data_model.optimization_settings.weight_fitting_energy))
			write_log(3) << "set weight_fitting_energy: " << data_model.optimization_settings.weight_fitting_energy << endl;


		ImGui::Separator();
		ImGui::Text("convergence settings:");

		ImGui::InputFloat("fitting threshold", &data_model.optimization_settings.convergence_fitting_threshold);
		ImGui::InputFloat("developability threshold", &data_model.optimization_settings.convergence_developablity_threshold);
	}

	ImGui::Unindent();

	if (ImGui::CollapsingHeader("Output", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Text(("average edge length: " + to_string(data_model.target.average_edge_length)).c_str());
		ImGui::InputFloat("outlier threshold", &data_model.optimization_settings.outlier_threshold);
		ImGui::InputFloat("outlier normals", &data_model.optimization_settings.outlier_normal_threshold);
		ImGui::InputFloat("coverage threshold", &data_model.optimization_settings.coverage_threshold);

		ImGui::Spacing();

		ImGui::Text(("# selected surfaces: " + std::to_string(data_model.selected_ruled_vertex_indices.size())).c_str());
		ImGui::TextWrapped(("indices: " + list_to_string(data_model.selected_ruled_vertex_indices)).c_str());
		
		ImGui::Text(("# found holes: " + std::to_string(data_model.hole_average_vertex.size())).c_str());
		ImGui::TextWrapped(("indices: " + list_to_string(data_model.hole_average_vertex)).c_str());

		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		int label_width = 80;

		ImGui::Text("-- approximation errors (absolute) --");
		ImGui::Text("(hausdorff // sum // RMS)");
		ImGui::Spacing();

		ImGui::Text("diagonal: %.4f", data_model.target_diagonal);
		ImGui::Spacing();

		ImGui::Text("torsal:"); ImGui::SameLine(label_width);
		if (data_model.initialization_coverage != NULL)
		{
			ImGui::Text("%.4f // %.4f // %.4f", data_model.initialization_coverage->hausdorff_distance(), data_model.initialization_coverage->distance_sum(), data_model.initialization_coverage->root_mean_square_error());
			ImGui::Text(""); ImGui::SameLine(label_width); ImGui::Text("range: %.4f - %.4f", data_model.initialization_coverage->distances().minCoeff(), data_model.initialization_coverage->distances().maxCoeff());
		}

		ImGui::Text("DOGs:"); ImGui::SameLine(label_width);
		if (data_model.patch_coverage != NULL)
		{
			ImGui::Text("%.4f // %.4f // %.4f", data_model.patch_coverage->hausdorff_distance(), data_model.patch_coverage->distance_sum(), data_model.patch_coverage->root_mean_square_error());
			ImGui::Text(""); ImGui::SameLine(label_width); ImGui::Text("range: %.4f - %.4f", data_model.patch_coverage->distances().minCoeff(), data_model.patch_coverage->distances().maxCoeff());
		}

		ImGui::Text("result:"); ImGui::SameLine(label_width);
		if (data_model.result_coverage != NULL)
		{
			ImGui::Text("%.4f // %.4f // %.4f", data_model.result_coverage->hausdorff_distance(), data_model.result_coverage->distance_sum(), data_model.result_coverage->root_mean_square_error());
			ImGui::Text(""); ImGui::SameLine(label_width); ImGui::Text("range: %.4f - %.4f", data_model.result_coverage->distances().minCoeff(), data_model.result_coverage->distances().maxCoeff());
		}
		

		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();


		ImGui::Text("-- result metrics --");
		ImGui::Spacing();

		ImGui::Text("K:            avg = %.4f, (%.4f - %.4f)", data_model.developable_model.gaussian_curvature.average(), data_model.developable_model.gaussian_curvature.min(), data_model.developable_model.gaussian_curvature.max());
		ImGui::Text("angle defect: %.4f", data_model.developable_model.output_angle_defect);
		ImGui::Text("stiching:     %.4f", data_model.developable_model.output_stiching);


		label_width = 120;
		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		ImGui::Text("-- patch energies --");
		ImGui::Spacing();

		ImGui::Text("bending:");		ImGui::SameLine(label_width); ImGui::Text("%f", DataModel::output_bending_energy);
		ImGui::Text("isometry:");		ImGui::SameLine(label_width); ImGui::Text("%f", DataModel::output_isometry_energy);
		ImGui::Text("length:"); 		ImGui::SameLine(label_width); ImGui::Text("%f", DataModel::output_regularizer_energy);

		bool is_fitted = DataModel::output_fitting_energy < data_model.optimization_settings.convergence_fitting_threshold;
		ImGui::Text("fitting:");		ImGui::SameLine(label_width);	ImGui::Text("%f", DataModel::output_fitting_energy);
		ImGui::SameLine(210); 			ImGui::Checkbox("", &is_fitted);

		bool is_DOG = DataModel::output_DOG_objective < data_model.optimization_settings.convergence_developablity_threshold;
		ImGui::Text("DOG error:");		ImGui::SameLine(label_width);	ImGui::Text("%f", DataModel::output_DOG_objective);
		ImGui::SameLine(210); 			ImGui::Checkbox("", &is_DOG);

		const string local_minimum_info = DataModel::is_local_minimum ? "local MINIMUM: changing weights" : " ";
		ImGui::Text(local_minimum_info.c_str());

		const string constraints_step_info = "current step: " + std::to_string(DataModel::current_constraints_step);
		ImGui::Text(constraints_step_info.c_str());
		const string convergence_info = "last converged " + std::to_string(data_model.iterations_since_converge) + " iterations ago.";
		ImGui::Text(convergence_info.c_str());
	}

}

void MenuUpdater::update_ruled_geodesics_menu(RuledGeodesicsController& ruled_geodesics)
{
	ImGui::Spacing();
	if (!ImGui::CollapsingHeader("Ruled Geodesics initialization", ImGuiTreeNodeFlags_DefaultOpen))
		return;

	ImGui::Indent();


	ImGui::Spacing();


	//bool has_trace_direction_changed = false;
	//has_trace_direction_changed = has_trace_direction_changed | ImGui::Checkbox("k_min", &view_model.do_trace_kmin);
	//ImGui::SameLine();
	//has_trace_direction_changed = has_trace_direction_changed | ImGui::Checkbox("k_max", &view_model.do_trace_kmax);
	//
	//if (has_trace_direction_changed)
	//{
	//	if (view_model.do_trace_kmin && view_model.do_trace_kmax)
	//		data_model.geodesics_direction = GeodesicsDirection::K_BOTH;
	//	else if (view_model.do_trace_kmin)
	//		data_model.geodesics_direction = GeodesicsDirection::K_MIN;
	//	else if (view_model.do_trace_kmax)
	//		data_model.geodesics_direction = GeodesicsDirection::K_MAX;
	//
	//	//data_model.do_rebuild_geodesics = true;
	//}


	
	if (ImGui::CollapsingHeader("geodesics stopping settings:"))
	{
		ImGui::Checkbox("use flat detection?", &data_model.use_geodesics_flat_detetction);
		ImGui::InputInt("window size", &data_model.geodesics_stopping.window_size);
		ImGui::InputDouble("outlier factor", &data_model.geodesics_stopping.outlier_factor);
		ImGui::InputDouble("winding factor", &data_model.geodesics_stopping.winding_factor);
		ImGui::Indent();
		ImGui::InputDouble("absolute threshold", &data_model.geodesics_stopping.absolute_threshold);
		ImGui::InputDouble("min bounds (abs)", &data_model.geodesics_stopping.bounds_min_absolute);
		ImGui::Checkbox("upper bound only?", &data_model.geodesics_stopping.use_upper_bound_only);
		ImGui::Checkbox("adaptive bounds?", &data_model.geodesics_stopping.use_adaptive_bounds_update);
		ImGui::Unindent();
	}

	if (ImGui::Button("trace all geodesics", ImVec2(-1, 0)))
	{
		ruled_geodesics.build_geodesics();

		view_model.do_update_init_geodesics = true;
		write_log(4) << "view:: tracing all geodesics" << endl;
	}

	ImGui::Text(("# traced geodesics: " + std::to_string(data_model.geodesics_candidates.paths.size())).c_str());
	ImGui::Separator;


	ImGui::Spacing();

	ImGui::InputInt("number random surfaces", &data_model.number_random_points);
	//ImGui::InputFloat("spline smoothness", &data_model.spline_smoothness);
	ImGui::InputFloat("ruled width", &data_model.ruled_width);

	if (ImGui::Button("update random surfaces", ImVec2(-1, 0)))
	{
		if (data_model.geodesics_candidates.is_empty())
			ruled_geodesics.build_geodesics();

		ruled_geodesics.build_random_developables(data_model.number_random_points);
		ruled_geodesics.select_ruled_developables();

		view_model.do_update_init_geodesics = true;
		//view_model.are_init_surfaces_visible = true;
		write_log(4) << "view:: updating " << data_model.selected_ruled_vertex_indices.size() << " random developables" << endl;
	}
	ImGui::Text(("# random vertices: " + std::to_string(data_model.ruled_vertex_indices.size())).c_str());

	//if (ImGui::Button("update random vertices", ImVec2(-1, 0)))
	//{
	//	if(data_model.geodesics_candidates.is_empty())
	//		ruled_geodesics.build_geodesics();
	//
	//	data_model.ruled_vertex_indices = ruled_geodesics.get_random_vertices(data_model.number_random_points);
	//	ruled_geodesics.build_ruled_developables();
	//
	//	view_model.do_update_init_geodesics = true;
	//	//view_model.are_init_surfaces_visible = true;
	//	write_log(4) << "view:: updating "<< data_model.number_random_points << " random vertices" << endl;
	//}
	//ImGui::Text(("# random vertices: " + std::to_string(data_model.ruled_vertex_indices.size())).c_str());
	ImGui::Separator;



	//ImGui::Spacing();
	//
	//ImGui::InputFloat("spline smoothness", &data_model.spline_smoothness);
	//ImGui::InputFloat("ruled width", &data_model.ruled_width);
	//if (ImGui::Button("update ruled developables", ImVec2(-1, 0)))
	//{
	//	if (data_model.geodesics_candidates.is_empty())
	//		ruled_geodesics.build_geodesics();
	//
	//	if(data_model.ruled_vertex_indices.size() != data_model.number_random_points)
	//		data_model.ruled_vertex_indices = ruled_geodesics.get_random_vertices(data_model.number_random_points);
	//
	//	ruled_geodesics.build_ruled_developables();
	//
	//	view_model.do_update_init_geodesics = true;
	//	//view_model.are_init_surfaces_visible = true;
	//	write_log(4) << "view:: updating " << data_model.number_random_points << " ruled developables" << endl;
	//}
	//ImGui::Text(("# random surfaces: " + std::to_string(data_model.ruled_developables.size())).c_str());
	//ImGui::Separator;



	ImGui::Spacing();

	ImGui::InputInt("label selection smoothness", &data_model.label_selection_smoothness);
	if (ImGui::Button("update selected developables", ImVec2(-1, 0)))
	{
		if (data_model.geodesics_candidates.is_empty())
			ruled_geodesics.build_geodesics();

		//if (data_model.ruled_vertex_indices.size() != data_model.number_random_points)
		//	data_model.ruled_vertex_indices = ruled_geodesics.get_random_vertices(data_model.number_random_points);

		//if(data_model.ruled_developables.size() != data_model.number_random_points)
		//	ruled_geodesics.build_ruled_developables();

		if (data_model.ruled_developables.size() < 1)
			ruled_geodesics.build_random_developables(data_model.number_random_points);

		ruled_geodesics.select_ruled_developables();

		view_model.do_update_init_geodesics = true;
		view_model.target_view_settings.has_changed = true;
		write_log(4) << "view:: updating " << data_model.number_random_points << " ruled developables" << endl;
	}
	ImGui::Text(("# selected surfaces: " + std::to_string(data_model.selected_ruled_vertex_indices.size())).c_str());
	ImGui::TextWrapped(("indices: " + list_to_string(data_model.selected_ruled_vertex_indices)).c_str());
	ImGui::Separator;


	ImGui::Spacing();

	ImGui::Text("show");

	if (ImGui::Checkbox("all", &view_model.show_all_geodesics))
	{
		if (view_model.show_all_geodesics)
		{
			view_model.show_random_geodesics = false;
			view_model.show_selected_geodesics = false;
		}
		view_model.do_update_init_geodesics = true;
	}

	ImGui::SameLine();
	if (ImGui::Checkbox("random", &view_model.show_random_geodesics))
	{
		if (view_model.show_random_geodesics)
		{
			view_model.show_all_geodesics = false;
			view_model.show_selected_geodesics = false;
		}
		view_model.do_update_init_geodesics = true;
	}

	ImGui::SameLine();
	if (ImGui::Checkbox("selected", &view_model.show_selected_geodesics))
	{
		if (view_model.show_selected_geodesics)
		{
			view_model.show_all_geodesics = false;
			view_model.show_random_geodesics = false;
		}
		view_model.do_update_init_geodesics = true;
	}

	ImGui::Indent();
	if(ImGui::Checkbox("points", &view_model.show_points_geodesics))
		view_model.do_update_init_geodesics = true;
	if (ImGui::Checkbox("geodesics", &view_model.show_curves_geodesics))
		view_model.do_update_init_geodesics = true;
	if (ImGui::Checkbox("surfaces", &view_model.show_surfaces_geodesics))
		view_model.do_update_init_geodesics = true;
	ImGui::Unindent();

	if (ImGui::InputInt("selected ruled developable", &view_model.selected_ruled_index))
	{
		view_model.selected_ruled_index = clip(view_model.selected_ruled_index, -1, data_model.selected_ruled_vertex_indices.size() - 1);
		view_model.do_update_init_geodesics = true;
		
		int print_vertex_index = view_model.selected_ruled_index >= 0 ? data_model.selected_ruled_vertex_indices[view_model.selected_ruled_index] : view_model.selected_ruled_index;
		write_log(3) << "selected ruled developable: " << print_vertex_index << endl;
	}

	if (ImGui::Button("print ruled data"))
	{
		if (view_model.selected_ruled_index >= 0 && view_model.selected_ruled_index < data_model.selected_ruled_vertex_indices.size())
		{
			int geodesic_vertex_index = data_model.selected_ruled_vertex_indices[view_model.selected_ruled_index];
			int index = index_of(geodesic_vertex_index, data_model.ruled_vertex_indices);
			
			RuledDevelopableSurface* surface = data_model.ruled_developables[index];
			surface->print_mathematica_data();
		}
	}

	ImGui::Unindent();
}

void MenuUpdater::update_geodesics_menu()
{
	ImGui::Spacing();

	if (ImGui::CollapsingHeader("Debug Geodesics"))
	{
		ImGui::Indent();

		if (ImGui::CollapsingHeader("geodesics settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Text("select walking geodesic:");

			//bool is_selecting_walking_geodesic = view_model.current_selection_mode == SelectionMode::WalkingGeodesic;
			//if (ImGui::Checkbox("select debug geodesic at vertex", &is_selecting_walking_geodesic))
			//{
			//	view_model.current_selection_mode = is_selecting_walking_geodesic ? SelectionMode::WalkingGeodesic : SelectionMode::NotSelecting;
			//}
			if (ImGui::Checkbox("select debug geodesic at vertex", &view_model.is_selecting_debug_geodesic))
			{
				view_model.current_selection_mode = view_model.is_selecting_debug_geodesic ? SelectionMode::VertexPicking : SelectionMode::NotSelecting;
			}

			if (ImGui::Checkbox("show debug geodesics", &view_model.do_render_debug_walking_geodesics))
				view_model.do_update_geodesics = true;

			//if (ImGui::Checkbox("show extended info", &view_model.show_debug_walking_geodesics_neighbors))
			//	view_model.do_update_geodesics = true;

			if (ImGui::Button("clear debug geodesics"))
			{
				view_model.debug_walking_geodesics.clear();
				view_model.do_update_geodesics = true;
			}
		}

		ImGui::Unindent();
	}
	else if (ImGui::CollapsingHeader("Debug Rulings", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ImGui::Checkbox("debug rulings at vertex", &view_model.is_debugging_ruled_surfaces))
		{
			view_model.current_selection_mode = view_model.is_debugging_ruled_surfaces ? SelectionMode::VertexPicking : SelectionMode::NotSelecting;

			if (!view_model.is_debugging_ruled_surfaces)
			{
				view_model.selected_debug_ruled_surface_index = -1;
				view_model.debug_ruled_surface = NULL;
				view_model.do_update_debug_ruled_surfaces = true;
			}
		}

		if(ImGui::RadioButton("auto", &view_model.use_debug_rulings_type, 0))
			view_model.do_update_debug_ruled_surfaces = true;
		if (ImGui::RadioButton("analytic", &view_model.use_debug_rulings_type, 1))
			view_model.do_update_debug_ruled_surfaces = true;
		if (ImGui::RadioButton("discrete", &view_model.use_debug_rulings_type, 2))
			view_model.do_update_debug_ruled_surfaces = true;
		if (ImGui::RadioButton("compound", &view_model.use_debug_rulings_type, 3))
			view_model.do_update_debug_ruled_surfaces = true;

		if (ImGui::Button("print data") && view_model.debug_ruled_surface != NULL)
			view_model.debug_ruled_surface->print_mathematica_data();

		if(ImGui::InputDouble("smooth", &GlobalSettings::rulings_spline_smoothness))
			view_model.do_update_debug_ruled_surfaces = true;
		if(ImGui::InputDouble("sampling", &GlobalSettings::rulings_spline_sample_factor)) //number of samples per avg_edge_length
			view_model.do_update_debug_ruled_surfaces = true;
		if(ImGui::InputInt("degree", &GlobalSettings::rulings_spline_degree))
			view_model.do_update_debug_ruled_surfaces = true;
	}
}

void MenuUpdater::update_preparataion_menu()
{
	/*
	ImGui::Spacing();

	if (!ImGui::CollapsingHeader("Prepare optimization"))
		return;

	update_surface_analysis_menu();
	*/
}

void MenuUpdater::update_surface_analysis_menu()
{
	ImGui::Spacing();
	if (!ImGui::CollapsingHeader("Overlay surface properties"))
		return;

	//if (ImGui::CollapsingHeader("Overlay surface properties"), ImGuiTreeNodeFlags_DefaultOpen)
	//{

		ImGui::Text("curvatures:");

		int choice = view_model.overlay_visualization;
		//int choice = view_model.current_target_visualization;

		ImGui::RadioButton("none", &choice, MeshVisualization::Standard);

		ImGui::RadioButton("vertex labelling", &choice, MeshVisualization::VertexLabelling);
		ImGui::RadioButton("ruled assignment", &choice, MeshVisualization::LabelAssignment);
		ImGui::RadioButton("error distances", &choice, MeshVisualization::ErrorDistances);

		ImGui::RadioButton("Gauss K", &choice, MeshVisualization::GaussianCurvature);
		ImGui::Text(get_curvature_range_info(data_model.target.surface_features().K).c_str());

		ImGui::RadioButton("Gauss Kp", &choice, MeshVisualization::GaussianCurvaturePrincipal);
		ImGui::Text(get_curvature_range_info(data_model.target.surface_features().Kp).c_str());

		ImGui::RadioButton("principal k", &choice, MeshVisualization::PrincipalCurvature);

		ImGui::RadioButton("mean H", &choice, MeshVisualization::MeanCurvature);
		ImGui::Text(get_curvature_range_info(data_model.target.surface_features().H).c_str());

		/*
		if (ImGui::CollapsingHeader("weighted surface properties"))
		{
			ImGui::RadioButton("weighted K + H", &choice, TargetVisualization::WeightedCurvature);
			ImGui::Text(get_curvature_range_info(data_model.target.surface_features().weighted_curvature).c_str());

			if (ImGui::SliderFloat("w_K", &view_model.w_K, 0.0, 1.0))
			{
				view_model.w_H = 1.0 - view_model.w_K;
				surface::compute_weighted_features(data_model.target, view_model.w_K, 1.0, view_model.w_H);
				view_model.target_view_settings.has_changed = true;
			}
			if (ImGui::SliderFloat("w_H", &view_model.w_H, 0.0, 1.0))
			{
				view_model.w_K = 1.0 - view_model.w_H;
				surface::compute_weighted_features(data_model.target, view_model.w_K, 1.0, view_model.w_H);
				view_model.target_view_settings.has_changed = true;
			}
		}
		*/
	
		//if (choice != view_model.current_target_visualization)
		//	view_model.target_view_settings.has_changed = true;
		//view_model.current_target_visualization = static_cast<TargetVisualization>(choice);
	
		if (choice != view_model.overlay_visualization)
			view_model.has_overlay_changed = true;
		view_model.overlay_visualization = static_cast<MeshVisualization>(choice);


		if (ImGui::CollapsingHeader("adjust visualization"), ImGuiTreeNodeFlags_DefaultOpen)
		{
			//if (ImGui::SliderFloat("filter top X%", &view_model.filter_curvature, 0.0, 1.0))
			//	view_model.target_view_settings.has_changed = true;

			if (ImGui::InputFloat("zero value", &view_model.visualized_min))
			{
				//view_model.visualized_min = view_model.visualized_max * -1;
				//view_model.target_view_settings.has_changed = true;
				view_model.target_view_settings.has_changed = true; //TODO delete!!

				view_model.overlay_settings.visualized_max = 1e6;
				view_model.overlay_settings.has_changed = true;
			}
			//if (ImGui::InputFloat("visualized value", &view_model.visualized_max))
			//{
			//	//view_model.visualized_min = view_model.visualized_max * -1;
			//	//view_model.target_view_settings.has_changed = true;
			//	view_model.target_view_settings.has_changed = true; //TODO delete!!

			//	view_model.overlay_settings.visualized_min = view_model.overlay_settings.visualized_max * -1;
			//	view_model.overlay_settings.has_changed = true;
			//}

			if (ImGui::Checkbox("filtered: show points?", &view_model.show_filtered_curvature_points))
				view_model.overlay_settings.has_changed = true;
			if (ImGui::Checkbox("filtered: show values?", &view_model.show_filtered_curvature_values))
				view_model.overlay_settings.has_changed = true;
		}

		if (ImGui::CollapsingHeader("creases"))
		{
			ImGui::Checkbox("show", &view_model.show_creases);

			ImGui::Indent();
			ImGui::Checkbox("label vertices", &view_model.label_crease_vertex);
			ImGui::Checkbox("label value", &view_model.label_crease_value);
			ImGui::Unindent();

			if (ImGui::SliderFloat("crease %", &view_model.normalized_creases_threshold, 0.0, 1.0))
			{
				double value;
				convert_distribution_percent_to_value(data_model.target.surface_features().crease_vertices, data_model.target.surface_features().crease_vertices_lookup, view_model.normalized_creases_threshold, value);

				view_model.angle_crease_threshold = value * 180.0 / igl::PI;
				view_model.overlay_settings.has_changed = true;
			}
			if (ImGui::SliderFloat("angle difference (deg)", &view_model.angle_crease_threshold, 0.0, 180.0))
			{
				const double angle_rad = view_model.angle_crease_threshold / 180.0 * igl::PI;
				double percent;
				convert_distribution_value_to_percent(data_model.target.surface_features().crease_vertices, data_model.target.surface_features().crease_vertices_lookup, angle_rad, percent);

				view_model.normalized_creases_threshold = percent;
				view_model.overlay_settings.has_changed = true;
			}
		}
	//}
}

void MenuUpdater::update_features_menu(InteractiveSegmentSelector& segment_selection)
{
	/*
	ImGui::Spacing();
	if (!ImGui::CollapsingHeader("Points of interest"))
		return;


	bool is_selecting_features = view_model.current_selection_mode == SelectionMode::DefineFeatures;
	if (ImGui::Checkbox("select features (POIs)", &is_selecting_features))
	{
		view_model.current_selection_mode = is_selecting_features ? SelectionMode::DefineFeatures : SelectionMode::NotSelecting;
		segment_selection.toggle_segmentation(is_selecting_features && view_model.current_interaction_mode == InteractionMode::Select);
	}

	ImGui::Checkbox("show 'points of interest'", &view_model.do_render_POI);

	ImGui::Text(("POI clusters: " + std::to_string(data_model.clustered_feature_vertices.size())).c_str());
	if (ImGui::InputInt("selected segment", &segment_selection.segment_index))
		segment_selection.segment_index = clip(segment_selection.segment_index, -1, data_model.clustered_feature_vertices.size() - 1);


	if (ImGui::Button("finalize crease"))
		segment_selection.finalize_segment();
	ImGui::SameLine();

	if (ImGui::Button("delete crease"))
		segment_selection.delete_segment();
	ImGui::SameLine();

	if (ImGui::Button("add boundary"))
	{
		//TODO add hashing to check if this was added already https://ideone.com/tieHbd
		vector<vector<int>> boundaries;
		igl::boundary_loop(data_model.target.F, boundaries);

		for(auto& boundary_indices : boundaries)
			segment_selection.add_segment(boundary_indices);
	}


	if (ImGui::Button("print segmentation", ImVec2(-1, 0)))
		write_log(0) << segment_selection.to_string();
	*/
}

void MenuUpdater::mesh_view_options(MeshViewSettings& settings)
{
	string id = "##" + to_string(settings.view_index);

	if (ImGui::Checkbox(("show" + id).c_str(), &settings.show_mesh))
		settings.has_changed = true;


	if (!settings.show_mesh)
		deactivate_elements();


	ImGui::SameLine();
	if(ImGui::Checkbox(("wire" + id).c_str(), &settings.show_wireframe))
		settings.has_changed = true;
	ImGui::SameLine();
	if (ImGui::Checkbox(("fill" + id).c_str(), &settings.show_faces))
		settings.has_changed = true;
	ImGui::SameLine();
	if (ImGui::Checkbox(("boundary" + id).c_str(), &settings.show_boundary))
		settings.has_changed = true;

	if (ImGui::Checkbox(("label V" + id).c_str(), &settings.label_vertices))
		settings.has_changed = true;
	ImGui::SameLine();
	if (ImGui::Checkbox(("label F" + id).c_str(), &settings.label_faces))
		settings.has_changed = true;


	if (!settings.show_mesh)
		activate_elements();
}

void MenuUpdater::patch_view_options(PatchViewSettings& settings)
{
	mesh_view_options(settings);

	if (settings.has_changed && !settings.show_mesh)
	{
		settings.label_quad_faces = false;
		settings.show_position_constraints = false;
		settings.show_target_constraints = false;
	}

	ImGui::SameLine();
	ImGui::Checkbox("label quad F", &settings.label_quad_faces);
	ImGui::Checkbox("show position constraints", &settings.show_position_constraints);
	ImGui::Checkbox("show target constraints", &settings.show_target_constraints);
}

void MenuUpdater::deactivate_elements()
{
	ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
	ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
}

void MenuUpdater::activate_elements()
{
	ImGui::PopItemFlag();
	ImGui::PopStyleVar();
}

std::string get_curvature_range_info(const Eigen::VectorXd& curvature)
{
	if (!curvature.rows())
		return "range: undef.";

	return "range: " + to_string(curvature.minCoeff()) + " - " + to_string(curvature.maxCoeff());
}
