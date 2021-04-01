#include "Serializer.h"

#include <igl/serialize.h>
#include <igl/xml/serialize_xml.h>

#include "DataModel.h"
#include "ViewModel.h"
#include "PatchModel.h"
#include "MeshController.h"

#include <functional>

//bool(*serialize_function_pointer)(const double& obj, const std::string& filename) = &igl::serialize<double>;
//bool(*serialize_function_pointer)(const std::string& obj, const std::string& filename) = &igl::serialize<std::string>;
//bool(*serialize_function_pointer)(const double& obj, const std::string& filename);
//bool(*serialize_function_pointer)(const std::string& obj, const std::string& filename);

//template <typename T>
//void persist(T& obj, std::string& obj_name, const std::string& filename, std::function<bool(const T&, const std::string&, const std::string&)> serialize_fcn)
//{
//	serialize_fcn(obj, obj_name, filename);
//}



//void serialize_feature_geodesics(const std::string& filename, const FeatureGeodesics& geodesics_candidates, const std::string& postfix);
//void deserialize_feature_geodesics(const std::string& filename, FeatureGeodesics& geodesics_candidates, const std::string& postfix);

void serialize_geodesics(const std::string& filename, const GeodesicCandidates& geodesics_candidates, const std::string& postfix);
void deserialize_geodesics(const std::string& filename, GeodesicCandidates& geodesics_candidates, const std::string& postfix);

void serialize_patch(const std::string& filename, const Patch& patch, const std::string postfix);
void deserialize_patch(const std::string& filename, Patch*& patch, const std::string postfix, DataModel& data_model);

void serialize(const std::string& filename, const DataModel& data_model)
{
	std::string file = data_model.models_folder + filename;

	//use overwrite = true for the first serialization to create or overwrite an existing file
	//igl::xml::serialize_xml(data_model.target_filename, "target_filename", file, false, true);
	igl::serialize(data_model.target_filename, "target_filename", file, true);
	igl::serialize(data_model.target_max_dimension, "target_max_dimension", file);
	//igl::serialize(data_model.target.scale, "target_scale", file);

	////save surface features
	//if(data_model.do_persist_features)
	//	igl::serialize(data_model.clustered_feature_vertices, "clustered_feature_vertices", file);
	//if (data_model.do_persist_graph)
	//{
	//	igl::serialize(data_model.number_connections, "number_connections", file);
	//	igl::serialize(data_model.feature_connectivity_graph, "feature_connectivity_graph", file);
	//}


	//save geodesics
	serialize_geodesics(file, data_model.geodesics_candidates, "");
	igl::serialize(data_model.geodesics_directions, "geodesics_directions", file);


	//save random ruled developables
	//igl::serialize(data_model.use_weighted_frame_error, "use_weighted_frame_error", file);
	//igl::serialize(data_model.max_frame_error, "max_frame_error", file);
	igl::serialize(data_model.ruled_width, "ruled_width", file);
	igl::serialize(data_model.label_selection_smoothness, "label_selection_smoothness", file);

	int number_random_ruled_developables = data_model.ruled_vertex_indices.size();
	igl::serialize(number_random_ruled_developables, "number_random_ruled_developables", file);
	igl::serialize(data_model.ruled_vertex_indices, "ruled_vertex_indices", file);


	//save selected ruled developables
	int number_selected_ruled_developables = data_model.selected_ruled_vertex_indices.size();
	igl::serialize(number_selected_ruled_developables, "number_selected_ruled_developables", file);
	igl::serialize(data_model.selected_ruled_vertex_indices, "selected_ruled_vertex_indices", file);

	
	//save DOGs
	int number_patches = data_model.patches.size();
	igl::serialize(number_patches, "number_patches", file);

	for (int i = 0; i < number_patches; i++)
		serialize_patch(file, *(data_model.patches[i]), std::to_string(i));
}


void deserialize(const std::string& filename, DataModel& data_model, ViewModel& view_model, bool is_absolute_path)
{
	std::string file = is_absolute_path ? filename : data_model.models_folder + filename;

	igl::deserialize(data_model.target_filename, "target_filename", file);

	////load features
	//if (data_model.do_persist_features)
	//	igl::deserialize(data_model.clustered_feature_vertices, "clustered_feature_vertices", file);
	////load graph
	//if (data_model.do_persist_graph)
	//{
	//	igl::deserialize(data_model.number_connections, "number_connections", file);
	//	igl::deserialize(data_model.feature_connectivity_graph, "feature_connectivity_graph", file);
	//}


	//load geodesics
	deserialize_geodesics(file, data_model.geodesics_candidates, "");
	igl::deserialize(data_model.geodesics_directions, "geodesics_directions", file);


	//load random ruled developables
	//igl::deserialize(data_model.use_weighted_frame_error, "use_weighted_frame_error", file);
	//igl::deserialize(data_model.max_frame_error, "max_frame_error", file);
	igl::deserialize(data_model.ruled_width, "ruled_width", file);
	igl::deserialize(data_model.label_selection_smoothness, "label_selection_smoothness", file);

	igl::deserialize(data_model.number_random_points, "number_random_ruled_developables", file);
	igl::deserialize(data_model.ruled_vertex_indices, "ruled_vertex_indices", file);


	//load selected ruled developables
	int number_selected_ruled_developables;
	igl::deserialize(number_selected_ruled_developables, "number_selected_ruled_developables", file);
	igl::deserialize(data_model.selected_ruled_vertex_indices, "selected_ruled_vertex_indices", file);


	//load target
	//igl::deserialize(data_model.target.scale, "target_scale", file);
	igl::deserialize(data_model.target_max_dimension, "target_max_dimension", file);
	meshhelper::add_target(data_model.target_filename, data_model, view_model.target_view_settings.view_index, view_model.target_developable_view_settings.view_index);

	//load patches
	int number_patches;
	igl::deserialize(number_patches, "number_patches", file);

	for (int i = 0; i < number_patches; i++)
	{
		Patch* patch;
		deserialize_patch(file, patch, std::to_string(i), data_model);
		data_model.patches.push_back(patch);

		meshhelper::add_mesh(data_model.viewer, patch->wrapper.V, patch->wrapper.F, false, false);

		//TODO manage indices and colors
		data_model.viewer.data().set_colors(Eigen::RowVector3d(0.1, 0.9, 0.1));
		data_model.wrapper_view_indices.push_back(data_model.viewer.selected_data_index);

		data_model.viewer.data().point_size = 10.0;
		data_model.viewer.data().line_width = 0.5;
	}

	//concatenate patches
	meshhelper::concatenate_wrappers(data_model, data_model.concatenated_patches_V, data_model.concatenated_patches_F);
}

void serialize_geodesics(const std::string& filename, const GeodesicCandidates& geodesics_candidates, const std::string& postfix)
{
	igl::serialize(geodesics_candidates.paths, "geodesics_candidates" + postfix + ".paths", filename);
	igl::serialize(geodesics_candidates.lengths, "geodesics_candidates" + postfix + ".lengths", filename);
	igl::serialize(geodesics_candidates.lengths_lookup, "geodesics_candidates" + postfix + ".lengths_lookup", filename);
	igl::serialize(geodesics_candidates.number_geodesics, "geodesics_candidates" + postfix + ".number_geodesics", filename);
}

void deserialize_geodesics(const std::string& filename, GeodesicCandidates& geodesics_candidates, const std::string& postfix)
{
	igl::deserialize(geodesics_candidates.paths, "geodesics_candidates" + postfix + ".paths", filename);
	igl::deserialize(geodesics_candidates.lengths, "geodesics_candidates" + postfix + ".lengths", filename);
	igl::deserialize(geodesics_candidates.lengths_lookup, "geodesics_candidates" + postfix + ".lengths_lookup", filename);
	igl::deserialize(geodesics_candidates.number_geodesics, "geodesics_candidates" + postfix + ".number_geodesics", filename);
}

void serialize_patch(const std::string& filename, const Patch& patch, const std::string postfix)
{
	igl::serialize(patch.target_curve , "patch" + postfix + ".target_curve", filename);
	igl::serialize(patch.wrapper_curve , "patch" + postfix + ".wrapper_curve", filename);
	igl::serialize(patch.geodesic_vertex , "patch" + postfix + ".geodesic_vertex", filename);
	//igl::serialize(patch.current_strategy , "patch" + postfix + ".current_strategy", filename);

	igl::serialize(patch.wrapper.V , "patch" + postfix + ".wrapper.V", filename);
	igl::serialize(patch.wrapper.F , "patch" + postfix + ".wrapper.F", filename);
	igl::serialize(patch.wrapper.F_quad , "patch" + postfix + ".wrapper.F_quad", filename);
	igl::serialize(patch.wrapper.initial_x0 , "patch" + postfix + ".wrapper.initial_x0", filename);
	igl::serialize(patch.wrapper.quad_width , "patch" + postfix + ".wrapper.quad_width", filename);
	igl::serialize(patch.wrapper.quad_height , "patch" + postfix + ".wrapper.quad_height", filename);
	igl::serialize(patch.wrapper.quad_topology , "patch" + postfix + ".wrapper.quad_topology", filename);
	
	igl::serialize(patch.DOG_error , "patch" + postfix + ".DOG_error", filename);
	igl::serialize(patch.fitting_error , "patch" + postfix + ".fitting_error", filename);
	igl::serialize(patch.is_local_minimum , "patch" + postfix + ".is_local_minimum", filename);
}

void deserialize_patch(const std::string& filename, Patch*& patch, const std::string postfix, DataModel& data_model)
{
	Eigen::MatrixXd target_curve, wrapper_curve;
	int geodesic_vertex;
	igl::deserialize(target_curve, "patch" + postfix + ".target_curve", filename);
	igl::deserialize(wrapper_curve, "patch" + postfix + ".wrapper_curve", filename);
	igl::deserialize(geodesic_vertex, "patch" + postfix + ".geodesic_vertex", filename);

	WrapperMesh wrapper;
	igl::deserialize(wrapper.V, "patch" + postfix + ".wrapper.V", filename);
	igl::deserialize(wrapper.F, "patch" + postfix + ".wrapper.F", filename);
	igl::deserialize(wrapper.F_quad, "patch" + postfix + ".wrapper.F_quad", filename);
	igl::deserialize(wrapper.initial_x0, "patch" + postfix + ".wrapper.initial_x0", filename);
	igl::deserialize(wrapper.quad_width, "patch" + postfix + ".wrapper.quad_width", filename);
	igl::deserialize(wrapper.quad_height, "patch" + postfix + ".wrapper.quad_height", filename);
	igl::deserialize(wrapper.quad_topology, "patch" + postfix + ".wrapper.quad_topology", filename);

	PositionConstraints constraints(data_model.target, wrapper, data_model.optimization_settings);
	patch = new Patch(data_model, target_curve, wrapper_curve, wrapper, geodesic_vertex, constraints);

	igl::deserialize(patch->DOG_error, "patch" + postfix + ".DOG_error", filename);
	igl::deserialize(patch->fitting_error, "patch" + postfix + ".fitting_error", filename);
	igl::deserialize(patch->is_local_minimum, "patch" + postfix + ".is_local_minimum", filename);
}
