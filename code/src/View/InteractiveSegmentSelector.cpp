#include "InteractiveSegmentSelector.h"

#include <igl/dijkstra.h>

#include "Utils.h"
#include "Logger.h"

using namespace std;


void InteractiveSegmentSelector::toggle_segmentation(bool is_selecting)
{
	write_log(5) << "toggle feature selection: " <<
		"is selecting: " << boolalpha << is_selecting <<
		",number of features: " << selected_vertices.size() <<
		", current index: " << segment_index;
	if (selected_vertices.size() > 0)
		write_log(5) << ", feature size: " << selected_vertices.back().size() << endl;
	else
		write_log(5) << endl;


	if (selected_vertices.size() < 1)
	{
		selected_vertices.push_back(vector<int>());
		segment_index = 0;
		return;
	}

	if (is_selecting && segment_index < 0 && selected_vertices.size() > 0)
		segment_index = selected_vertices.size() - 1;

	std::vector<int> *current_loop = &selected_vertices[segment_index];

	if (!is_selecting && current_loop->size() < 1)
	{
		selected_vertices.pop_back();
		segment_index--;
		return;
	}

	if (is_selecting && current_loop->size() > 1 && current_loop->back() == end_signifier)
	{
		selected_vertices.push_back(vector<int>());
		segment_index++;
		return;
	}


	if (is_selecting && current_loop->size() < 1)
		return;

	if (is_selecting && current_loop->front() == current_loop->back()) //if current loop is closed
		finalize_segment();
}

void InteractiveSegmentSelector::finalize_segment()
{
	std::vector<int> *current_loop = &selected_vertices[segment_index];

	if (current_loop->size() < 1)
		return;

	if (current_loop->front() == current_loop->back() && current_loop->size() > 1) //if current loop is closed
		current_loop->pop_back();

	if (current_loop->back() != end_signifier)
		current_loop->push_back(end_signifier); //signifies that this crease is final

	selected_vertices.push_back(vector<int>());
	segment_index++;
}

bool InteractiveSegmentSelector::add_path()
{
	if (data_model.pre_selected_vertex < 0)
		return false;

	std::vector<int> *current_loop = &selected_vertices[segment_index];
	int loop_size = current_loop->size();

	if (loop_size < 1) //if segment is empty, add selected vertex
	{
		current_loop->push_back(data_model.pre_selected_vertex);

		write_log(0) << "adding first index to loop" << endl;
		log_list(0, *current_loop, "current_loop: ", false);

		return false;
	}

	if (pre_segment_path.size() <= 0)
		return false;

	reverse(pre_segment_path.begin(), pre_segment_path.end());
	current_loop->insert(current_loop->end(), pre_segment_path.begin() + 1, pre_segment_path.end());
	pre_segment_path.clear();

	write_log(0) << endl << "ADD regularly" << endl;

	int first = current_loop->front();
	int last = current_loop->back();

	if (first == last) //is closing loop
	{
		log_list(0, *current_loop, "current_loop: ", false);
		write_log(0) << endl << "CLOSED segmentation loop!" << endl << endl;

		finalize_segment();

		return true;
	}

	log_list(0, *current_loop, "current_loop: ", false);
	write_log(0) << endl;

	return true;
}

void InteractiveSegmentSelector::add_segment()
{
	selected_vertices.push_back(vector<int>());
	segment_index++;
}

void InteractiveSegmentSelector::add_segment(std::vector<int>& vertex_indices)
{
	selected_vertices.push_back(vertex_indices);
	segment_index = selected_vertices.size() - 1;
	finalize_segment();
}

void InteractiveSegmentSelector::delete_segment()
{
	selected_vertices.erase(selected_vertices.begin() + segment_index);
	segment_index--;
}

void InteractiveSegmentSelector::find_edge_path()
{
	//if (data_model.pre_selected_vertex < 0)
	//	return;

	//int loop_index = selected_vertices.size() - 1;
	//if (loop_index < 0)
	//	return;

	//std::vector<int> *current_loop = &selected_vertices[loop_index];
	//int loop_size = current_loop->size();

	//if (loop_size < 1)
	//	return;

	//int source = current_loop->at(loop_size - 1);
	//write_log(6) << "segmentation: find path from " << source << " to " << data_model.pre_selected_vertex;

	//std::set<int> targets{ data_model.pre_selected_vertex };
	//Eigen::VectorXd min_distance;
	//Eigen::VectorXi previous;
	//igl::dijkstra_compute_paths(source, targets, data_model.target.adjacency_VV, min_distance, previous);
	//igl::dijkstra_get_shortest_path_to(data_model.pre_selected_vertex, previous, pre_segment_path);

	//log_list(6, pre_segment_path, "  --> dijkstra path: ", false);
}

Eigen::MatrixXd InteractiveSegmentSelector::get_segement_points_at(int index)
{
	if (index < 0 || index > selected_vertices.size() - 1)
		return Eigen::MatrixXd();

	auto& loop = selected_vertices[index];
	if (loop.size() < 1)
		return Eigen::MatrixXd();

	int size = loop.back() == end_signifier ? loop.size() - 1 : loop.size();

	Eigen::MatrixXd segment_points(size, 3);
	for (int j = 0; j < size; j++)
		segment_points.row(j) = data_model.target.V.row(loop[j]);

	return segment_points;
}

Eigen::MatrixXd InteractiveSegmentSelector::get_preview_segement_points()
{
	Eigen::MatrixXd preview_segement_points;
	index_to_value(pre_segment_path, data_model.target.V, preview_segement_points);

	return preview_segement_points;
}

 void InteractiveSegmentSelector::get_preview_segement_points(Eigen::MatrixXd& out_preview_segement_points)
{
	 index_to_value(pre_segment_path, data_model.target.V, out_preview_segement_points);
}

 std::string InteractiveSegmentSelector::to_string()
 {
	 stringstream stream;
	 stream << endl << "segments on target: " << endl;

	 for (int i = 0; i < selected_vertices.size(); i++)
		 log_list(0, selected_vertices[i], "selected_vertices: ", false);

	 stream << endl;
	 return stream.str();
 }
