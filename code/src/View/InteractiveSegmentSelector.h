#pragma once

#include "DataModel.h"

class InteractiveSegmentSelector
{
public:
	InteractiveSegmentSelector(DataModel& data_model) : data_model(data_model) {};

	int segment_index = -1;
	std::vector<std::vector<int>> selected_vertices;

	void toggle_segmentation(bool is_selecting);
	void finalize_segment();
	
	bool add_path();

	void find_edge_path();
	void add_segment();
	void add_segment(std::vector<int>& vertex_indices);
	void delete_segment();

	Eigen::MatrixXd get_segement_points_at(int index);
	Eigen::MatrixXd get_preview_segement_points();
	void get_preview_segement_points(Eigen::MatrixXd& out_preview_segement_points);

	std::string to_string();


private:
	DataModel& data_model;
	std::vector<int> pre_segment_path;
	const int end_signifier = -1;

};