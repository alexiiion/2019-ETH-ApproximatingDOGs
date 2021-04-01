#pragma once

#include "GeodesicsController.h"
#include "MeritSelector.h"

class RuledGeodesicsController : public GeodesicsController
{
public:
	RuledGeodesicsController(DataModel& data_model) : GeodesicsController(data_model) {};
	
	virtual bool initialize() override;
	virtual bool get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex) override;

	//all part of initialization, only public to re-run steps interactively from UI
	void build_geodesics();
	void build_random_developables(int n);
	void select_ruled_developables();

	//bool try_get_surface(int& vertex, RuledDevelopableSurface*& out_surface);
	bool try_get_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices);


private:
	std::vector<Eigen::VectorXd> patch_vertex_distances;
	std::vector<Eigen::VectorXd> patch_vertex_normal_deviations;
	MeritSelector* selector;

	std::vector<int> flat_indices;
	std::vector<std::vector<int>> flat_area_indices;


	bool try_get_flat_surface_at(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices);
	bool try_get_curved_surface_at(const int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& assigned_indices);

	//bool try_get_ruled_developable(const int vertex, RuledDevelopableSurface*& surface);
	//bool try_get_flat_surface(int& vertex, RuledDevelopableSurface*& out_surface, std::vector<int>& out_flat_vertices);
	std::vector<RuledDevelopableSurface*> build_developables_at(std::vector<int>& ruled_indices);
	void create_ruled_constraints_selector();

	bool select_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex);
	bool select_geodesic(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex);
	bool select_hole(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex);

	bool have_checked_holes = false;
	void check_coverage_holes();
	std::vector<std::vector<int>> label_vertices_uncovered_areas();
	std::vector<Mesh> build_uncovered_area_meshes(const std::vector<std::vector<int>>& labeled_component_vertices);


	////TODO delete
	//std::vector<int> get_random_vertices(int n);
	//void build_ruled_developables();

	//get is bad naming if it doesn't return anything
	//void get_ruled_developable(const Eigen::MatrixXd& geodesic);
	//void get_random_ruled_developables(int n);
	//void get_ruled_developables(std::vector<int> vertices);
};
