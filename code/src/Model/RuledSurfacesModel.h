#pragma once

#include "MeshModel.h"
#include "GeodesicsModel.h"
#include "Rulings/RuledDevelopableSurface.h"


struct RuledSurfaceSettings
{
	float default_width = 15;
	int number_ruled_surfaces = 100;
	int label_selection_smoothness = 50;
	double area_scale = 1.5;

	//float spline_smoothness = 0.05;
};

class RuledSurfacesModel
{
public:
	RuledSurfacesModel(Mesh& target, GeodesicsModel& geodesics) 
		: target(target), geodesics(geodesics)
	{ };


	/* Methods */

	//void initialize();
	void clear();
	bool try_get_surface_at(int& vertex, RuledDevelopableSurface*& out_surface);


	/* Members */

	RuledSurfaceSettings settings;

	//all vertics for which we randomly created surfaces
	std::vector<int> vertex_indices; 
	std::vector<int> flat_vertex_indices; //only the flat ones, they are also contained above in *vertex_indices*

	//selected ruled surfaces
	std::vector<int> selected_vertex_indices;
	std::vector<RuledDevelopableSurface*> surfaces; //sparse vector of size #V, only has entries at ruled_vertex_indices 

	//assigned label for each vertex
	std::vector<int> vertex_ruled_labelling; // #target.V


private:

	Mesh& target;
	GeodesicsModel& geodesics;

	bool try_get_flat_surface_at(int& vertex, RuledDevelopableSurface*& out_surface);
	bool try_get_curved_surface_at(const int& vertex, RuledDevelopableSurface*& out_surface);
};