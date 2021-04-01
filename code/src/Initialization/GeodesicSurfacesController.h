#pragma once

class Mesh;
class GeodesicsModel;
class RuledSurfacesModel;

class GeodesicSurfacesController
{

public:

	//GeodesicSurfacesController(Mesh& target)
	//	: target(target)
	//{ };
	GeodesicSurfacesController(Mesh& target, GeodesicsModel& geodesics, RuledSurfacesModel& ruled_surfaces)
		: target(target), geodesics(geodesics), ruled_surfaces(ruled_surfaces)
	{ };

	void build_surfaces(int n);
	void clear();
	//select developables
	//get next constraints

private:

	Mesh& target;
	GeodesicsModel& geodesics;
	RuledSurfacesModel& ruled_surfaces;

	void create_all_surfaces();
	void create_random_surfaces(int n);

	bool try_add_surface(int vertex);
};
