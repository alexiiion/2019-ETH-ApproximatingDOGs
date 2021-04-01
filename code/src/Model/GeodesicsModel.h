#pragma once

#include <vector>

#include "MeshModel.h"
#include "Geodesic.h"
#include "GeodesicWalker.h"

#include "DataModel.h"

class GeodesicsModel
{
public:
	//GeodesicsModel(Mesh& target)
	//	: target(target), walker(GeodesicWalker(target))
	//{
	//	initialize();
	//};
	GeodesicsModel(DataModel& data_model)
		: target(data_model.target), walker(GeodesicWalker(data_model.target, data_model.geodesics_stopping))
	{
		initialize();
	};

	//tracing settings
	// - search directions
	// - stopping criteria


	void clear();
	void initialize();
	Geodesic* at(int vertex);
	//std::vector<Geodesic*> get_all();

	int number_geodesics();

	/* Members */
	int current_candidates_index = -1;


private:

	Mesh& target;
	std::vector<Geodesic*> geodesics;
	GeodesicWalker walker;

	void trace_at(int vertex);
};