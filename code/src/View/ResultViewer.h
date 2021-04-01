#pragma once

#include "DevelopableModel.h"
#include "DevelopableOptimizationSettings.h"

class ResultViewer
{
public:
	ResultViewer(Developable& developable_model, DevelopableOptimizationSettings& optimization_settings)
		: developable_model(developable_model), optimization_settings(optimization_settings)
	{};
	~ResultViewer() 
	{};

	void update_view();
	void update_menu();
	void update_debug_menu();

private:
	Developable& developable_model;
	DevelopableOptimizationSettings& optimization_settings;

	/*
	list of mesh indices
	std::vector<MeshViewSettings*> meshes;
	overlay_view_index; 
	*/
};
