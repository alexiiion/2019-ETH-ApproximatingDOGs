#pragma once

#include "DataModel.h"
#include "View.h"

class ApplicationController
{
public:
	ApplicationController(DataModel& data_model, View& view) : data_model(data_model), view(view) {};

private:
	DataModel& data_model;
	View& view;

	//void load_mesh();

};