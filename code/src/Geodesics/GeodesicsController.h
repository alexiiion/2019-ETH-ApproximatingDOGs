#pragma once

#include "DataModel.h"

class GeodesicsController
{
public:
	GeodesicsController(DataModel& data_model) : data_model(data_model) {};

	virtual bool initialize() = 0;
	virtual bool get_next_constraints(Eigen::MatrixXd& out_target_curve, Eigen::MatrixXd& out_wrapper_curve, int& vertex) = 0;


protected:
	DataModel& data_model;
	bool is_initialzed = false;

};
