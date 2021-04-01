#pragma once

#include "DataModel.h"
#include "PositionConstraintsModel.h"

class ConstraintsStrategy;

class Patch
{

public:
	Patch(DataModel& data_model, Eigen::MatrixXd& target_curve, Eigen::MatrixXd& wrapper_curve, WrapperMesh& wrapper, int& vertex, PositionConstraints& constraints)
		: data_model(data_model), target_curve(target_curve), wrapper_curve(wrapper_curve), geodesic_vertex(vertex), wrapper(wrapper), constraints(constraints)
	{};

	Patch(DataModel& data_model, Eigen::MatrixXd& target_curve, Eigen::MatrixXd& wrapper_curve, int& vertex);
	Patch(DataModel& data_model, Eigen::MatrixXd& target_curve, int& vertex);

	~Patch();

	//if multiple curves for _one_ patch (like orthogonal geodesics), then use list here!
	Eigen::MatrixXd target_curve;
	Eigen::MatrixXd wrapper_curve;
	int geodesic_vertex;

	WrapperMesh wrapper;
	PositionConstraints constraints;
	
	std::vector<ConstraintsStrategy*> constraints_strategies;
	int current_strategy;

	double fitting_error;
	double DOG_error;
	bool is_local_minimum;
	bool has_strategy_changed;

	void update();
	//bool should_renew_constraints();

	bool is_done_approximating();
	bool has_converged();
	bool is_last_strategy();

	ConstraintsStrategy* get_current_strategy()
	{
		return constraints_strategies[current_strategy];
	}


private:
	DataModel& data_model;

	void add_wrapper();
	void create_wrapper();
	void align_wrapper();
	void resample_constraint_curves(double segment_length);


	void renew_position_constraints();
	void advance_constraints_strategy();

};