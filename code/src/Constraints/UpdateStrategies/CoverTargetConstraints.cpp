#include "CoverTargetConstraints.h"

#include "PairingModes.h"
#include "Logger.h"
#include "Utils.h"

#include "MeshController.h"

using namespace std;

void CoverTargetConstraints::configure()
{
	write_log(3) << endl << "    -- configuring cover target constraints " << endl;

	do_update = true;

	//constraints.set_pairing_direction(PairingDirection::DOG_to_TARGET);
	//constraints.add(*constraints.pair_from_V);

	constraints.should_limit_face_constraints = true;
	constraints.should_remove_outliers = true;
	constraints.should_remove_normal_outliers = true;
	constraints.is_keeping_pairs_constant = false;
	constraints.use_previous_constraints = false;
	constraints.should_foster_stretch = true;


	previous_weights = constraints.optimization_settings;
	constraints.optimization_settings.weight_isometry_energy = 0.0;
	constraints.optimization_settings.weight_regularizing_energy = 2.0;
	//constraints.optimization_settings.weight_fitting_energy = 0.5;
	//constraints.optimization_settings.weight_bending_energy *= 2;
	//constraints.optimization_settings.outlier_threshold *= 2;

	

	constraints.pairing_target = PairingTarget::to_SURFACE;
	constraints.pairing_direction = PairingDirection::TARGET_to_DOG;
	//update pairing direction?


	//constraints.add(patch->target_curve);
	//constraints.add(*constraints.pair_from_V);

	Eigen::MatrixXd c = *constraints.pair_from_V;
	append_matrix(patch->target_curve, c);
	constraints.add(c);

	//stopping criteria setup
	current_iteration = 0;
	max_iterations = 20;

	previous_error = 0;
	previous_derivative = 0;
	min_derivative = 0.1;
	//min_derivative = 0.8;

	minimum_iterations = 0;
	//minimum_window = 100000;
	minimum_window = 3;
 

	//what is this??
	patch->target_curve.resize(0, 3);

	DataModel::current_constraints_step = 0.0;
}

void CoverTargetConstraints::update()
{
	write_log(3) << endl << "    -- update cover target constraints " << endl;
	write_log(3)         << "       update:: current_iteration = " << current_iteration << endl << endl;

	if (constraints.paired_points_from.rows() < 1)
	{
		do_update = false;
		constraints.optimization_settings = previous_weights;
		return;
	}

	//TODO measure covered area to see if it is still streching? or normalize the error by the area

	double current_error = current_error_squared();

	double derivative = previous_error - current_error;
	double derivative2 = previous_derivative - derivative;
	write_log(0) << "         iteration: " << current_iteration << ", previous: " << previous_error << ", current: " << current_error << ", der1: " << derivative << " (der2: " << derivative2 << ")" << endl;
	previous_error = current_error;
	previous_derivative = derivative;

	if (std::abs(derivative) < min_derivative)
		minimum_iterations++;
	else
		minimum_iterations = 0;

	if (minimum_iterations >= minimum_window || current_iteration >= max_iterations - 1)
		do_update = false;
	 

	Eigen::MatrixXd c = *constraints.pair_from_V;
	append_matrix(patch->target_curve, c);
	constraints.set_current_constraints(c);

	if (constraints.paired_points_from.rows() < 1)
		do_update = false;

	//constraints.set_current_constraints(*constraints.pair_from_V);
	current_iteration++;
	DataModel::current_constraints_step++;
}

bool CoverTargetConstraints::can_update()
{
	write_log(4) << endl << "    -- can update cover target constraints? " << boolalpha << do_update << endl;

	if (!do_update)
		constraints.optimization_settings = previous_weights;

	return do_update;
}

double CoverTargetConstraints::current_error_squared()
{
	double current_constraints_error_squared = 0.0;
	for (int i = 0; i < constraints.paired_points_from.rows(); i++)
		current_constraints_error_squared += (constraints.paired_points_to.row(i) - constraints.paired_points_from.row(i)).squaredNorm();

	return current_constraints_error_squared;
}
