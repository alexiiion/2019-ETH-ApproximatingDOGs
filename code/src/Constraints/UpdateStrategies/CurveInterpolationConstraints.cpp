#include "CurveInterpolationConstraints.h"

#include <igl/procrustes.h>
#include <igl/per_vertex_normals.h>

#include "CurveInterpolation.h"
#include "PairingModes.h"
#include "Logger.h"
#include "Utils.h"

using namespace std;


void CurveInterpolationConstraints::configure()
{
	write_log(4) << endl << "    -- configuring curve interpolation constraints " << endl;

	constraints.should_remove_outliers = false;
	constraints.should_remove_normal_outliers = false;
	constraints.should_limit_face_constraints = false;
	constraints.is_keeping_pairs_constant = true;
	constraints.should_foster_stretch = false;
	constraints.pairing_target = PairingTarget::to_SURFACE;

	previous_should_converge = constraints.optimization_settings.should_converge;
	constraints.optimization_settings.should_converge = false;

	current_timestep = 0.0;
	timestep_increment = 0.01;
	//timestep_increment = 0.02;
	//timestep_increment = 0.2;

	constraints.add(patch->wrapper_curve);

	current_timestep = timestep_increment;
	DataModel::current_constraints_step = current_timestep;
}

void CurveInterpolationConstraints::update()
{
	write_log(4) << endl << "    -- update curve interpolation constraints " << endl;
	DataModel::current_constraints_step = current_timestep;

	Eigen::MatrixXd constraints_from;
	interpolate_geodesic_constraints(constraints_from);
	constraints.set_current_constraints(constraints_from);

	current_timestep += timestep_increment;
}

bool CurveInterpolationConstraints::can_update()
{
	const bool can_interpolate_curve = current_timestep < 1.0 + timestep_increment;
	write_log(4) << endl << "    -- can update curve interpolation constraints? " << boolalpha << can_interpolate_curve << endl;

	if (!can_interpolate_curve)
	{
		constraints.optimization_settings.should_converge = previous_should_converge;
		//TODO check if should invert normals because orientation is wrong, outlier removal based on normals causes errors
	}

	return can_interpolate_curve;
}

void CurveInterpolationConstraints::interpolate_geodesic_constraints(Eigen::MatrixXd& out_constraints_from)
{
	write_log(5) << endl << endl << "    - interpolate geodesic GEODESIC" << endl << endl;

	if (current_timestep >= 1.0) //adding target curve as last interpolation step because curve inerpolation is not robust
	{
		out_constraints_from = patch->target_curve;
		//return;
	}

	CurveInterpolation wrapper_curve(patch->wrapper_curve);
	CurveInterpolation target_curve(patch->target_curve);
	
	if (wrapper_curve.is_same_as(target_curve))
	{
		current_timestep = 1.0;
		write_log(3) << "    curve interpolation:: curves are THE SAME, fast forwarding to timestep 1.0" << endl;
	}

	CurveInterpolation interpolated_curve(wrapper_curve, target_curve, current_timestep); // At 0 its wrapperCurve, at 1 targetCurve
	write_log(3) << "       curve interpolation:: current_timestep = " << current_timestep << endl << endl;

	//target_curve.print_geometric_represenation();

	Eigen::RowVector3d translation = Eigen::RowVector3d::Zero();
	Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
	Eigen::MatrixXd interpolated_curve_coords = interpolated_curve.getCoords(translation, rotation);

	if (LOG_LEVEL >= 6)
	{
		cout << "target geodesic: " << endl << patch->target_curve << endl << endl;
		cout << "wrapper geodesic: " << endl << patch->wrapper_curve << endl << endl;
		cout << "interpolated_curve_coords: " << endl << interpolated_curve_coords << endl << endl;
	
		std::cout << "target: " << std::endl; target_curve.print_geometric_represenation(); std::cout << std::endl;
		std::cout << "interpolated_curve: " << std::endl; interpolated_curve.print_geometric_represenation(); std::cout << std::endl;
	}

	igl::procrustes(interpolated_curve_coords, patch->target_curve, false, false, rotation, translation);
	out_constraints_from = (interpolated_curve_coords * rotation).rowwise() + translation;

	write_log(5) << "curve interpolation:: translation: " << endl << translation << endl << endl;
	write_log(5) << "curve interpolation:: rotation: " << endl << rotation << endl << endl;

	write_log(6) << "curve interpolation:: points from: " << endl << constraints.paired_points_from << endl << endl;
	write_log(6) << "curve interpolation:: points to: " << endl << constraints.paired_points_to << endl << endl;
}
