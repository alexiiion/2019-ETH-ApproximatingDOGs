#pragma once

struct GlobalSettings
{

	static double flat_curvature_threshold; // the curvature up to which we consider something to be flat (i.e. 1 / osculating circle radius; should be dependent on the mesh size)
	static double flat_angle_threshold; // the angle up to which we consider something to be flat
	static double crease_angle_threshold; // the angle from when we consider something to be a crease

	static double max_average_edge;



	//Ruling settings (MOVE to own settings)
	static double rulings_spline_smoothness;
	static double rulings_spline_sample_factor; //number of samples per avg_edge_length
	static int rulings_spline_degree;

	static double min_tangent_ruling_deviation;
	static double min_ruling_length;

};
