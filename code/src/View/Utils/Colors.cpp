#include "Colors.h"

#include <igl/rgb_to_hsv.h>
#include <igl/hsv_to_rgb.h>

Eigen::RowVector3d Colors::darker(Eigen::RowVector3d color)
{
	return color * 0.8;
}

Eigen::RowVector3d Colors::ligher(Eigen::RowVector3d color)
{
	return color * 1.2;
}

Eigen::RowVector3d Colors::value(double value, Eigen::RowVector3d color)
{
	Eigen::RowVector3d hsv;
	igl::rgb_to_hsv(color, hsv);
	hsv(2) = value;

	Eigen::RowVector3d rgb;
	igl::hsv_to_rgb(hsv, rgb);
	return rgb;
}

Eigen::RowVector3d Colors::saturation(double saturation, Eigen::RowVector3d color)
{
	Eigen::RowVector3d hsv;
	igl::rgb_to_hsv(color, hsv);
	hsv(1) = saturation;

	Eigen::RowVector3d rgb;
	igl::hsv_to_rgb(hsv, rgb);
	return rgb;
}