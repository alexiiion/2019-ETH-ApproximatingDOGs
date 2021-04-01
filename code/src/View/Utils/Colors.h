#pragma once

#include <Eigen/Core>


namespace Colors
{
	const Eigen::RowVector3d BLACK = Eigen::RowVector3d(0.001, 0.001, 0.001);
	const Eigen::RowVector3d WHITE = Eigen::RowVector3d(0.999, 0.999, 0.999);

	const Eigen::RowVector3d GRAY_DARK			= Eigen::RowVector3d(0.3, 0.3, 0.3);
	const Eigen::RowVector3d GRAY_MID			= Eigen::RowVector3d(0.5, 0.5, 0.5);
	const Eigen::RowVector3d GRAY_LIGHT			= Eigen::RowVector3d(0.7, 0.7, 0.7);
	const Eigen::RowVector3d GRAY_ULTRALIGHT	= Eigen::RowVector3d(0.85, 0.85, 0.85);

	const Eigen::RowVector3d RED   = Eigen::RowVector3d(0.999, 0.001, 0.001);
	const Eigen::RowVector3d GREEN = Eigen::RowVector3d(0.001, 0.999, 0.001);
	const Eigen::RowVector3d BLUE  = Eigen::RowVector3d(0.001, 0.001, 0.999);
	
	const Eigen::RowVector3d CYAN     = Eigen::RowVector3d(0.001, 0.999, 0.999);
	const Eigen::RowVector3d YELLOW   = Eigen::RowVector3d(0.999, 0.999, 0.001);
	const Eigen::RowVector3d MAGENTA  = Eigen::RowVector3d(0.999, 0.001, 0.999);


	//#E6E2DF
	const Eigen::RowVector3d TAUPE_LIGHT  = Eigen::RowVector3d(230.0/255.0, 226.0/255.0, 223.0/255.0);

	//the office blue: #00B0F0
	const Eigen::RowVector3d ACCENT			= Eigen::RowVector3d(000.0/255.0, 176.0/255.0, 240.0/255.0);
	const Eigen::RowVector3d ACCENT_LIGHT	= Eigen::RowVector3d(207.0/255.0, 233.0/255.0, 243.0/255.0);
	const Eigen::RowVector3d ACCENT_DARK	= Eigen::RowVector3d(000.0/255.0, 094.0/255.0, 128.0/255.0);
	const Eigen::RowVector3d ACCENT_GRAY	= Eigen::RowVector3d(170.0/255.0, 221.0/255.0, 243.0/255.0);


	Eigen::RowVector3d darker(Eigen::RowVector3d color);
	Eigen::RowVector3d ligher(Eigen::RowVector3d color);

	Eigen::RowVector3d value(double value, Eigen::RowVector3d color);
	Eigen::RowVector3d saturation(double value, Eigen::RowVector3d color);
};


////Michaels colors from previous papers
//Eigen::Vector4d diffuse; diffuse << 135. / 255, 206. / 255, 250. / 255, D.alpha;
//Eigen::Vector4d ambient; /*ambient = 0.05*diffuse;*/ ambient << 0.05, 0.05, 0.05, 1.;
//Eigen::Vector4d specular; specular << 0, 0, 0, 1.;// specular << 0.1,0.1,0.1,1.;
//viewer.data.uniform_colors(ambient, diffuse, specular);

//Eigen::Vector3d diffuse; diffuse << 135. / 255, 206. / 255, 250. / 255;
//Eigen::Vector3d ambient; /*ambient = 0.05*diffuse;*/ ambient << 0.05, 0.05, 0.05.;
//Eigen::Vector3d specular; specular << 0, 0, 0, 1.;// specular << 0.1,0.1,0.1;
