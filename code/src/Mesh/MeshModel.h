#pragma once

#include "Quad.h"
#include "SurfaceFeatures.h"



class NewMesh 
{
public:

	NewMesh() { };
	NewMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) : _V(V), _F(F) { };
	NewMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& NV) : _V(V), _F(F), _NV(NV) 
	{
		_keep_NV = NV.rows() > 1;
	};

	virtual ~NewMesh() {}


	Eigen::MatrixXd& V();
	void V(Eigen::MatrixXd& value);

	Eigen::MatrixXi& F();
	void F(Eigen::MatrixXi& value);

	Eigen::MatrixXd& NV();
	Eigen::MatrixXd& NF();

	std::vector<std::vector<int>>& adjacency_VV();
	std::vector<std::vector<int>>& adjacency_VF();
	Eigen::MatrixXi& adjacency_FF();

	SurfaceFeatures& surface_features();


private:
	Eigen::MatrixXd _V;
	Eigen::MatrixXi _F;

	Eigen::MatrixXd _NV;
	Eigen::MatrixXd _NF;

	std::vector<std::vector<int>> _adjacency_VV;
	std::vector<std::vector<int>> _adjacency_VF;
	Eigen::MatrixXi _adjacency_FF;

	SurfaceFeatures _surface_features;
	
	bool _keep_NV;
	void update();

};

class Mesh 
{
public:
	virtual ~Mesh() {}

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	Eigen::MatrixXd normals_vertices;
	Eigen::MatrixXd normals_faces;

	std::vector<std::vector<int>> adjacency_VV;
	std::vector<std::vector<int>> adjacency_VF;
	Eigen::MatrixXi adjacency_FF;

	double average_edge_length = 0.0;
	double surface_area = 0.0;

	SurfaceFeatures& surface_features() 
	{
		if (!_surface_features.are_built)
			_surface_features.build(*this);
		return _surface_features;
	};

private:
	SurfaceFeatures _surface_features;

};

////TODO remove
//class Mesh : Mesh
//{
//public:
//	bool is_loaded = false;
//	//double scale = 1.0;
//};

class WrapperMesh : public Mesh
{
public:
	virtual ~WrapperMesh() {}

	Eigen::MatrixXi F_quad;
	QuadTopology quad_topology;
	Eigen::VectorXd initial_x0;
	
	int quad_width = 0;
	int quad_height = 0;
};
