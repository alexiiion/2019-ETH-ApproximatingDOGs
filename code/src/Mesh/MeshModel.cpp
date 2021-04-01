#include "MeshModel.h"

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>

using namespace std;


Eigen::MatrixXd& NewMesh::V()
{
	return _V;
}

void NewMesh::V(Eigen::MatrixXd& value)
{
	_V = value;
	update();
}

Eigen::MatrixXi& NewMesh::F()
{
	return _F;
}
void NewMesh::F(Eigen::MatrixXi& value)
{
	_F = value;
	update();
}

Eigen::MatrixXd& NewMesh::NV()
{ 
	if(_NV.rows() < 1)
		igl::per_vertex_normals(_V, _F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, _NV);

	return _NV;
}

Eigen::MatrixXd& NewMesh::NF()
{ 
	if (_NF.rows() < 1)
		igl::per_face_normals(_V, _F, _NF);

	return _NF;
}

std::vector<std::vector<int>>& NewMesh::adjacency_VV()
{
	if(_adjacency_VV.size() < 1)
		igl::adjacency_list(_F, _adjacency_VV);

	return _adjacency_VV;
}

std::vector<std::vector<int>>& NewMesh::adjacency_VF()
{ 
	if (_adjacency_VF.size() < 1)
	{ 
		std::vector<std::vector<int>> VFi_unused;
		igl::vertex_triangle_adjacency(_V, _F, _adjacency_VF, VFi_unused);
	}

	return _adjacency_VF;
}
Eigen::MatrixXi& NewMesh::adjacency_FF()
{
	if (_adjacency_FF.size() < 1)
			igl::triangle_triangle_adjacency(_F, _adjacency_FF);

	return _adjacency_FF;
}

SurfaceFeatures& NewMesh::surface_features()
{
	//if (!_surface_features.are_built)
	//	_surface_features.build(*this);

	return _surface_features;
};


void NewMesh::update()
{ 
	if(!_keep_NV)
		_NV.resize(0, Eigen::NoChange);
	_NF.resize(0, Eigen::NoChange);

	_adjacency_VV.clear();
	_adjacency_VF.clear();
	_adjacency_FF.resize(0, Eigen::NoChange);

	_surface_features.reset();
}