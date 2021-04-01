#pragma once

#include <Eigen/Core>
#include <vector>

#include "CurveHelper.h"

class Geodesic
{
public:

	Geodesic(int source_vertex) : _source_vertex(source_vertex) {};
	//Geodesic(Eigen::MatrixXd& path, double length, int source_vertex)
	//	: _path(path), _length(length), _source_vertex(source_vertex)
	//{
	//	_is_flat = curve::is_flat(_path);
	//};
	//Geodesic(Eigen::MatrixXd& path, int source_vertex)
	//	: _path(path), _source_vertex(source_vertex)
	//{
	//	_is_flat = curve::is_flat(_path);
	//	_length = curve::compute_length(_path);
	//};

	
	/* Properties (borrowed from C#) */

	int source_vertex() { return _source_vertex; };
	bool is_traced() { return _is_traced; };

	Eigen::MatrixXd path() { return _path; };
	void path(Eigen::MatrixXd& value) 
	{ 
		_path = value;
		_is_traced = true;
	};

	double length()
	{
		if (_length == -1 && _path.rows() > 1)
			_length = curve::compute_length(_path);
		return _length;
	};
	void length(double value) {	_length = value; };
	
	bool is_flat()
	{
		if(_is_flat == -1)
			_is_flat = curve::is_flat(_path);
		return _is_flat;
	}
	
	//might be that due to the current settings & stopping conditions, no geodesic can be traced at this vertex
	bool is_valid() { return _is_removed && _is_traced && _path.rows() > 2; }; 


	/* Methods */
	
	void remove() 
	{
		_is_removed = true; 
	};

	void clear()
	{
		_path.resize(0, 0);
		_length = -1;
		_is_flat = -1; 
		_is_traced = false;
	};


private:

	int _source_vertex = -1;
	Eigen::MatrixXd _path;
	double _length = -1;

	bool _is_valid = true;
	bool _is_traced = false;
	int _is_flat = -1; // -1 = "undefined", 0 = "false", 1 = "true"

	bool _is_removed = false;
};
