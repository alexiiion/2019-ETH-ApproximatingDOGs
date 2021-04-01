#include "GeodesicWalker.h"

#include <iostream>
#include <limits>
#include <tuple>
#include <fstream>
#include <iterator>

#include <igl/rotation_matrix_from_directions.h>
#include <igl/per_edge_normals.h>
#include <igl/barycentric_coordinates.h>

#include "CurveHelper.h"
#include "CoordinateConverter.h"
#include "Utils.h"
#include "Logger.h"

#include "GlobalSettings.h"

using namespace std;


double GeodesicWalker::frame_error(const int vertex, const int vertex_compare, bool use_weighted_error)
{
	Eigen::Vector3d n = target.normals_vertices.row(vertex);
	Eigen::Vector3d n_compare = target.normals_vertices.row(vertex_compare);
	Eigen::Vector3d k = target.surface_features().principal_k_max_direction.row(vertex).normalized(); 
	Eigen::Vector3d k_compare = target.surface_features().principal_k_max_direction.row(vertex_compare).normalized();

	const Eigen::Matrix3d rotation = igl::rotation_matrix_from_directions(n, n_compare);
	Eigen::Vector3d rotated_k = rotation * k;
	rotated_k.normalize();
	
	double error = 1 - abs(rotated_k.dot(k_compare));

	//if principal curvatures min & max are the same (difference 0), then their direction is not conclusive
	double weight = 1.0;
	if(use_weighted_error)
		weight = target.surface_features().principal_k_min(vertex_compare) * target.surface_features().principal_k_max(vertex_compare); //weight is gauss curvature, i.e., the flatter, the less important
		//weight = abs(target.surface_features.principal_k_min(vertex_compare) - target.surface_features.principal_k_max(vertex_compare)); 
	
	error *= weight;
	return error;
}

Eigen::VectorXd ortho_project(const Eigen::VectorXd projection_target, const Eigen::VectorXd to_project)
{
	//return to_project - (projection_target.dot(to_project) / projection_target.dot(projection_target) * projection_target);
	return (to_project.dot(projection_target) * projection_target) / projection_target.norm();
}

Eigen::Vector2d edgeIntersect(const Eigen::Vector2d& b, const Eigen::Vector2d& d, const Eigen::Vector2d& e0, const Eigen::Vector2d& e)
{
	const double det = -d(0) * e(1) + d(1) * e(0);

	if (abs(det) < 1e-8)
	{
		assert(0 && "no intersection found");
		return Eigen::Vector2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	}

	return Eigen::Vector2d((e0(1) - b(1)) * e(0) - (e0(0) - b(0)) * e(1), d(0) * (e0(1) - b(1)) - d(1) * (e0(0) - b(0))) / det;
}

int whichEdge(const Eigen::Vector2d& b, const Eigen::Vector2d& d, const Eigen::Matrix<double, 3, 2>& t, double& p)
{
	for (int i = 0; i < 3; ++i)
	{
		const Eigen::Vector2d e = t.row((i + 2) % 3) - t.row((i + 1) % 3);
		Eigen::Vector2d edg = edgeIntersect(b, d, t.row((i + 1) % 3), e);

		p = edg(1);

		if (edg(0) > 1e-6 && p > 0 && p < 1) 
			return i;

	}

	return -1;
}

Eigen::Vector2d triangeVertexCoordinates(const int eid, const Eigen::Matrix3d& t)
{
	const Eigen::Vector3d e2 = (t.row(eid) - t.row((eid + 1) % 3));
	const Eigen::Vector3d e = (t.row((eid + 2) % 3) - t.row((eid + 1) % 3)).normalized();

	const double c0 = e2.dot(e);

	return Eigen::Vector2d(c0, (c0 * e - e2).norm());
}

void GeodesicWalker::get_geodesic_kmax_at(const int vertex_index, Eigen::MatrixXd& out_geodesic, double& out_length)
{
	Eigen::RowVector3d direction = target.surface_features().principal_k_max_direction.row(vertex_index);
	get_geodesic_at(vertex_index, direction, out_geodesic, out_length);
}

void GeodesicWalker::get_geodesic_kmin_at(const int vertex_index, Eigen::MatrixXd& out_geodesic, double& out_length)
{
	Eigen::RowVector3d direction = target.surface_features().principal_k_min_direction.row(vertex_index);
	get_geodesic_at(vertex_index, direction, out_geodesic, out_length);
}

void GeodesicWalker::get_geodesic_at(const int vertex_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_geodesic, double& out_length)
{
	Eigen::MatrixXi vertices;
	Eigen::VectorXi faces;
	get_geodesic_at(vertex_index, direction, out_geodesic, out_length, vertices, faces);
}

Eigen::MatrixXd GeodesicWalker::merge_parts(const Eigen::MatrixXd& part_positive, const Eigen::MatrixXd& part_negative)
{
	if (part_positive.rows() < 1 && part_negative.rows() > 0)
		return part_negative;
	if (part_negative.rows() < 1 && part_positive.rows() > 0)
		return part_positive;
	if (part_negative.rows() < 1 && part_positive.rows() < 1)
		return Eigen::MatrixXd(0,0);

	int number_points = part_positive.rows() + part_negative.rows() - 1;
	Eigen::MatrixXd merged(number_points, part_positive.cols());

	Eigen::MatrixXd part_negative_reversed = part_negative.colwise().reverse();
	merged << part_negative_reversed.block(0, 0, part_negative_reversed.rows() - 1, part_negative_reversed.cols()),
			  part_positive;

	return merged;
}

Eigen::MatrixXi GeodesicWalker::merge_parts(const Eigen::MatrixXi& part_positive, const Eigen::MatrixXi& part_negative)
{
	int number_points = part_positive.rows() + part_negative.rows() - 1;
	Eigen::MatrixXi merged(number_points, part_positive.cols());

	Eigen::MatrixXi part_negative_reversed = part_negative.colwise().reverse();
	merged << part_negative_reversed.block(0, 0, part_negative_reversed.rows() - 1, part_negative_reversed.cols()),
			  part_positive;

	return merged;
}

Eigen::VectorXi GeodesicWalker::merge_parts(const Eigen::VectorXi& part_positive, const Eigen::VectorXi& part_negative)
{
	int number_points = part_positive.rows() + part_negative.rows() - 1;
	Eigen::VectorXi merged(number_points);

	Eigen::VectorXi part_negative_reversed = part_negative.colwise().reverse();
	merged << part_negative_reversed.block(0, 0, part_negative_reversed.rows() - 1, part_negative_reversed.cols()),
			  part_positive;

	return merged;
}

void GeodesicWalker::get_geodesic_at(const int vertex_index, const Eigen::RowVector3d& direction, Eigen::MatrixXd& out_geodesic, double& out_length, Eigen::MatrixXi& out_path_adjacent_vertices, Eigen::VectorXi& out_path_faces)
{
	out_geodesic.resize(0, 0);
	out_length = 0.0;

	source_vertex = vertex_index;

	Eigen::Vector3d dir(direction);

	Eigen::MatrixXd path_positive;
	Eigen::MatrixXi vertices_positive;
	Eigen::VectorXi faces_positive;
	double length_positive;

	write_log(6) << linebreak << "trace geodesic at " << source_vertex << " in POSITIVE direction" << endl;
	trace_geodesic(dir, path_positive, length_positive, vertices_positive, faces_positive);


	Eigen::MatrixXd path_negative;
	Eigen::MatrixXi vertices_negative;
	Eigen::VectorXi faces_negative;
	double length_negative;

	dir = -dir;
	write_log(6) << linebreak << "trace geodesic at " << source_vertex << " in NEGATIVE direction" << endl;
	trace_geodesic(dir, path_negative, length_negative, vertices_negative, faces_negative);


	if (path_positive.rows() + path_negative.rows() < 1)
		return;

	if (path_positive.rows() > 0 && path_negative.rows() <= 0)
	{
		out_length = length_positive;
		out_geodesic = path_positive;
		out_path_adjacent_vertices = vertices_positive;
		out_path_faces = faces_positive;
	}
	else if (path_positive.rows() <= 0 && path_negative.rows() > 0)
	{
		out_length = length_negative;
		out_geodesic = path_negative;
		out_path_adjacent_vertices = vertices_negative;
		out_path_faces = faces_negative;
	}
	else
	{

		bool starts_at_crease = false;

		if (path_positive.rows() > 0 && path_negative.rows() > 0)
		{
			auto normal_positive = get_surface_normal(path_positive.row(1), target.V, target.F, target.normals_vertices, target.adjacency_VF);
			auto normal_negative = get_surface_normal(path_negative.row(1), target.V, target.F, target.normals_vertices, target.adjacency_VF);
			//write_log(0) << "source normal: " << target.normals_vertices.row(source_vertex) << ", normal_positive: " << normal_positive.transpose() << ", normal_negative: " << normal_negative.transpose() << endl;

			double dot = normal_positive.dot(normal_negative);
			double cos_angle = clip(dot, -1, 1);
			double angle = acos(cos_angle);
			//write_log(0) << "dot: " << dot << ", cos_angle: " << cos_angle << ", angle: " << angle << endl;

			starts_at_crease = angle >= begin_crease_threshold;
		}

		if (starts_at_crease)
		{
			if (length_positive > length_negative)
			{
				out_length = length_positive;
				out_geodesic = path_positive;
				out_path_adjacent_vertices = vertices_positive;
				out_path_faces = faces_positive;
			}
			else
			{
				out_length = length_negative;
				out_geodesic = path_negative;
				out_path_adjacent_vertices = vertices_negative;
				out_path_faces = faces_negative;
			}
		}
		else
		{
			out_length = length_positive + length_negative;
			out_geodesic = merge_parts(path_positive, path_negative);
			out_path_adjacent_vertices = merge_parts(vertices_positive, vertices_negative);
			out_path_faces = merge_parts(faces_positive, faces_negative);
		}
	}

	double segment_length = out_length / out_geodesic.rows();
	Eigen::MatrixXd resampled_geodesic;
	double resampled_length;
	curve::resample_uniformly(out_geodesic, segment_length, resampled_geodesic, resampled_length);
	
	write_log(5) << "geodesic: length change after resampling = " << (resampled_length - out_length)  << ", one segment? " << segment_length - (resampled_length - out_length) << endl;

	out_geodesic = resampled_geodesic;
	out_length = resampled_length;


	write_log(5) << "joined geodesic: " << endl << out_geodesic << endl;
	write_log(6) << "path_positive: " << endl << path_positive << endl;
	write_log(6) << "path_negative: " << endl << path_negative << endl;
}

int GeodesicWalker::get_source_face(const Eigen::Vector3d& direction, const Eigen::Vector3d& e0, const Eigen::Vector3d& e1)
{
	int fid = -1;

	Eigen::Vector2d n2d(e0.dot(direction), e1.dot(direction));
	n2d.normalize();

	/* project adjacent triangles and find correct one */
	for (int fi : target.adjacency_VF[source_vertex])
	{
		int j = 0;
		for (; j < 3; ++j)
		{
			if (target.F(fi, j) == source_vertex)
				break;
		}

		if (j == 3)
		{
			assert(0);
			break;
		}

		Eigen::Vector3d v = target.V.row(target.F(fi, (j + 1) % 3)) - target.V.row(target.F(fi, j));
		Eigen::Vector2d v0(e0.dot(v), e1.dot(v));

		v = target.V.row(target.F(fi, (j + 2) % 3)) - target.V.row(target.F(fi, j));
		Eigen::Vector2d v1(e0.dot(v), e1.dot(v));

		v0.normalize();
		v1.normalize();

		// is the projected triangle flipped?
		if (v0(0) * v1(1) - v1(0) * v0(1) < 0)
			break;

		// we are in the right triangle if n2d splits it into two positively oriented triangles
		if (v0(0) * n2d(1) - n2d(0) * v0(1) > 0 &&
			n2d(0) * v1(1) - v1(0) * n2d(1) > 0)
		{
			fid = fi;
			break;
		}
	}

	return fid;


	Eigen::RowVector3d offset = direction.normalized() * 0.001;
	Eigen::RowVector3d point = target.V.row(source_vertex) + offset;
	
	int my_fid = -1;
	for (int fi : target.adjacency_VF[source_vertex])
	{
		Eigen::RowVectorXi face = target.F.row(fi);
		
		Eigen::RowVector3d barycentric_weights;
		igl::barycentric_coordinates(point, target.V.row(face(0)), target.V.row(face(1)), target.V.row(face(2)), barycentric_weights);

		if (barycentric_weights(0) >= 0.0 && barycentric_weights(0) <= 1.0 && barycentric_weights(1) >= 0.0 && barycentric_weights(1) <= 1.0 && barycentric_weights(2) >= 0.0 && barycentric_weights(2) <= 1.0)
		{
			my_fid = fi;
			break;
		}
	}
	//write_log(0) << "my_fid: " << my_fid << endl;
	//write_log(0) << "fid: " << fid << endl;

	return my_fid;
}

void GeodesicWalker::trace_geodesic(Eigen::Vector3d& dir, Eigen::MatrixXd& path, double& length, Eigen::MatrixXi& path_adjacent_vertices, Eigen::VectorXi& path_faces)
{
	WalkedData traced;
	stopping_criteria.reset();


	const int l = 6;
	write_log(l) << linebreak << "trace geodesic:: " << endl;



	Eigen::Matrix<double, 3, 2> curr;

	/* project 'dir' onto tangent plane and find correct triangle. */
	/* build tangent space basis */
	Eigen::Vector3d n = target.normals_vertices.row(source_vertex);
	write_log(l) << "  vertex [" << source_vertex << "] " << target.V.row(source_vertex) << ", normal = " << n.transpose() << endl;
	write_log(l) << "  direction = " << dir.transpose() << endl;

	Eigen::Vector3d e0 = abs(n(0)) > .9 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0);
	e0 -= e0.dot(n) * n;
	e0.normalize();
	Eigen::Vector3d e1 = n.cross(e0);
	write_log(l) << "  e0 = " << e0.transpose() << ", e1 = " << e1.transpose() << endl;



	//check if direction is (0,0,0), if so choose another direction 
	if (abs(dir.sum()) < 1e-6)
	{
		dir = 0.98*e0 + 0.02*e1; 
		write_log(l) << "    adapted direction = " << dir.transpose() << endl;
	}

	int fid = get_source_face(dir, e0, e1);
	if (fid == -1)
	{
		//write_log(3) << "warn at line: " << __LINE__ << ". no valid starting face found at vertex " << source_vertex << endl;
		return;
	}


	/* place first triangle */
	n = target.normals_faces.row(fid).normalized();
	e0 = e0 = abs(n(0)) > .9 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0);
	e0 -= e0.dot(n) * n;
	e0.normalize();
	e1 = n.cross(e0);

	write_log(l) << linebreak << "  first face id [" << fid << "] " << target.F.row(fid) << ", normal = " << n.transpose() << endl;
	write_log(l) << "  e0 = " << e0.transpose() << ", e1 = " << e1.transpose() << endl;

	for (int i = 0; i < 3; ++i)
	{
		curr(i, 0) = e0.dot(target.V.row(target.F(fid, i)) - target.V.row(target.F(fid, 0)));
		curr(i, 1) = e1.dot(target.V.row(target.F(fid, i)) - target.V.row(target.F(fid, 0)));
	}
	write_log(l) << "  triangle 'curr': " << linebreak << curr << endl;




	Eigen::Vector2d d(e0.dot(dir), e1.dot(dir));
	d.normalize();
	write_log(l) << "  d = " << d.transpose() << endl;


	int j = 0;
	for (; j < 3; ++j)
	{
		if (target.F(fid, j) == source_vertex) 
			break;
	}

	Eigen::Vector2d b = curr.row(j).transpose() + 1e-9 * d;
	write_log(l) << "  b = " << b.transpose() << endl;


	traced.path.push_back(target.V.row(source_vertex));
	traced.faces.push_back(fid);
	traced.vertices.push_back(pair<int, int>(source_vertex, source_vertex));

	bool is_done_walking = false;

	while (!is_done_walking)
	{
		double p;
		int eid = whichEdge(b, d, curr, p);
		write_log(l) << "    eid = " << eid << endl;

		if (eid == -1)
		{
			write_log(l) << "geodesic walker at " << source_vertex << ": edge ID not found!" << endl;
			break;
		}

		int fidOld = fid;
		//fid_previous = fidOld;

		fid = target.adjacency_FF(fidOld, (eid + 1) % 3);
		//fid_current = fid;
		write_log(l) << "    fid = " << eid << endl;

		if (fid == -1)
		{
			//write_log(3) << "warn line: " << __LINE__ << ". current face is invalid. previous face was " << fidOld << endl;
			
			//TODO this hits the border, ADD last point on border!!
			//path.row(added_segments) = (target.V.row(target.F(fidOld, (eid + 1) % 3)) + p * (target.V.row(target.F(fidOld, (eid + 2) % 3)) - target.V.row(target.F(fidOld, (eid + 1) % 3))));
			//added_segments++;
			break;
		}

		// layout new triangle
		int eid2 = -1;
		for (int k = 0; k < 3; ++k)
		{
			if (target.adjacency_FF(fid, (k + 1) % 3) == fidOld)
			{
				eid2 = k;
				break;
			}
		}

		//if (eid2 == -1) cout << "error line: " << __LINE__ << " hit a border" << endl;
		if (eid2 == -1)
		{
			write_log(3) << "info at line: " << __LINE__ << ". reached a border." << endl;

			//path.row(added_segments) = (target.V.row(target.F(previous_face, (eid + 1) % 3)) + p * (target.V.row(target.F(previous_face, (eid + 2) % 3)) - target.V.row(target.F(previous_face, (eid + 1) % 3))));
			//added_segments++;

			break;
		}

		Eigen::Matrix3d t2;
		for (int k = 0; k < 3; ++k)
			t2.row(k) = target.V.row(target.F(fid, k));

		Eigen::Vector2d coords = triangeVertexCoordinates(eid2, t2);

		Eigen::Vector2d p0 = curr.row((eid + 2) % 3);
		Eigen::Vector2d p1 = curr.row((eid + 1) % 3);

		Eigen::Vector2d edge = (p1 - p0).normalized();
		curr.row(eid2) = p0 + coords(0) * edge + coords(1) * Eigen::Vector2d(-edge(1), edge(0));
		curr.row((eid2 + 1) % 3) = p0;
		curr.row((eid2 + 2) % 3) = p1;

		b = p1 + p * (p0 - p1);
		auto point = target.V.row(target.F(fidOld, (eid + 1) % 3)) + p * (target.V.row(target.F(fidOld, (eid + 2) % 3)) - target.V.row(target.F(fidOld, (eid + 1) % 3)));

		auto last_point = traced.path.back();
		traced.path.push_back(point);
		traced.length += (last_point - point).norm();

		int edge_vertex_index1 = target.F(fid, (eid2 + 1) % 3);
		int edge_vertex_index2 = target.F(fid, (eid2 + 2) % 3);
		traced.vertices.push_back(pair<int,int>(edge_vertex_index1, edge_vertex_index2));
		traced.faces.push_back(fid);

		stopping_criteria.update(traced, target);

		is_done_walking = stopping_criteria.should_stop();
	}

	////TODO add last point

	//write_log(0) << "  done tracing:: faces:" << linebreak << list_to_string(traced.faces) << linebreak << " vertices: " << linebreak;
	//for (int i = 0; i < traced.vertices.size(); i++)
	//	write_log(0) << "(" << traced.vertices[i].first << ", " << traced.vertices[i].second << ")" << linebreak;
	//write_log(0) << endl;


	length = traced.length;
	
	path.resize(traced.path.size(), 3);
	for (int i = 0; i < traced.path.size(); i++)
		path.row(i) = traced.path[i];

	path_faces = Eigen::Map<Eigen::VectorXi>(traced.faces.data(), traced.faces.size());

	path_adjacent_vertices.resize(traced.vertices.size(), 2);
	for (int i = 0; i < traced.vertices.size(); i++)
	{
		path_adjacent_vertices(i, 0) = traced.vertices[i].first;
		path_adjacent_vertices(i, 1) = traced.vertices[i].second;
	}
}

/*
void GeodesicWalker::trace_geodesic(Eigen::Vector3d& dir, Eigen::MatrixXd& path, double& length, Eigen::MatrixXi& path_adjacent_vertices, Eigen::VectorXi& path_faces)
{
	path.resize(target.F.rows(), 3);
	path_adjacent_vertices.resize(target.F.rows(), 2);
	path_faces.resize(target.F.rows());
	length = 0.0;
	
	//reset_stopping_criteria_data();
	//
	////TODO
	//cluster.reset();
	//cluster.absolute_threshold = GlobalSettings::crease_angle_threshold;
	//cluster.outlier_factor = 3.0; //data_model.stopping_outlier_factor
	//cluster.window_size = 5; //data_model.stopping_window_size

	stopping_criteria.reset();


	Eigen::Matrix<double, 3, 2> curr;

	// project 'dir' onto tangent plane and find correct triangle.
	// build tangent space basis 
	Eigen::Vector3d n = target.normals_vertices.row(source_vertex);
	Eigen::Vector3d e0 = n(0) > .9 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0);
	e0 -= e0.dot(n) * n;
	e0.normalize();
	Eigen::Vector3d e1 = n.cross(e0);

	//Eigen::Vector3d n = target.normals_vertices.row(source_vertex);
	//Eigen::RowVector3d neighbor_vertex = target.V.row(target.adjacency_VV[source_vertex][0]);
	//Eigen::Vector3d e0 = (neighbor_vertex - target.V.row(source_vertex)).transpose().normalized();
	//Eigen::Vector3d e1 = n.cross(e0);

	//check if direction is (0,0,0), if so choose another direction 
	if (abs(dir.sum()) < 1e-6)
		dir = 0.95*e0 + 0.05*e1; 

	int fid = get_source_face(dir, e0, e1);
	if (fid == -1)
	{
		write_log(3) << "warn at line: " << __LINE__ << ". no valid starting face found at vertex " << source_vertex << endl;
		path.resize(0, 0);
		
		return;
	}

	// place first triangle 
	n = target.normals_faces.row(fid).normalized();
	e0 = n(0) > .9 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0);
	e0 -= e0.dot(n) * n;
	e0.normalize();
	e1 = n.cross(e0);

	for (int i = 0; i < 3; ++i)
	{
		curr(i, 0) = e0.dot(target.V.row(target.F(fid, i)) - target.V.row(target.F(fid, 0)));
		curr(i, 1) = e1.dot(target.V.row(target.F(fid, i)) - target.V.row(target.F(fid, 0)));
	}

	Eigen::Vector2d d(e0.dot(dir), e1.dot(dir));
	d.normalize();

	int j = 0;
	for (; j < 3; ++j)
	{
		if (target.F(fid, j) == source_vertex) 
			break;
	}

	Eigen::Vector2d b = curr.row(j).transpose() + 1e-9 * d;
	path.row(added_segments) = target.V.row(source_vertex);

	path_faces(added_segments) = fid;
	path_adjacent_vertices(added_segments, 0) = source_vertex;
	path_adjacent_vertices(added_segments, 1) = source_vertex;

	added_segments++;

	//fid_current = fid;
	//double angle_sum = 0.0;
	//bool is_done_winding = false;
	//bool is_done_walking_ = false;
	//
	//while (!is_done_walking_ && !is_done_winding)

	bool is_done_walking = false;
	while (!is_done_walking)
	{
		if (added_segments >= path.rows())
			break;


		double p;
		int eid = whichEdge(b, d, curr, p);

		if (eid == -1)
			break;

		int fidOld = fid;
		//fid_previous = fidOld;

		fid = target.adjacency_FF(fidOld, (eid + 1) % 3);
		//fid_current = fid;

		if (fid == -1)
		{
			write_log(3) << "warn line: " << __LINE__ << ". current face is invalid. previous face was " << fidOld << endl;
			
			//TODO this hits the border, ADD last point on border!!
			//path.row(added_segments) = (target.V.row(target.F(fidOld, (eid + 1) % 3)) + p * (target.V.row(target.F(fidOld, (eid + 2) % 3)) - target.V.row(target.F(fidOld, (eid + 1) % 3))));
			//added_segments++;
			break;
		}

		// layout new triangle
		int eid2 = -1;
		for (int k = 0; k < 3; ++k)
		{
			if (target.adjacency_FF(fid, (k + 1) % 3) == fidOld)
			{
				eid2 = k;
				break;
			}
		}

		//if (eid2 == -1) cout << "error line: " << __LINE__ << " hit a border" << endl;
		if (eid2 == -1)
		{
			write_log(3) << "info at line: " << __LINE__ << ". reached a border." << endl;

			//path.row(added_segments) = (target.V.row(target.F(previous_face, (eid + 1) % 3)) + p * (target.V.row(target.F(previous_face, (eid + 2) % 3)) - target.V.row(target.F(previous_face, (eid + 1) % 3))));
			//added_segments++;

			break;
		}

		Eigen::Matrix3d t2;
		for (int k = 0; k < 3; ++k)
			t2.row(k) = target.V.row(target.F(fid, k));

		Eigen::Vector2d coords = triangeVertexCoordinates(eid2, t2);

		Eigen::Vector2d p0 = curr.row((eid + 2) % 3);
		Eigen::Vector2d p1 = curr.row((eid + 1) % 3);

		Eigen::Vector2d edge = (p1 - p0).normalized();
		curr.row(eid2) = p0 + coords(0) * edge + coords(1) * Eigen::Vector2d(-edge(1), edge(0));
		curr.row((eid2 + 1) % 3) = p0;
		curr.row((eid2 + 2) % 3) = p1;

		b = p1 + p * (p0 - p1);
		auto point = target.V.row(target.F(fidOld, (eid + 1) % 3)) + p * (target.V.row(target.F(fidOld, (eid + 2) % 3)) - target.V.row(target.F(fidOld, (eid + 1) % 3)));

		path.row(added_segments) = point;
		length += (path.row(added_segments - 1) - point).norm();


		//update_stopping_criteria(path);
		stopping_criteria.update(path, target, added_segments);


		path_faces(added_segments) = fid;

		int edge_vertex_index1 = target.F(fid, (eid2 + 1) % 3);
		int edge_vertex_index2 = target.F(fid, (eid2 + 2) % 3);
		path_adjacent_vertices(added_segments, 0) = edge_vertex_index1;
		path_adjacent_vertices(added_segments, 1) = edge_vertex_index2;

		//write_log(4) << "fid = " << fid << ", fidOld = " << fidOld << ", eid = " << eid << ", eid2 = " << eid2 << endl;
		//write_log(4) << "target.F.row(fid): " << target.F.row(fid) << endl;
		//write_log(4) << "target.F.row(fidOld): " << target.F.row(fidOld) << endl;
		//write_log(4) << "target.adjacency_FF.row(fid): " << target.adjacency_FF.row(fid) << endl;
		//write_log(4) << "target.adjacency_FF.row(fidOld): " << target.adjacency_FF.row(fidOld) << endl;

		added_segments++;

		////update stopping conditions
		//angle_sum += last_deviation;
		//is_done_winding = angle_sum >= igl::PI;
		//write_log(5) << "is_done_winding? " << boolalpha << is_done_winding << ", angle_sum: " << angle_sum << ", last_deviation: " << last_deviation << endl;
		//
		//is_done_walking_ = is_done_walking();

		is_done_walking = stopping_criteria.should_stop();
	}

	//TODO add last point

	path.conservativeResize(added_segments, Eigen::NoChange);
	path_faces.conservativeResize(added_segments, Eigen::NoChange);
	path_adjacent_vertices.conservativeResize(added_segments, Eigen::NoChange);
}
*/

/*
bool GeodesicWalker::is_done_walking()
{
	if (stopping_mode == WalkingStop::NormalsDeviation)
		return current_normals_deviation >= max_normals_deviation;
	else if (stopping_mode == WalkingStop::FrameDeviation)
		return error_derivative >= *max_frame_error || added_segments >= 50;
	else if (stopping_mode == WalkingStop::MaxLength)
		return current_length >= *max_length / 2;
	else if (stopping_mode == WalkingStop::AdaptiveNormalsDeviation)
		return cluster.is_value_within(last_deviation) == false;
}

void GeodesicWalker::update_stopping_criteria(const Eigen::MatrixXd& path)
{
	auto current_point = path.row(added_segments);
	auto previous_point = path.row(added_segments-1);

	if (stopping_mode == WalkingStop::NormalsDeviation)
	{
		//TODO use barycentric normal instead of closest vertices

		//previous_normals_deviation = current_normals_deviation;
		int previous_closest_vertex = get_closest_vertex(target.V, previous_point);
		int closest_vertex = get_closest_vertex(target.V, current_point);
		current_normals_deviation = normals_deviation(closest_vertex, previous_closest_vertex, target.normals_vertices);
	}
	else if (stopping_mode == WalkingStop::FrameDeviation)
	{
		int previous_closest_vertex = get_closest_vertex(target.V, previous_point);
		int closest_vertex = get_closest_vertex(target.V, current_point);

		if (closest_vertex != previous_closest_vertex)
		{
			double error = frame_error(closest_vertex, previous_closest_vertex, use_weighted_frame_error);
			error_integral += error;
			error_derivative = abs(error_previous - error);
			error_previous = error;
		}
	}
	else if (stopping_mode == WalkingStop::MaxLength)
	{
		current_length += (current_point - previous_point).norm();
	}
	else if (stopping_mode == WalkingStop::AdaptiveNormalsDeviation)
	{
		last_deviation = normals_deviation(fid_current, fid_previous, target.normals_faces);
	}
}

void GeodesicWalker::reset_stopping_criteria_data()
{
	added_segments = 0;

	error_previous = 0.0;
	error_derivative = 0.0;
	error_integral = 0.0;

	current_length = 0.0; 
	current_normals_deviation = 0.0; 

	last_deviation = 0.0;
	fid_previous = -1;
	fid_current = -1;
}
*/