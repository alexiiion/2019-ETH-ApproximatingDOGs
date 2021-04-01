#include "RuledSurfacesModel.h"

#include "Utils.h"
#include "MeshController.h"
#include "MeshLabelling.h"

#include "Logger.h"

using namespace std;


void RuledSurfacesModel::clear()
{
	vertex_indices.clear();
	flat_vertex_indices.clear();
	selected_vertex_indices.clear();
	vertex_ruled_labelling.clear(); 

	for (auto surface : surfaces)
		delete surface;
	surfaces.clear();
}

bool RuledSurfacesModel::try_get_surface_at(int& vertex, RuledDevelopableSurface*& out_surface)
{
	if (vertex < 0 || vertex >= geodesics.number_geodesics())
		return false;

	Geodesic* geodesic = geodesics.at(vertex);
	if (!geodesic->is_valid())
		return false;

	if (geodesic->is_flat())
		return try_get_flat_surface_at(vertex, out_surface);

	return try_get_curved_surface_at(vertex, out_surface);
}

bool RuledSurfacesModel::try_get_curved_surface_at(const int& vertex, RuledDevelopableSurface*& out_surface)
{
	Geodesic* geodesic = geodesics.at(vertex);
	if (geodesic->path().rows() < 6) //otherwise can't create spline of degree 5
		return false;

	out_surface = new RuledDevelopableSurface(target, *&(geodesic->path()), vertex, RulingsType::Analytic);
	out_surface->create(settings.default_width);

	if (out_surface->developable.F.rows() < 10) //TODO extract parameter!
		return false;

	return true;
}

bool RuledSurfacesModel::try_get_flat_surface_at(int& vertex, RuledDevelopableSurface*& out_surface)
{
	const int loglevel = 5;

	//check one ring neighborhood first, to see if starts at crease
	auto vertex_normal = target.normals_vertices.row(vertex);
	for (int fi : target.adjacency_VF[vertex])
	{
		auto n = target.normals_faces.row(fi);
		if (angle(n, vertex_normal) > GlobalSettings::flat_angle_threshold)
			return false;
	}


	vector<int> flat_vertices = label::dfs_flat_area(vertex, target.adjacency_VV, target.normals_vertices);
	if (flat_vertices.size() < 1)
		return false;

	//using the center vertex of the flat area for better coverage
	int center_vertex = get_central_vertex(flat_vertices, target.V);
	int original_vertex = vertex; //for debug only
	vertex = center_vertex;

	//invalidate all other vertices
	for (int vi : flat_vertices)
	{
		if (vi != vertex)
			geodesics.at(vi)->remove();
	}

	Geodesic* geodesic = geodesics.at(vertex);
	if (!geodesic->is_valid())
		return false;
	if (geodesic->path().rows() < 3) //otherwise can't create spline of degree 5
		return false;

	//compute width of ruled surface based on size of flat area
	vector<int> flat_faces = meshhelper::get_faces_from_vertices(flat_vertices, target.adjacency_VF);
	double area = meshhelper::compute_face_area(flat_faces, target.V, target.F);
	double length = geodesic->length();
	double width = (area * settings.area_scale) / length;

	//create surface
	out_surface = new RuledDevelopableSurface(target, *&(geodesic->path()), vertex, RulingsType::Discrete);
	out_surface->create(width);
	
	//check that ruled surface is valid
	if (out_surface->developable.F.rows() < 2)
		return false;

	flat_vertex_indices.push_back(vertex);

	//debug output
	if (loglevel <= LOG_LEVEL) out_surface->print_mathematica_data();
	write_log(loglevel) << "flat_vertices: " << list_to_string(flat_vertices) << endl;
	write_log(loglevel) << "(v" << original_vertex << " to center v" << vertex << ") area: " << area << ", length: " << length << ", width = (area * 1.5) / length: " << width << endl;
	//write_log(loglevel) << "F_areas: " << endl << F_areas << endl;
	//write_log(loglevel) << "flat_faces: "; log_list(loglevel, flat_faces, "", false);

	return true;
}
