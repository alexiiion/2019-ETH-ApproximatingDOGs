#include "GeodesicSurfacesController.h"

#include <igl/Timer.h>

#include "MeshModel.h"
#include "GeodesicsModel.h"
#include "RuledSurfacesModel.h"

#include "Utils.h"
#include "Logger.h"

using namespace std;


void GeodesicSurfacesController::clear()
{
	ruled_surfaces.clear();
	geodesics.clear();
}

void GeodesicSurfacesController::build_surfaces(int n)
{
	write_log(3) << endl << "create ~" << n << " random ruled developable surfaces..." << endl;

	igl::Timer timer;
	double init_time = timer.getElapsedTime();

	//TODO check!!
	clear();
	
	int number_vertices = target.V.rows();
	if (n >= number_vertices)
		create_all_surfaces();
	else
		create_random_surfaces(n);

	ruled_surfaces.settings.number_ruled_surfaces = ruled_surfaces.surfaces.size();

	////TODO move to view!
	//Mesh concatenated_developables = meshhelper::concatenate_ruled_developables(data_model.ruled_developables);
	//meshhelper::add_mesh(data_model.viewer, concatenated_developables.V, concatenated_developables.F, false, false);
	//data_model.ruled_view_index = data_model.viewer.selected_data_index;


	double t = timer.getElapsedTime();
	write_log(3) << "...done creating ruled developables -- ELAPSED TIME: " << t - init_time << endl << endl;

	write_log(4) << linebreak << "  random indices for valid ruled developables: " << list_to_string(ruled_surfaces.vertex_indices) << linebreak << endl;
}

void GeodesicSurfacesController::create_all_surfaces()
{
	int number_vertices = target.V.rows();
	
	for(int i = 0; i < number_vertices; i++)
		bool is_added = try_add_surface(i);
}

void GeodesicSurfacesController::create_random_surfaces(int n)
{
	int number_vertices = target.V.rows();
	vector<int> invalid_indices; // that previously gave invalid surface

	int added_indices = 0;
	int iterations = 0;
	int max_iterations = n * 3.0;

	while (added_indices < n && iterations < max_iterations)
	{
		iterations++;

		int random_index = std::rand() % number_vertices;

		if (is_contained(random_index, invalid_indices))
			continue;

		bool is_added = try_add_surface(random_index);
		if (is_added)
		{
			added_indices++;
			continue;
		}

		invalid_indices.push_back(random_index);
		write_log(5) << "invalid surface at " << random_index << endl;
	}
}

bool GeodesicSurfacesController::try_add_surface(int vertex)
{
	return false;

	if (is_contained(vertex, ruled_surfaces.vertex_indices))
		return false;

	RuledDevelopableSurface* surface;
	bool is_valid_surface = ruled_surfaces.try_get_surface_at(vertex, surface);


	if (!is_valid_surface)
		return false;

	ruled_surfaces.surfaces.push_back(surface);
	ruled_surfaces.vertex_indices.push_back(vertex);
}
