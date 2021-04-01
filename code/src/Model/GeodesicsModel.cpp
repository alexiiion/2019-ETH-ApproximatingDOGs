#include "GeodesicsModel.h"

#include "Logger.h"

using namespace std;



void GeodesicsModel::clear()
{
	for (Geodesic* geodesic : geodesics)
		geodesic->clear();
}

void GeodesicsModel::initialize()
{
	if (geodesics.size() > 0)
		clear();

	//initialize list of geodesics
	for (int i = 0; i < target.V.rows(); i++)
		geodesics.push_back(new Geodesic(i));

}

Geodesic* GeodesicsModel::at(int vertex)
{
	if (!geodesics[vertex]->is_traced())
		trace_at(vertex);

	return geodesics[vertex];
}

int GeodesicsModel::number_geodesics()
{
	return geodesics.size();
}

void GeodesicsModel::trace_at(int vertex)
{
	Eigen::RowVector3d k_max = target.surface_features().principal_k_max_direction.row(vertex);
	Eigen::RowVector3d direction = k_max;

	Eigen::MatrixXd current_geodesic;
	double current_length;

	walker.get_geodesic_at(vertex, direction, current_geodesic, current_length);
	write_log(5) << "v" << vertex << " geodesic.rows(): " << current_geodesic.rows() << ", current_length: " << current_length << endl << endl;
	
	geodesics[vertex]->path(current_geodesic);
}
