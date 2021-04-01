#include "MeshLabelling.h"

#include <stack>

#include "GlobalSettings.h"
#include "Utils.h"
#include "Logger.h"

/*
bool label::try_get_flat_area(const int& source_vertex, const Mesh& mesh, std::vector<int>& flat_vertices, int& center_vertex)
{
	const int loglevel = 5;

	flat_vertices.clear();
	flat_vertices.push_back(source_vertex);

	center_vertex = -1;

	//check one ring neighborhood first, to see if starts at crease
	auto vertex_normal = mesh.normals_vertices.row(source_vertex);
	write_log(loglevel - 1) << linebreak << "try get FLAT area at " << source_vertex << ", normal = " << vertex_normal << std::endl;

	for (int fi : mesh.adjacency_VF[source_vertex])
	{
		auto n = mesh.normals_faces.row(fi);
		double a = angle(n, vertex_normal);
		write_log(loglevel - 1) << a << ",";

		if (a > GlobalSettings::flat_angle_threshold)
		{
			write_log(4) << "try get flat surface: not valid b/c starts at crease" << std::endl;
			return false;
		}
	}
	write_log(loglevel - 1) << "does not start at crease, proceed with labeling" << std::endl;

	flat_vertices = label::dfs_flat_area(source_vertex, mesh.adjacency_VV, mesh.normals_vertices);
	if (flat_vertices.size() < 1)
	{
		write_log(4) << "try get flat surface: not valid b/c no labelling found" << std::endl;
		return false;
	}
	write_log(loglevel - 1) << "labelled " << flat_vertices.size() << "/" << mesh.V.rows() << " as 'flat'" << std::endl;

	//using the center vertex of the flat area for better coverage
	center_vertex = get_central_vertex(flat_vertices, mesh.V);
	return true;
}
*/

std::vector<int> label::dfs_flat_area(const int& start_index, const std::vector<std::vector<int>>& connectivity_list, const Eigen::MatrixXd& N)
{
	std::vector<int> component;

	// Initially mark all verices as not visited 
	std::vector<bool> visited(connectivity_list.size(), false);

	// Create a stack for DFS 
	std::stack<int> stack;

	// Push the current source node. 
	stack.push(start_index);
	Eigen::RowVectorXd mean_normal = N.row(start_index);


	int index;

	while (!stack.empty())
	{
		// Pop a vertex from stack 
		index = stack.top();
		stack.pop();

		// Stack may contain same vertex twice. So we need to print the popped item only if it is not visited. 
		if (!visited[index])
		{
			visited[index] = true;

			Eigen::RowVectorXd normal = N.row(index);
			double a = angle(normal.transpose(), mean_normal.transpose());
			if (a <= GlobalSettings::flat_angle_threshold)
			{
				component.push_back(index);
				mean_normal = (mean_normal + normal).normalized();
			}
			else
			{
				continue;
			}
		}

		// Get all adjacent vertices of the popped vertex s 
		// If a adjacent has not been visited, then push it to the stack. 
		for (auto i = connectivity_list[index].begin(); i != connectivity_list[index].end(); ++i)
		{
			if (!visited[*i])
				stack.push(*i);
		}
	}

	return component;
}
