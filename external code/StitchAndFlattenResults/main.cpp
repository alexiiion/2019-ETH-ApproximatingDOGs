#include <igl/opengl/glfw/Viewer.h>
#include <igl/readCSV.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/facet_components.h>
#include <igl/vertex_components.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_lengths.h>
#include <igl/avg_edge_length.h>
#include <igl/dijkstra.h>
#include <igl/remove_duplicate_vertices.h>

#include <igl/slim.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/MappingEnergyType.h>
#include <igl/flipped_triangles.h>

#include <stack>

#include "Utils.h"

using namespace std;
using namespace Eigen;

struct Mesh
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
};

string folder_name;
string model_name;
string output_prefix;

void stitch_mesh(const string& folder_name, const string& model_name);
void flatten_patches(const string& folder_name, const string& model_name);

std::string list_to_string(const std::vector<int>& list, const std::string delimiter = " ", const bool multi_line = false, const bool print_index = false)
{
	std::stringstream output;
	output << std::fixed;

	for (int i = 0; i < list.size(); i++)
	{
		if (print_index)
			output << "[" << i << "]: ";

		output << list[i];

		if (i >= list.size() - 1)
			continue;

		if (multi_line)
			output << '\n';
		else
			output << delimiter;
	}

	return output.str();
}

int index_of(const int& element, const Eigen::VectorXi& list)
{
	for (int i = 0; i < list.size(); i++)
		if (list(i) == element)
			return i;

	return -1;
}

vector<int> get_edge_incident_faces(int v1, int v2, const std::vector<std::vector<int>>& adjacency_VF)
{
	vector<int> face_current = adjacency_VF[v1];
	vector<int> face_next = adjacency_VF[v2];
	vector<int> edge_faces = intersect(face_current, face_next);

	//cout << "  VV(" << vi_current << "):   " << list_to_string(adjacency_VV[vi_current]) << endl;
	//cout << "  VV(" << vi_next << "):   " << list_to_string(adjacency_VV[vi_next]) << endl;
	//cout << "  face_current: " << list_to_string(face_current) << endl;
	//cout << "  face_next:    " << list_to_string(face_next) << endl;

	//cout << "  edge_faces:   " << list_to_string(edge_faces) << endl;
	//cout << "     F.row(" << edge_faces[0] << "):   " << F.row(edge_faces[0]) << endl;
	//cout << "     F.row(" << edge_faces[1] << "):   " << F.row(edge_faces[1]) << endl;

	return edge_faces;
}

//void get_oriented_incident_face_at(const int& index, const Eigen::MatrixXd& cut_path, const Eigen::MatrixXi& F, const std::vector<std::vector<int>>& adjacency_VF, int& out_valid_face, int& out_unused_face)
//{
//	if (index >= cut_path.rows() - 1)
//		return;
//
//	int vi_current = cut_path[index];
//	int vi_next = cut_path[index + 1];
//	cout << endl << vi_current << " - " << vi_next << endl;
//
//	vector<int> cutline_faces = get_edge_incident_faces(vi_current, vi_next, adjacency_VF);
//	cout << "  cutline_faces(" << vi_current << "): " << list_to_string(cutline_faces) << endl;
//
//	bool has_found_face = false;
//	int fi = 0;
//
//	for (int i = 0; i < cutline_faces.size() && !has_found_face; i++)
//	{
//		int fi = cutline_faces[i];
//		cout << "    F.row(" << fi << "): " << F.row(fi) << endl;
//
//		int vi1 = index_of(vi_current, F.row(fi));
//		int vi2 = (vi1 + 1) % 3;
//
//		if (F(fi, vi2) != vi_next)
//			continue;
//
//		has_found_face = true;
//	}
//
//	out_valid_face = fi;
//	out_
//}

int get_oriented_incident_face_at(const int& query_face, const vector<int>& cut_path_indices, const Eigen::MatrixXi& F, vector<int>& contained_cut_indices)
{
	contained_cut_indices.clear();

	for (int k = 0; k < F.row(query_face).size(); k++)
	{
		int vi_f = F(query_face, k);
		int ci_f = index_of(vi_f, cut_path_indices);
		if (ci_f >= 0)
			contained_cut_indices.push_back(ci_f);
	}

	if (contained_cut_indices.size() != 2)
		return -1;


	int i1 = min(contained_cut_indices[0], contained_cut_indices[1]);
	int i2 = max(contained_cut_indices[0], contained_cut_indices[1]);
	int v1 = cut_path_indices[i1];
	int v2 = cut_path_indices[i2];

	//cout << "    F.row(" << query_face << "): " << F.row(query_face) << endl;

	int vi1 = index_of(v1, F.row(query_face));
	int vi2 = (vi1 + 1) % 3;

	if (F(query_face, vi2) == v2)
		return query_face;

	return -1;
}

int get_oriented_incident_face_at(vector<int> cutline_faces, const int& v1, const int& v2, const Eigen::MatrixXi& F)
{
	if (cutline_faces.size() < 1)
		return -1;

	int fi = -1;
	for (int i = 0; i < cutline_faces.size(); i++)
	{
		fi = cutline_faces[i];
		//cout << "    F.row(" << fi << "): " << F.row(fi) << endl;

		int vi1 = index_of(v1, F.row(fi));
		int vi2 = (vi1 + 1) % 3;

		if (F(fi, vi2) == v2)
			return fi;
	}

	return -1;
}


Mesh sub_mesh(std::vector<int>& vertex_indices, std::vector<int>& face_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original)
{
	vertex_indices = remove_duplicates(vertex_indices);

	Eigen::MatrixXi mapped_faces(face_indices.size(), 3);
	int added_faces = 0;

	for (int i = 0; i < face_indices.size(); i++)
	{
		int fi = face_indices[i];
		Eigen::RowVectorXi face_vertices_global = F_original.row(fi);
		//write_log(4) << i << ": global face_vertices for fi=" << fi << ": " << face_vertices_global << endl;

		vector<int> face_vertices_local;
		for (int j = 0; j < face_vertices_global.size(); j++)
		{
			int vi = face_vertices_global(j);
			face_vertices_local.push_back(index_of(vi, vertex_indices));
		}

		if (is_contained(-1, face_vertices_local))
		{
			//write_log(4) << "   NOT including this face!" << endl;
			continue;
		}

		mapped_faces.row(added_faces) = Eigen::Map<Eigen::RowVectorXi>(face_vertices_local.data(), face_vertices_local.size());
		added_faces++;

		//write_log(4) << "local component_faces.row(" << i << "): " << component_faces.row(i) << endl;
	}

	mapped_faces.conservativeResize(added_faces, 3);
	//write_log(4) << "local component_faces: " << endl << component_faces << endl;

	Eigen::MatrixXd V;
	index_to_value(vertex_indices, V_original, V);

	Mesh mesh;
	mesh.V = V;
	mesh.F = mapped_faces;

	return mesh;
}

Mesh sub_mesh(std::vector<int>& vertex_indices, const Eigen::MatrixXd& V_original, const Eigen::MatrixXi& F_original, const std::vector<std::vector<int>>& adjacency_VF)
{
	vector<int> face_indices;
	for each (int vertex in vertex_indices)
	{
		vector<int> faces = adjacency_VF[vertex];
		face_indices.insert(face_indices.end(), faces.begin(), faces.end());
	}

	std::sort(face_indices.begin(), face_indices.end());
	face_indices.erase(std::unique(face_indices.begin(), face_indices.end()), face_indices.end());


	Mesh mesh = sub_mesh(vertex_indices, face_indices, V_original, F_original);
	return mesh;
}

void read_config()
{
	std::cout << endl << "READING config settings: " << endl << endl;

	//get config file path
	//string file = __FILE__;
	//auto index = file.find_last_of("/\\");
	//auto path = file.substr(0, index + 1);
	//auto config_filepath = path + "config.txt";
	string config_filepath = "./config.txt";

	// Open the File
	std::ifstream in(config_filepath.c_str());

	// Check if object is valid
	if (!in)
	{
		std::cerr << "Cannot open the config file : " << config_filepath << std::endl;
		exit(1);
	}

	std::string line;
	// Read the next line from File untill it reaches the end.
	while (std::getline(in, line))
	{
		if (line.find('#') == 0)
			continue;
		if (line.empty())
			continue;

		auto colon_index = line.find_first_of(":");
		string qualifier = line.substr(0, colon_index);
		string content = line.substr(colon_index + 2, line.length() - colon_index + 2);

		if (content.empty())
			continue;

		if (qualifier == "folder_name")
		{
			folder_name = content.c_str();
			cout << "  folder_name: " << folder_name << endl;
		}
		else if (qualifier == "model_name")
		{
			model_name = content.c_str();
			cout << "  model_name: " << model_name << endl;
		}
		else if (qualifier == "output_prefix")
		{
			output_prefix = content.c_str();
			cout << "  output_prefix: " << output_prefix << endl;
		}
	}

	//Close The File
	in.close();
	cout << endl << "DONE reading config settings. " << endl << endl;
}

int main(int argc, char* argv[]) {
	// if (argc < 3) {
	//  cout << "Usage: example_bin folder_name model_name" << endl;
	//  exit(1);
	//}
	//const string folder_name = argv[1];
	//const string model_name = argv[2];

	read_config();

	//stitch_mesh(folder_name, model_name);
	flatten_patches(folder_name, model_name);

	return 0;
}

void stitch_mesh(const string& folder_name, const string& model_name) {
	cout << "Reading model " << model_name << " from folder " << folder_name << endl;

	const string mesh_path = folder_name + string("/") + model_name + string("__concatenated.obj");
	const string seams_path = folder_name + string("/") + model_name + string("__concatenated_seams.txt");

	MatrixXd V; MatrixXi F, C; 
	igl::readOBJ(mesh_path, V, F); 
	igl::facet_components(F, C);
	cout << "Input has " << V.rows() << " vertices and " << F.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << endl;
	
	MatrixXd seams; 
	igl::readCSV(seams_path, seams);
	cout << "seams contains " << seams.rows() << " pairs" << endl;
	
	// snap vertices to each other
	for (auto i = 0; i < seams.rows(); i++) 
		V.row(seams(i, 0)) = V.row(seams(i, 1));

	MatrixXd newV; MatrixXi newVI, newVJ, newF;
	igl::remove_duplicate_vertices(V, F, 1e-7, newV, newVI, newVJ, newF);

	// count connected components of the stitched mesh
	igl::facet_components(newF, C);
	cout << "Output has " << newV.rows() << " vertices and " << newF.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << endl;
	igl::writeOBJ(folder_name + string("/") + model_name + string("__stitched.obj"), newV, newF);
}

void flatten_single_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& uv)
{
	Eigen::MatrixXi cut_F = F;
	Eigen::MatrixXd cut_V = V;

	vector<vector<int>> boundaries;
	igl::boundary_loop(F, boundaries);



	bool do_cut_mesh = false;
	//if(boundaries.size() == 2 && boundaries[0] > )

	//cut mesh
	if (boundaries.size() > 1)
	{
		cout << "boundaries[0].size() = " << boundaries[0].size() << endl;
		cout << "boundaries[1].size() = " << boundaries[1].size() << endl;



		vector<int> flat_boundary_indices;
		for (vector<int> loop : boundaries)
			flat_boundary_indices.insert(flat_boundary_indices.end(), loop.begin(), loop.end());


		std::vector<std::vector<int>> adjacency_VV;
		igl::adjacency_list(F, adjacency_VV);


		int point_index = boundaries[0][0];
		Eigen::RowVectorXd point = V.row(point_index);


		std::set<int> target_other_boundary(boundaries[1].begin(), boundaries[1].end());
		Eigen::VectorXd min_distance;
		Eigen::VectorXi previous;
		int target_index = igl::dijkstra(point_index, target_other_boundary, adjacency_VV, min_distance, previous);

		vector<int> cut_path_indices;
		igl::dijkstra(target_index, previous, cut_path_indices);
		cout << "shortest path at " << point_index << ": " << list_to_string(cut_path_indices) << endl;

		//dijkstra sometimes gives edge on boundary back, so remove that
		int n = cut_path_indices.size();
		if (is_contained(cut_path_indices[n - 1], flat_boundary_indices) && is_contained(cut_path_indices[n - 2], flat_boundary_indices))
			cut_path_indices.pop_back();


		Eigen::MatrixXd cut_path;
		index_to_value(cut_path_indices, V, cut_path);
		cut_path = cut_path.rowwise() + Eigen::RowVector3d(1.0, 0, 0);




		std::vector<std::vector<int>> VFi_unused;
		std::vector<std::vector<int>> adjacency_VF;
		igl::vertex_triangle_adjacency(V, F, adjacency_VF, VFi_unused);

		Eigen::MatrixXi adjacency_FF;
		igl::triangle_triangle_adjacency(F, adjacency_FF);
		//cout << "adjacency_FF: \n" << adjacency_FF  << endl;



		Eigen::MatrixXi cut_path_F(F.rows(), F.cols());
		cut_path_F.setZero();

		int n_V = V.rows();
		cut_V.resize(n_V + cut_path.rows(), 3);
		cut_V << V, cut_path;


		vector<int> faces_offset;

		int vi_current = cut_path_indices[0];
		int vi_next = cut_path_indices[1];
		//cout << endl << vi_current << " - " << vi_next << endl;

		vector<int> cutline_faces = get_edge_incident_faces(vi_current, vi_next, adjacency_VF);
		//cout << "  cutline_faces(" << vi_current << "): " << list_to_string(cutline_faces) << endl;

		int fi = get_oriented_incident_face_at(cutline_faces, vi_current, vi_next, F);
		//vector<int> contained_cut_indices;
		//int fi = get_oriented_incident_face_at(cutline_faces, cut_path_indices, F, contained_cut_indices);

		//visited: add all seen
		//stack: add viable ones
		std::vector<int> visited;
		//visited.push_back(unused_face);

		std::stack<int> stack;
		stack.push(fi);
		int f;

		while (!stack.empty())
		{
			// Pop a vertex from stack 
			f = stack.top();
			stack.pop();

			if (!is_contained(f, visited))
				visited.push_back(f);


			//check if face is valid 
			if (f < 0)
				continue;
			if (is_contained(f, faces_offset))
				continue;

			vector<int> contained_cut_indices;
			int oriented_fi = get_oriented_incident_face_at(f, cut_path_indices, F, contained_cut_indices);

			if (contained_cut_indices.size() < 1)
				continue;

			if (contained_cut_indices.size() == 2 && oriented_fi != f)
				continue;


			//add neighbors of valid face
			faces_offset.push_back(f);
			//cout << "faces_offset: " << list_to_string(faces_offset) << endl;

			auto adjacent_faces = adjacency_FF.row(f);
			//cout << "      adjacent_faces(" << f << "): " << adjacent_faces << endl;

			for (int j = 0; j < adjacent_faces.size(); j++)
				if (!is_contained(adjacent_faces(j), visited))
					stack.push(adjacent_faces(j));
		}


		//vector<int> face_current = adjacency_VF[vi_current];
		//vector<int> face_next = adjacency_VF[vi_next];
		//vector<int> edge_faces = intersect(face_current, face_next);
		//
		//cout << "  VV(" << vi_current << "):   " << list_to_string(adjacency_VV[vi_current]) << endl;
		//cout << "  VV(" << vi_next << "):   " << list_to_string(adjacency_VV[vi_next]) << endl;
		//cout << "  face_current: " << list_to_string(face_current) << endl;
		//cout << "  face_next:    " << list_to_string(face_next) << endl;
		//
		//cout << "  edge_faces:   " << list_to_string(edge_faces) << endl;
		//cout << "     F.row(" << edge_faces[0] << "):   " << F.row(edge_faces[0]) << endl;
		//cout << "     F.row(" << edge_faces[1] << "):   " << F.row(edge_faces[1]) << endl;
		////cout << "     edge_faces[0]:   " << F.row(edge_faces[0]) << endl;
		////cout << "     edge_faces[1]:   " << F.row(edge_faces[1]) << endl;

	/*
	for (int cut_index = 0; cut_index < cut_path_indices.size() - 1; cut_index++)
	{
		int vi_current = cut_path_indices[cut_index];
		int vi_next = cut_path_indices[cut_index+1];
		cout << endl << vi_current << " - " << vi_next << endl;

		vector<int> cutline_faces = get_edge_incident_faces(vi_current, vi_next, adjacency_VF);
		cout << "  cutline_faces(" << vi_current << "): " << list_to_string(cutline_faces) << endl;

		if (cutline_faces.size() < 2)
		{
			cout << "\nERROR : cutline_faces.size() != 2\n" << endl;
			faces_offset.push_back(f);

			break;
		}

		bool has_found_face = false;

		for (int i = 0; i < 2; i++)
		{
			if (has_found_face)
				continue;

			int fi = cutline_faces[i];
			cout << "    F.row(" << fi << "): " << F.row(fi) << endl;

			int vi1 = index_of(vi_current, F.row(fi));
			int vi2 = (vi1 + 1) % 3;

			if (F(fi, vi2) != vi_next)
				continue;

			has_found_face = true;

			if (!is_contained(fi, faces_offset))
				faces_offset.push_back(fi);


			int unused_face = cutline_faces[(i + 1) % 2];
			auto adjacent_faces = adjacency_FF.row(fi);
			cout << "      unused_face(" << fi << "): " << unused_face << endl;
			cout << "      adjacent_faces(" << fi << "): " << adjacent_faces << endl;

			for (int j = 0; j < adjacent_faces.size(); j++)
			{
				int f = adjacent_faces(j);
				if (f == unused_face || f < 0)
					continue;

				//store in list for debug
				if (!is_contained(f, faces_offset))
					faces_offset.push_back(f);
			}

		}

		cout << "faces_offset: " << list_to_string(faces_offset) << endl;

	}
	*/

	//cout << endl << endl << "adapt vertex indices for boundary, all that lie on cut line" << endl << endl;
	//adapt vertex indices for boundary, all that lie on cut line
		for (int i = 0; i < faces_offset.size(); i++)
		{
			int fi = faces_offset[i];
			//cout << "        F.row(" << fi << "): " << F.row(fi) << endl;

			for (int k = 0; k < F.row(fi).size(); k++)
			{
				int vi_f = F(fi, k);
				int ci_f = index_of(vi_f, cut_path_indices);
				if (ci_f >= 0)
					cut_F(fi, k) = ci_f + n_V;
			}
			//cout << "    cut_F.row(" << fi << "): " << cut_F.row(fi) << endl;

			if (cut_F.maxCoeff() >= cut_V.rows())
				cout << "\nERROR : cut_F.maxCoeff() >= cut_V.rows()\n" << endl;
		}
	}

	vector<vector<int>> boundaries2;
	igl::boundary_loop(cut_F, boundaries2);
	if (boundaries2.size() > 1)
		cout << "\nERROR : boundaries2.size() > 1\n" << endl;

	Eigen::VectorXi bnd;
	Eigen::MatrixXd bnd_uv;
	igl::boundary_loop(cut_F, bnd);


	//hacking around errors
	uv = cut_V;
	V = cut_V;
	F = cut_F;
	if (cut_F.rows() == 1)
	{
		cout << "--> only 1 triangle; return directly" << endl;
		return;
	}
	if (cut_V.rows() == bnd.rows())
	{
		cout << "--> everything is on boundary, no inner vertices; return directly" << endl;
		return;
	}

	igl::map_vertices_to_circle(cut_V, bnd, bnd_uv);

	igl::harmonic(cut_V, cut_F, bnd, bnd_uv, 1, uv);
	if (igl::flipped_triangles(uv, cut_F).size() != 0) {
		igl::harmonic(cut_F, bnd, bnd_uv, 1, uv); // use uniform laplacian
	}

	igl::SLIMData sData;
	sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

	Eigen::VectorXi b;
	Eigen::MatrixXd bc;
	igl::slim_precompute(cut_V, cut_F, uv, sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);
	igl::slim_solve(sData, 20); // 20 iters
	uv = sData.V_o;
	
	V = cut_V;
	F = cut_F;
	cout << "energy after flattening: " << sData.energy << endl;


	//Eigen::VectorXd lengths_3D;
	//Eigen::VectorXd lengths_2D;
	//igl::edge_lengths(cut_V, cut_F, lengths_3D);
	//igl::edge_lengths(uv, cut_F, lengths_2D);
	//double avg_edge_length = igl::avg_edge_length(cut_V, cut_F);

	//Eigen::VectorXd lengths_difference = (lengths_3D - lengths_2D).cwiseAbs();
	//double max_difference = lengths_difference.maxCoeff();
	//
	//
	////cout << "lengths_3D: \n" << lengths_3D << endl;
	////cout << "lengths_2D: \n" << lengths_2D << endl;
	////cout << "lengths_difference: \n" << lengths_difference << endl;

	//cout << "max_difference:          " << max_difference << " (" << max_difference/avg_edge_length << "%)" << endl;

}

void flatten_patches(const string& folder_name, const string& model_name)
{
	const string mesh_path = folder_name + string("/") + model_name; // +string("__concatenated.obj");
	MatrixXd V_concatenated;
	MatrixXi F_concatenated;
	igl::readOBJ(mesh_path, V_concatenated, F_concatenated);

	const string output_path = folder_name + string("/") + output_prefix;

	std::vector<std::vector<int>> VFi_unused;
	std::vector<std::vector<int>> adjacency_VF;
	igl::vertex_triangle_adjacency(V_concatenated, F_concatenated, adjacency_VF, VFi_unused);

	 MatrixXi C; 
	 igl::facet_components(F_concatenated, C);
	 ofstream fileCF("facet_components.txt");
	fileCF << C << endl;
	 fileCF.close();

	//Eigen::VectorXi CV;
	//igl::vertex_components(F_concatenated, CV);
	//ofstream fileCV("vertex_components.txt");
	//	fileCV << CV << endl;
	//fileCV.close();


	int number_parts = C.maxCoeff() + 1;
	std::vector<std::vector<int>> face_components(number_parts);
	std::vector<std::vector<int>> vertex_components(number_parts);

	for (int fi = 0; fi < C.rows(); fi++)
	{
		int index = C(fi);
		face_components[index].push_back(fi);

		Eigen::RowVectorXi v = F_concatenated.row(fi);
		vector<int> vertices(v.data(), v.data() + v.size());
		vertex_components[index].insert(vertex_components[index].end(), vertices.begin(), vertices.end());
	}


	for (int i = 0; i < number_parts; i++)
	{
		cout << endl << endl << endl << "FLATTEN PART " << i << endl << endl;

		//Mesh part_mesh;
		//part_mesh.V = V_concatenated;
		//index_to_value(face_components[i], F_concatenated, part_mesh.F);

		Mesh part_mesh = sub_mesh(vertex_components[i], face_components[i], V_concatenated, F_concatenated);
		if (part_mesh.F.rows() < 1)
		{
			//part_start = part_end;
			cout << "\nERROR : part_mesh.F.rows() < 1, skipping this part\n" << endl;

			continue;
		}


		igl::writeOBJ(output_path + string("part") + to_string(i) + ".obj", part_mesh.V, part_mesh.F);

		Eigen::MatrixXd subUV;
		flatten_single_mesh(part_mesh.V, part_mesh.F, subUV);
		igl::writeOBJ(output_path + string("part-cut") + to_string(i) + ".obj", part_mesh.V, part_mesh.F);

		// To save it as an OBJ we need to have 3 columns, i.e. x,y,z coordinates, so we will just save the 'z' as 0
		Eigen::MatrixXd subUV3d(subUV.rows(), 3);
		subUV3d.setZero();
		subUV3d.col(0) = subUV.col(0);
		subUV3d.col(1) = subUV.col(1);

		igl::writeOBJ(output_path + string("flat") + to_string(i) + ".obj", subUV3d, part_mesh.F);
	}
}