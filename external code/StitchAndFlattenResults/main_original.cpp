#include <igl/opengl/glfw/Viewer.h>
#include <igl/readCSV.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/facet_components.h>
#include <igl/remove_duplicate_vertices.h>

#include <igl/slim.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/MappingEnergyType.h>
#include <igl/flipped_triangles.h>

using namespace std;
using namespace Eigen;
using namespace igl;

void stitch_mesh(const string& folder_name, const string& model_name);
void flatten_patches(const string& folder_name, const string& model_name);

int main(int argc, char *argv[]) {
   if (argc < 3) {
    cout << "Usage: example_bin folder_name model_name" << endl;
    exit(1);
  }
  const string folder_name = argv[1];
  const string model_name = argv[2];

  //stitch_mesh(folder_name, model_name);
  flatten_patches(folder_name, model_name);
  
  return 0;
}

void stitch_mesh(const string& folder_name, const string& model_name) {
  cout << "Reading model " << model_name << " from folder " << folder_name << endl;

  const string mesh_path = folder_name + string("/") + model_name + string("__concatenated.obj");
  const string seams_path = folder_name + string("/") + model_name + string("__concatenated_seams.txt");

  MatrixXd V; MatrixXi F,C; readOBJ(mesh_path, V, F); facet_components(F,C);
  cout << "Input has " << V.rows() << " vertices and " << F.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << endl;
  MatrixXd seams; readCSV(seams_path, seams);
  cout << "seams contains " << seams.rows() << " pairs" << endl; 
  // snap vertices to each other
  for (auto i = 0; i < seams.rows(); i++) V.row(seams(i,0)) = V.row(seams(i,1));
  MatrixXd newV; MatrixXi newVI,newVJ,newF;
  remove_duplicate_vertices(V,F,1e-7,newV,newVI,newVJ,newF);

  // count connected components of the stitched mesh
  facet_components(newF,C);
  cout << "Output has " << newV.rows() << " vertices and " << newF.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << endl;
  writeOBJ(folder_name +string("/") + model_name + string("__stitched.obj"),newV,newF);
}

void flatten_single_mesh(const MatrixXd& V, const MatrixXi& F, MatrixXd& uv) {
  VectorXi bnd; MatrixXd bnd_uv;
  boundary_loop(F,bnd);
  map_vertices_to_circle(V,bnd,bnd_uv);

  harmonic(V,F,bnd,bnd_uv,1,uv);
  if (flipped_triangles(uv,F).size() != 0) {
    harmonic(F,bnd,bnd_uv,1,uv); // use uniform laplacian
  }

  SLIMData sData;
  sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
  Eigen::VectorXi b; Eigen::MatrixXd bc;
  slim_precompute(V,F,uv,sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b,bc,0);
  slim_solve(sData,20); // 20 iters
  uv = sData.V_o;
}

void flatten_patches(const string& folder_name, const string& model_name) {
  const string mesh_path = folder_name + string("/") + model_name + string("__concatenated.obj");
  MatrixXd V; 
  MatrixXi F,C; 
  igl::readOBJ(mesh_path, V, F); 
  facet_components(F,C);

  int comp_num = C.maxCoeff() + 1;
  for (int comp_i = 1; comp_i <= comp_num; comp_i++) {
    MatrixXd subV, subUV; MatrixXi subF; 
    string submesh_path = folder_name + string("/") + model_name + string("__parts") + to_string(comp_i) + ".obj";
    string submesh_uv_path = folder_name + string("/") + model_name + string("__parts") + to_string(comp_i) + "-uv.obj";
    readOBJ(submesh_path, subV, subF);


    flatten_single_mesh(subV, subF, subUV);

    // To save it as an OBJ we need to have 3 columns, i.e. x,y,z coordinates, so we will just save the 'z' as 0
    MatrixXd subUV3d(subUV.rows(),3);
    subUV3d.setZero(); subUV3d.col(0) = subUV.col(0); subUV3d.col(1) = subUV.col(1);
    writeOBJ(submesh_uv_path, subUV3d, subF);
  }
}