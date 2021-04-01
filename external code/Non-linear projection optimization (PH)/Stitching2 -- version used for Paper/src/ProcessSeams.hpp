#pragma once
#include <Eigen/Dense>
#include <vector>
#include <array>

class ProcessSeams
{
public:
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
   
    // part labels
    Eigen::VectorXi L;
    
    // original vertex index
    Eigen::VectorXi I;
    

    ProcessSeams(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_, const Eigen::VectorXi& L_);
    
    void smoothLabels();
    
    void seamEdges(Eigen::MatrixXi& E);

    void smoothSeams(const Eigen::MatrixXi& E);
    
    void splitMeshes(const double smoothFactor = 0.);
    void stitchMesh(const Eigen::MatrixXd& current_V,  Eigen::MatrixXd& out_V, Eigen::MatrixXi& out_F);

    std::vector<std::array<int, 2>> getPairs();
    std::vector<std::vector<int>> getSeamVertices();

};
