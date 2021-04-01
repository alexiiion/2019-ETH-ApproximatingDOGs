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
    
    void splitMeshes(const double smoothFactor = 0.);

    void seamEdges(Eigen::MatrixXi& E);
    
    void smoothSeams(const Eigen::MatrixXi& E);
    
    std::vector<std::array<int, 2>> getPairs();    
};
