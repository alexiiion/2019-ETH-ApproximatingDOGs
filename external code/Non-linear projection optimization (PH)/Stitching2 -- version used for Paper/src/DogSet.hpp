#pragma once

#include <Eigen/Dense>
#include <vector>
#include <igl/AABB.h>


class DogSet
{
    int n; //number of DOGs
    
    std::vector<Eigen::MatrixXi> F;
    std::vector<Eigen::MatrixXd> V;
    std::vector<igl::AABB<Eigen::MatrixXd, 3>> trees;
    
    
    void init();
    
public:
    
    //concatentated DOGs
    Eigen::MatrixXd Vall;
    Eigen::MatrixXi Fall;
    
    Eigen::VectorXi I;
    std::vector<Eigen::MatrixXd> projected;
        
    DogSet(std::string fname);
    DogSet(const Eigen::MatrixXd& V0, const Eigen::MatrixXi& F0);
    
	static double weight_normals;
	static bool use_distance_smoothness;
    
    void projectPoint(Eigen::RowVector3d& p, const int l);
    
    void projectPoints(const Eigen::MatrixXd& V, const Eigen::VectorXi& L, Eigen::MatrixXd& VP) const;
    
    void computeDistances(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin);
    void computeDistancesAlongMeshNormals(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin);
    void computeNormalDeviations(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin);

    void getPointsFromLabels(const Eigen::VectorXi& L, Eigen::MatrixXd& Vout);
    
    void minimalAssignement(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin);//, Eigen::MatrixXd& Vout);
    
    void graphCutAssignement(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin, const double smoothnessScale = 2.0);//, Eigen::MatrixXd& Vout);

	//void set_weight_normals(double& weight);
};