#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <igl/AtA_cached.h>
#include "DogSet.hpp"

class OptimizationObjectives
{
    // precompute Jacobian for smoothness
    Eigen::SparseMatrix<double> Jsmooth, Jsmooth2;
    
    // precompute Jacobian for seam smoothness
    Eigen::SparseMatrix<double> JseamSmooth, JseamSmooth2;
    
    // precompute Jacobian for pairs
    Eigen::SparseMatrix<double> Jpairs, Jpairs2;
    
    // precompute Jacobian for boundary position
    Eigen::SparseMatrix<double> Jbpos, Jbpos2;
    
    // precompute Jacobian for position
    Eigen::SparseMatrix<double> Jpos, Jpos2;
    
    // precompute data structures for precomputed sparse matrix construction
    Eigen::SparseMatrix<double> Jdev, Jdev2;
    
    // Identity
    Eigen::SparseMatrix<double> Id;
    
    Eigen::VectorXi devJPrecompute;
    igl::AtA_cached_data devHCache;
    
    // solver data
    int iter = 0;

    
    bool choleskyInitialized = false;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
    
    // a position each vertex should stay close to
    Eigen::MatrixXd posConstraints;
    
    // position for the boundary. Vertices should stay close to initial projection.
    Eigen::MatrixXd boundaryPosConstraints;
    
    // mesh triangles
    const Eigen::MatrixXi& F;
    const Eigen::VectorXi& L;
    const Eigen::VectorXi& I;
    
    // projector
    const DogSet& projector;
    
    // flag indicating boundary
    std::vector<char> boundary;
    std::vector<char> merged;
    
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> areaNormalization;
    std::vector<double> invSqrtArea;
    
public:
    double wSmooth = 10;
    double wSeamProx = 100000;
    double wSeamSmooth =  1;
    double wProx = 5;
    double wPair = 1;
    double wDev = 10000;
    
    OptimizationObjectives(const Eigen::MatrixXd& V0, const Eigen::VectorXi& L, const Eigen::MatrixXi& F, const Eigen::VectorXi& I, const DogSet& project, const std::vector<std::array<int, 2>>& pairs);
    
    // current hessian
    Eigen::SparseMatrix<double> H;
    
    // current grad
    Eigen::VectorXd grad;
        
    void checkCurvatureDerivative(const Eigen::MatrixXd& V);
    
    void angleDefect(const Eigen::MatrixXd& V, Eigen::VectorXd& vals, const bool excludeBoundary = true);
    
    // add objectives
    double gaussianCurvatureObjective(const Eigen::MatrixXd& V, const double w = .0, bool sqrtAreaNormalization = true);
    
    double smoothnessObjective(const Eigen::MatrixXd& V, double w = .0);
    
    double seamPairObjective(const Eigen::MatrixXd& V, double w = .0);

    double proximityObjective(const Eigen::MatrixXd& V, double w = .0);
    
    double seamProximityObjective(const Eigen::MatrixXd& V, double w = .0);
    
    double seamSmoothnessObjective(const Eigen::MatrixXd& V, double w = .0);
    
    double objective(const Eigen::MatrixXd& V);
    
    double setHessianAndGradient(const Eigen::MatrixXd& V);
    
    void addLevenbergMarquardtCorrection(const double t);
    
    void mergeParts(const Eigen::MatrixXd& V, const double thres = 0.01);
    
    int optimize(Eigen::MatrixXd& V, const double t0, const double lm);
};
