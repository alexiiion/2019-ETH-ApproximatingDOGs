#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include "objectives.hpp"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/cat.h>
#include <igl/sparse_cached.h>
#include <unsupported/Eigen/SparseExtra>
#include <igl/slice_into.h>
#include <igl/boundary_loop.h>
#include <igl/writeOFF.h>
#include <igl/adjacency_list.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/adjacency_matrix.h>

#include <igl/invert_diag.h>

//Eigen::VectorXd V0_mass;

OptimizationObjectives::OptimizationObjectives(const Eigen::MatrixXd& V0,
                                               const Eigen::VectorXi& L_,
                                               const Eigen::MatrixXi& F_,
                                               const Eigen::VectorXi& I_,
                                               const DogSet& projector_,
                                               const std::vector<std::array<int, 2>>& pairs)
: L(L_), F(F_), projector(projector_), I(I_)
{
    const int nv = V0.rows();
    grad.resize(3 * nv);
    grad.setZero();
    H.resize(3 * nv, 3 * nv);
    
    projector.projectPoints(V0, L, posConstraints);
   // posConstraints = VP;
    
    // save boundary flag
    boundary.resize(nv, 0);
    merged.resize(nv, -1);
    
    std::vector<std::vector<int>> loops;
    igl::boundary_loop(F, loops);
    for(auto& l : loops) 
        for(int i : l) 
            boundary[i] = 1;
    
    for(auto p : pairs)
    {
        boundary[p[0]] = 2;
        boundary[p[1]] = 2;
    }
    
    
    // boundary proximity
    int nb = 0;
    for(int i : boundary) 
        if(i == 1) 
            ++nb;
    boundaryPosConstraints.resize(nb, 3);
    
    int cnt = 0;
    std::vector<Eigen::Triplet<double>> trip;
    
    for(int i = 0; i < nv; ++i)
        if(boundary[i] == 1)
        {
            boundaryPosConstraints.row(cnt) = posConstraints.row(i);
           
            for(int k = 0; k < 3; ++k)
                trip.emplace_back(3 * cnt + k, i + k * nv, 1.);
            
            ++cnt;
        }
    
    Jbpos.resize(3 * cnt, 3 * nv);
    Jbpos.setFromTriplets(trip.begin(), trip.end());
    trip.clear();
    
    Jbpos2 = Jbpos.transpose() * Jbpos;
    

    //compute area weight for gaussian curvature normalization
    updateGaussianCurvatureNormalization(V0);


    
    //compute area weights -- still used???
    invSqrtArea.resize(nv, 0);
    for(int i = 0; i < F.rows(); ++i)
    {
        const Eigen::Vector3d e0 = V0.row(F(i, 1)) - V0.row(F(i, 0));
        const Eigen::Vector3d e1 = V0.row(F(i, 2)) - V0.row(F(i, 0));
        const double ai = 0.5 * e0.cross(e1).norm();
        
        for(int j = 0; j < 3; ++j)
            invSqrtArea[F(i, j)] += ai;
    }
    
    for(auto& x : invSqrtArea) 
        x = sqrt(x / 3.);
       

    areaNormalization.resize(3 * nv);
    for(int k = 0; k < 3; ++k)
        areaNormalization.diagonal().block(k * nv, 0, nv, 1) = Eigen::Map<Eigen::VectorXd>(invSqrtArea.data(), nv);




    
    // initialize smoothness Jacobian
    Eigen::SparseMatrix<double> L, M, A;
   // igl::cotmatrix(V0, F, L);
 
    
    {
        igl::adjacency_matrix(F, L);
        // sum each row
        Eigen::SparseVector<double> Asum;
        igl::sum(L,1,Asum);
        // Convert row sums into diagonal of sparse matrix
        Eigen::SparseMatrix<double> Adiag;
        igl::diag(Asum,Adiag);
        // Build uniform laplacian
        Eigen::SparseMatrix<double> U;
        L = Adiag - L;
    }
    
    igl::massmatrix(V0, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    A = M.cwiseSqrt() * L;

    Jsmooth.resize(3 * nv, 3 * nv);
    
    Eigen::VectorXi R(nv);
    for(int i = 0; i < nv; ++i) 
        R(i) = i;
    
    igl::slice_into(A, R, R, Jsmooth);
    
    R.array() += nv;
    igl::slice_into(A, R, R, Jsmooth);
      
    R.array() += nv;
    igl::slice_into(A, R, R, Jsmooth);
    
    Jsmooth = areaNormalization * Jsmooth;
    Jsmooth2 = Jsmooth.transpose() * Jsmooth;
    

    // initialize proximity Jacobian
    Jpos.resize(3 * nv, 3 * nv);
    Jpos = areaNormalization;
    Jpos2 = Jpos.transpose() * Jpos;
   
    
    // initialize seam Jacobian
    trip.clear();
    cnt = 0;
    
    for(const auto& p : pairs)
    {
        const double w = std::pow(std::pow(invSqrtArea[p[0]], -2) + std::pow(invSqrtArea[p[1]], -2), -0.5);
        
        for(int k = 0; k < 3; ++k)
        {
            trip.emplace_back(cnt, p[0] + k * nv, -1. / w);
            trip.emplace_back(cnt, p[1] + k * nv,  1. / w);
                      
            ++cnt;
        }
    }
    
    Jpairs.resize(3 * pairs.size(), 3 * nv);
    Jpairs.setFromTriplets(trip.begin(), trip.end());
    Jpairs2 = Jpairs.transpose() * Jpairs;
    
    // boundary smoothness
    
    std::vector<std::vector<int>> adj;
    igl::adjacency_list(F, adj);
    
    for(auto p : pairs)
    {
        auto merged = adj[p[0]];
        std::copy(adj[p[1]].begin(), adj[p[1]].end(), std::back_inserter(merged));
        std::sort(merged.begin(), merged.end());
        merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
        adj[p[0]] = adj[p[1]] = merged;
    }
    


    std::vector<char> boundaryFlag(nv, 0);
    
    for(const auto& l : loops)
    {
        for(auto i : l)
        {
            boundaryFlag[i] = 1;
        }
    }
    
    for(auto p : pairs)
    {
        boundaryFlag[p[0]] = 2;
        boundaryFlag[p[1]] = 2;
    
        for(int j : adj[p[0]])
        {
            if(boundaryFlag[j] == 1) boundaryFlag[p[0]] = 3;
        }
        
        for(int j : adj[p[1]])
        {
            if(boundaryFlag[j] == 1) boundaryFlag[p[1]] = 3;
        }
    }
    
    trip.clear();
    
    for(int i = 0; i < nv; ++i)
    {
        const char flag = boundaryFlag[i];
        
        if(flag)
        {
            std::vector<int> nbh;
            for(int j : adj[i])
                if(boundaryFlag[j] || (flag == 2 && boundaryFlag[j] == 2))
                    nbh.push_back(j);
            
            if(!nbh.empty())
            {
                const double w = -1. / nbh.size();
                
                for(int k = 0; k < 3; ++k)
                {
                    for(int j : nbh)
                        trip.emplace_back(nv * k + i, nv * k + j, w);
                    
                    trip.emplace_back(nv * k + i, nv * k + i, 1.);
                }
            }
        }
    }
        
    JseamSmooth.resize(3 * nv, 3 * nv);
    JseamSmooth.setFromTriplets(trip.begin(), trip.end());
    JseamSmooth2 = JseamSmooth.transpose() * JseamSmooth;
    
    Id.resize(3 * nv, 3 * nv);
    Id.setIdentity();
}


void OptimizationObjectives::updateGaussianCurvatureNormalization(const Eigen::MatrixXd& V)
{
    Eigen::SparseMatrix<double> M, Minv;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);

    gaussianCurvatureNormalization.resize(V.rows());
    gaussianCurvatureNormalization = Minv.diagonal();
}

void OptimizationObjectives::mergeParts(const Eigen::MatrixXd& V, const double thres)
{
    const int nv = V.rows();
    
    Eigen::VectorXd ad;
    angleDefect(V, ad, false);
    
    std::map<int, double> adb;
    
    for(int i = 0; i < nv; ++i)
    {
        if(boundary[i] == 2)
        {
            adb[I(i)] += ad(i);
        }
    }
    
    std::vector<char> flag(nv, 0);
    for(auto x : adb)
    {
        if(std::abs(x.second - 2 * M_PI) < thres)
            flag[x.first] = 1;
    }
    
    
    std::fill_n(merged.begin(), nv, -1);
    
    for(int i = 0; i < nv; ++i)
        if(flag[I(i)]) 
            merged[i] = I(i);
}

void OptimizationObjectives::angleDefect(const Eigen::MatrixXd& V, Eigen::VectorXd& vals, const bool excludeBoundary)
{
     const int nf = F.rows();
     const int nv = V.rows();
    
    vals.resize(nv);

    for(int i = 0; i < nv; ++i)
    {
        if(boundary[i])
        {
            vals(i) = .0;
        } else vals[i] = -2.0 * M_PI;
    }
    
    for(int i = 0; i < nf; ++i)
    {
        const Eigen::Vector3d t[3]{V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2))};
        
        for(int j = 0; j < 3; ++j)
        {
            Eigen::Vector3d u = t[(j+1)%3] - t[j];
            Eigen::Vector3d v = t[(j+2)%3] - t[j];
            
            const double u2 = u.dot(u);
            const double v2 = v.dot(v);
            
            if(!excludeBoundary || boundary[F(i,j)] == 0)
            {
                double val = acos(std::max(-1., std::min(1., u.dot(v) / sqrt(u2 * v2))));
                if(std::isnan(val)) 
                    std::cout << "NAN value\n";

                vals(F(i, j)) += val;
            }
        }
    }
}

void OptimizationObjectives::addLevenbergMarquardtCorrection(const double t)
{
    if(t) H += t * Id;
}

double OptimizationObjectives::gaussianCurvatureObjective(const Eigen::MatrixXd& V, const double w, bool sqrtAreaNormalization)
{
    
    const int nf = F.rows();
    const int nv = V.rows();
   
    Eigen::VectorXd vals(nv);
    
    for(int i = 0; i < nv; ++i)
    {
        if(boundary[i]) 
            vals(i) = 0;
        else 
            vals[i] = -2.0 * M_PI;
    }
        
    std::vector<Eigen::Triplet<double>> trip;
    Eigen::Vector3d n;
        
    for(int i = 0; i < nf; ++i)
    {
        const Eigen::Vector3d t[3]{V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2))};
        if(w) n = (t[1] - t[0]).cross(t[2] - t[0]).normalized();
    
        for(int j = 0; j < 3; ++j)
        {
            if(boundary[F(i,j)] == 0)
            {
                Eigen::Vector3d u = t[(j+1)%3] - t[j];
                Eigen::Vector3d v = t[(j+2)%3] - t[j];
                
                const double u2 = u.dot(u);
                const double v2 = v.dot(v);
                
              //  bool degenerate = u2 * v2 < 1e-12 || n != n;
                bool degenerate = false;
                
                if(!degenerate)
                {
                    double val = acos(std::max(-1., std::min(1., u.dot(v) / sqrt(u2 * v2))));
                    if(std::isnan(val)) std::cout << "NAN value\n";
                    
                    vals(F(i, j)) += val;
                }
                        
                if(w)
                {
                    Eigen::Matrix3d d;
                    
                    if(!degenerate)
                    {
                        d.col(1) = n.cross(u) / -u2;
                        d.col(2) = n.cross(v) / v2;
                        d.col(0) = -(d.col(1) + d.col(2));
                    } else d.setZero();
                    
                    const int ids[3]{F(i, j), F(i, (j+1)%3), F(i, (j+2)%3)};
                    
                    for(int l = 0; l < 3; ++l)
                    {
                        for(int k = 0; k < 3; ++k)
                        {
                            trip.emplace_back(ids[0], k * nv + ids[l], d(k, l));
                        }
                    }
                }
            }
        }
    }
 
    if(sqrtAreaNormalization)
    {
        for(int i = 0; i < nv; ++i)
            vals(i) *= gaussianCurvatureNormalization(i);

        if(w)
        {
            for(auto& t : trip)
            {
                t = Eigen::Triplet<double>(t.row(), t.col(), t.value() * gaussianCurvatureNormalization(t.row()));
            }
        }
    }
    
    if(w)
    {
        if(0)
        {
            if(!Jdev.size())
            {
                Jdev.resize(nv, 3 * nv);
                igl::sparse_cached_precompute(trip, devJPrecompute, Jdev);
            } else
            {
                igl::sparse_cached(trip, devJPrecompute, Jdev);
            }
            
            if(!Jdev2.size())
            {
                igl::AtA_cached_precompute(Jdev, devHCache, Jdev2);
            } else
            {
                igl::AtA_cached(Jdev, devHCache, Jdev2);
            }
        } else
        {
            Jdev.resize(nv, 3 * nv);
            Jdev.setFromTriplets(trip.begin(), trip.end());
            Jdev2 = Jdev.transpose() * Jdev;
        }
        
        H += 0.5 * w * Jdev2;
        grad += w * (Jdev.transpose() * vals);
    }
    

    //std::cout << "    gaussianCurvatureObjective: |max| =  " << vals.cwiseAbs().maxCoeff() << std::endl;
    currentDevelopabilityObjective.resize(vals.rows());
    currentDevelopabilityObjective = vals;
    return vals.squaredNorm();
}

double OptimizationObjectives::seamPairObjective(const Eigen::MatrixXd& V, const double w)
{
    currentSeamProximityObjective.resize(V.rows());

    if (w)
    {
        const Eigen::VectorXd vals = Jpairs * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size());

        grad += w * Jpairs.transpose() * vals;
        H += 0.5 * w * Jpairs2;

        currentSeamProximityObjective = vals;
        return vals.squaredNorm();
    }

    return (Jpairs * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size())).squaredNorm();
}

double OptimizationObjectives::seamProximityObjective(const Eigen::MatrixXd& V, double w)
{
    currentBoundaryProximityObjective.resize(V.rows());

    const int n = boundaryPosConstraints.size();
    Eigen::VectorXd vals(n);
    
    int cnt = 0;
    
    for(int i = 0; i < boundary.size(); ++i)
        if(boundary[i] == 1)
        {
            for(int k = 0; k < 3; ++k)
            {
                vals(3 * cnt + k) = V(i, k) - boundaryPosConstraints(cnt, k);
            }
            
            ++cnt;
        }
    
    if(w)
    {
        grad += w * Jbpos.transpose() * vals;
        H += 0.5 * w * Jbpos2;
    }
    
    currentBoundaryProximityObjective = vals;
    return vals.squaredNorm();
}

double OptimizationObjectives::proximityObjective(const Eigen::MatrixXd& V, const double w)
{
    currentProximityObjective.resize(V.rows());

    Eigen::VectorXd vals = areaNormalization * (Eigen::Map<const Eigen::VectorXd>(V.data(), V.size()) -  Eigen::Map<const Eigen::VectorXd>(posConstraints.data(), posConstraints.size()));
    
    if(w)
    {
        grad += w * Jpos.transpose() * vals;
        H += 0.5 * w * Jpos;
    }
    
    currentProximityObjective = vals;
    return vals.squaredNorm();
}

double OptimizationObjectives::seamSmoothnessObjective(const Eigen::MatrixXd& V, const double w)
{
    currentSeamSmoothnessObjective.resize(V.rows());

    if(w)
    {
        const Eigen::VectorXd vals = JseamSmooth * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size());
        grad += w * JseamSmooth.transpose() * vals;
        H += 0.5 * w * JseamSmooth2;
        
        currentSeamSmoothnessObjective = vals;
        return vals.squaredNorm();
    }
    
    return (JseamSmooth * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size())).squaredNorm();
}

double OptimizationObjectives::smoothnessObjective(const Eigen::MatrixXd& V, const double w)
{
    currentSmoothnessObjective.resize(V.rows());

    if(w)
    {
        const Eigen::VectorXd vals = Jsmooth * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size());
        grad += w * Jsmooth.transpose() * vals;
        H += 0.5 * w * Jsmooth2;
        
        currentSmoothnessObjective = vals;
        return vals.squaredNorm();
    }
    
    return (Jsmooth * Eigen::Map<const Eigen::VectorXd>(V.data(), V.size())).squaredNorm();
}


double OptimizationObjectives::objective(const Eigen::MatrixXd& V)
{
    double oDev = gaussianCurvatureObjective(V);
    double oPair = seamPairObjective(V);
    double oProx = proximityObjective(V);
    double oSeamProx = seamProximityObjective(V);
    double oSmooth = smoothnessObjective(V);
    double oSeamSmooth = seamSmoothnessObjective(V);
    
    return wDev * oDev + wPair * oPair + wProx * oProx + wSmooth * oSmooth + wSeamSmooth * oSeamSmooth + wSeamProx * oSeamProx;
}

double OptimizationObjectives::setHessianAndGradient(const Eigen::MatrixXd& V)
{
    H.setZero();
    grad.setZero();
    
    double oDev = gaussianCurvatureObjective(V, wDev);
    double oPair = seamPairObjective(V, wPair);
    double oProx = proximityObjective(V, wProx);
    double oSeamProx = seamProximityObjective(V, wSeamProx);
    double oSmooth = smoothnessObjective(V, wSmooth);
    double oSeamSmooth = seamSmoothnessObjective(V, wSeamSmooth);
    
    return wDev * oDev + wPair * oPair + wProx * oProx + wSmooth * oSmooth + wSeamSmooth * oSeamSmooth + wSeamProx * oSeamProx;

}

void OptimizationObjectives::checkCurvatureDerivative(const Eigen::MatrixXd& V)
{
    const int nv = V.rows();
    Eigen::MatrixXd V0 = V;
    
    double eps = 1e-5;
    double o0 = gaussianCurvatureObjective(V0);
    
    Eigen::VectorXd fd(3 * nv);
    
    for(int j = 0; j < 3; ++j)
    {
        
        for(int i = 0; i < nv; ++i)
        {
            V0(i, j) += eps;
            
            
            double o1 = gaussianCurvatureObjective(V0);
            
            fd(j * nv + i ) = (o1 - o0) / eps;
            
            V0(i, j) -= eps;
        }
    }
    
    grad.setZero();
    gaussianCurvatureObjective(V, 1.);
    
    
    {
        std::ofstream file("../d");
        file << std::setprecision(64);
        file << grad;
        file.close();
    }
    
    {
        std::ofstream file("../fd");
        file << std::setprecision(64);
        file << fd;
        file.close();
    }
}

int OptimizationObjectives::optimize(Eigen::MatrixXd& V, const double t0, const double levMaBeta)
{
    //update mass matrix for area normalization of gaussian curvature
    updateGaussianCurvatureNormalization(V);


    // build gradient and hession
    
 //   posConstraints = V;
    projector.projectPoints(V, L, posConstraints);
    double obj = setHessianAndGradient(V);
    
    //Eigen::VectorXd ad;
   // angleDefect(V, ad);
   // obj = ad.squaredNorm();
    
    addLevenbergMarquardtCorrection(levMaBeta);
    
    // solve system
    if(!choleskyInitialized)
    {
        chol.analyzePattern(H);
        choleskyInitialized = true;
    }
    
    chol.factorize(H);
    if(chol.info()) std::cout << "Cholesky factorization did not succeed : " << chol.info() << std::endl;
    Eigen::VectorXd dir = chol.solve(grad);
    
    // line search
    double step = t0;
    const int nv = V.rows();

    Eigen::MatrixXd Vi = V - step * Eigen::Map<Eigen::MatrixXd>(dir.data(), nv, 3);
    
    double obji;
    int cnt = 24;
    
    while(1)
    {
        if(cnt == 0)
        {
            V.swap(Vi);
            break;
        }
        
        obji = objective(Vi);
     
       // angleDefect(V, ad);
        //obji = ad.squaredNorm();
          
        
        if(obji < obj) break;
        
        Vi += (step - step / 2.) * Eigen::Map<Eigen::MatrixXd>(dir.data(), nv, 3);
        step /= 2.;
        
        --cnt;
    }
    
    V.swap(Vi);
    
    // update weights
    ++iter;
    
    //print current energies
    std::cout << "Energies:\n";
    std::cout << "  developable = " << currentDevelopabilityObjective.norm() << ",  |max| = " << currentDevelopabilityObjective.cwiseAbs().maxCoeff() << "\n";
    std::cout << "  proximity   = " << currentProximityObjective.norm() << "\n";
    std::cout << "  smooth      = " << currentSmoothnessObjective.norm() << "\n";
    std::cout << "  seam pair   = " << currentSeamProximityObjective.norm() << "\n";
    std::cout << "  seam smooth = " << currentSeamSmoothnessObjective.norm() << "\n";
    std::cout << "  boundary    = " << currentBoundaryProximityObjective.norm() << "\n";


    return 24 - cnt;
}
