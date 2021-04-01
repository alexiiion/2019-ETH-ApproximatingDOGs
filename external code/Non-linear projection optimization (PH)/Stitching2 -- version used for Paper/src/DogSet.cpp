#include "DogSet.hpp"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/facet_components.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>
#include <igl/edges.h>
#include <igl/upsample.h>

#include <igl/PI.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <iostream>

#define GCO_ENERGYTYPE long long
//#undef GCO_ENERGYTYPE32
#include "GCoptimization.h"

Eigen::MatrixXd D;
Eigen::MatrixXd N;
Eigen::MatrixXd DN; //currently unused
std::vector<Eigen::MatrixXd> P; //COPY of projected, unite!


//#define CUTOLD


#ifdef CUTOLD
double scale = 10000000;
double cover_smoothness = 50;
double max_value;
double max_distance = 30;
#else
double scale = 1000000;
double cover_smoothness = 1000;
//double escale = 2 * cover_smoothness;
double escale = 0.1 * cover_smoothness;
#endif

//double weight_normals = 0.0;

GCO_ENERGYTYPE smoothnessTerm(const int n1, const int n2, const int l1, const int l2)
{
    if (l1 == l2)
        return 0;
   
	if(!DogSet::use_distance_smoothness)
		return cover_smoothness;
    
    
    
    
    //double data = 2000000;// + std::pow(std::max(D(n1, l1), D(n2, l2)), 4);
    
    Eigen::RowVector3d projected_point1 = P[l1].row(n1);
    Eigen::RowVector3d projected_point2 = P[l2].row(n2);
    double squared_distance = (projected_point1 - projected_point2).squaredNorm();
    squared_distance *= escale;

    long data = cover_smoothness + (1 / (squared_distance)) * scale;

    if (data < 0 || data > GCO_MAX_ENERGYTERM)
        return GCO_MAX_ENERGYTERM;

    return data;

  
}

GCO_ENERGYTYPE dataTerm(const int node, const int label)
{
    double distance = D(node, label);
    double normal_deviation = N(node, label);

    //double distance_normal = DN(node, label);
    //double d = (distance + (distance_normal * distance_normal)) / (1*1 + 2*2);
    //double d = (distance + distance_normal) / (1 + 2);

    normal_deviation = std::pow(normal_deviation, 3);
	double wn = DogSet::weight_normals; //0.0;
    double d = (1.0-wn)*distance + wn*normal_deviation;
    long data = d * scale;
  
    //int data = distance * scale;
    //std::cout << "dataTerm() data = " << data << "\n";

    if (data < 0 || data > GCO_MAX_ENERGYTERM)
        return GCO_MAX_ENERGYTERM;

    return data;
    
    //double distance = D(node, label);
    //int data = distance * scale;
    //  
    //if (data < 0 || data > GCO_MAX_ENERGYTERM)
    //    return GCO_MAX_ENERGYTERM;
    //
    //return data;
}


void DogSet::init()
{
    Eigen::VectorXi I;
    igl::facet_components(Fall, I);
    
    n = 1 + I.maxCoeff();
    F.resize(n);
    V.resize(n);
    trees.resize(n);
    projected.resize(n);
    
    std::vector<int> cnts(n, 0);
    
    for(int i = 0; i < I.size(); ++i)
        ++cnts[I(i)];
    
    for(int i = 0; i < n; ++i)
        F[i].resize(cnts[i], 3);
    
    std::fill_n(cnts.data(), n, .0);
    
    for(int i = 0; i < I.size(); ++i)
    {
        int pos = cnts[I(i)]++;
        F[I(i)].row(pos) = Fall.row(i);
    }
    
    Eigen::VectorXi tmp1;
    for(int i = 0; i < n; ++i)
    {
        Eigen::MatrixXi Fi;
        Fi.swap(F[i]);
        igl::remove_unreferenced(Vall, Fi, V[i], F[i], tmp1);
        
        trees[i].init(V[i], F[i]);
    }
}

DogSet::DogSet(std::string fname)
{
    igl::readOBJ(fname, Vall, Fall);
    init();
}

DogSet::DogSet(const Eigen::MatrixXd& V0, const Eigen::MatrixXi& F0)
{
    Vall = V0;
    Fall = F0;
    init();
}

void DogSet::projectPoint(Eigen::RowVector3d& p, const int l)
{
    int j;
    Eigen::RowVector3d p2;
    trees[l].squared_distance(V[l], F[l], p, j, p2);
    
    p.swap(p2);
}


void DogSet::projectPoints(const Eigen::MatrixXd& Vin, const Eigen::VectorXi& L, Eigen::MatrixXd& VP) const
{
    const int nv = Vin.rows();
    assert(nv == L.rows());
    
    VP.resize(nv, 3);
    
    for(int i = 0; i < nv; ++i)
    {
        const int l = L(i);
        int j;
        Eigen::RowVector3d p;
        trees[l].squared_distance(V[l], F[l], Vin.row(i), j, p);
        
        VP.row(i) = p;
    }
}

void DogSet::computeDistances(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin)
{
    int nv = Vin.rows();
    D.resize(nv, n);
    
    Eigen::VectorXi tmp1;
    Eigen::MatrixXd tmp2;
    
    for(int i = 0; i < n; ++i)
    {
        Eigen::VectorXd Di;
        Eigen::MatrixXd VPi;
        
        igl::point_mesh_squared_distance(Vin, V[i], F[i], Di, tmp1, VPi);
        D.col(i).swap(Di);
        projected[i].swap(VPi);
    }
    
    P.clear();
    P = projected;

#ifdef CUTOLD
        D = D.cwiseSqrt();
        D /= D.maxCoeff();
        
        for (int i = 0; i < D.rows(); ++i)
        {
            for (int j = 0; j < D.cols(); j++)
            {
                        D(i, j) = std::min(std::exp(D(i, j) * 10), 100000.0);
              //  D(i, j) = std::min(std::exp(D(i, j) * 150), 1000000.0);
            }
        }
    
        cover_smoothness = 50;//smoothness;
        scale = 10000000;
        //max_distance = std::exp(1.50 * 10);
    
        max_value = D.maxCoeff();
        scale /= max_value;
    
#else
    
    D /= D.maxCoeff();
    
#endif
    
}

void DogSet::computeDistancesAlongMeshNormals(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin)
{
    int nv = Vin.rows();
    DN.resize(nv, n);

    Eigen::MatrixXd normals;
    igl::per_vertex_normals(Vin, Fin, normals);

    Eigen::VectorXi tmp1;
    Eigen::MatrixXd tmp2;

    std::vector<std::vector<bool>> is_intersecting(nv);

    for (int vi = 0; vi < nv; ++vi)
    {
        is_intersecting[vi] = std::vector<bool>(n, false);

        for (int i = 0; i < n; ++i)
        {
            igl::Hit hit_up;
            igl::Hit hit_down;
            bool found_hit_up = igl::ray_mesh_intersect(Vin.row(vi), normals.row(vi), V[i], F[i], hit_up);
            bool found_hit_down = igl::ray_mesh_intersect(Vin.row(vi), -normals.row(vi), V[i], F[i], hit_down);
            
            double Di = found_hit_up || found_hit_down ? std::min(hit_up.t, hit_down.t) : 0;
            is_intersecting[vi][i] = found_hit_up || found_hit_down;
            //std::cout << "ray hit found? " << std::boolalpha << (found_hit_up || found_hit_down) << " Di = " << Di << " (hit_up.t = " << hit_up.t << ", hit_down.t = " << hit_down.t << "\n";

            DN(vi, i) = Di;
        }
    }

    DN /= DN.maxCoeff();

    for (int vi = 0; vi < nv; ++vi)
    {
        for (int i = 0; i < n; ++i)
        {
            if (!is_intersecting[vi][i])
                DN(vi, i) = 2;
        }
    }
}

void DogSet::computeNormalDeviations(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin)
{
    int nv = Vin.rows();
    N.resize(nv, n);

    Eigen::MatrixXd normals_target;
    igl::per_vertex_normals(Vin, Fin, normals_target);

    for (int i = 0; i < n; ++i)
    {
        Eigen::VectorXi closest_faces;
        Eigen::VectorXd square_distance_unused;
        Eigen::MatrixXd closest_points_unused;
        igl::point_mesh_squared_distance(Vin, V[i], F[i], square_distance_unused, closest_faces, closest_points_unused);

        Eigen::MatrixXd normals_patch;
        igl::per_face_normals(V[i], F[i], normals_patch);

        for (int vi = 0; vi < nv; ++vi)
        {
            int npi = closest_faces(vi);
            Eigen::VectorXd np = normals_patch.row(npi);
            Eigen::VectorXd nt = normals_target.row(vi);

            N(vi, i) = acos(nt.dot(np));
            //std::cout << "N = " << N(vi, i) << "\n";
        }

    }

    //std::cout << "N.minCoeff() = " << N.minCoeff() << ", N.maxCoeff() = " << N.maxCoeff() << "\n";
    //N /= N.maxCoeff();
    N /= igl::PI;
    std::cout << "N.minCoeff() = " << N.minCoeff() << ", N.maxCoeff() = " << N.maxCoeff() << "\n";
}


void DogSet::getPointsFromLabels(const Eigen::VectorXi& L, Eigen::MatrixXd& Vout)
{
    Vout.resize(L.rows(), 3);
    
    for(int i = 0; i < Vout.rows(); ++i)
    {
        Vout.row(i) = projected[L(i)].row(i);
    }
}

void DogSet::minimalAssignement(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin)
{
    const int nv = Vin.rows();
    assert(D.rows() == nv);
    I.resize(nv);
    
    for(int i = 0; i < D.rows(); ++i)
    {
        int m;
        D.row(i).minCoeff(&m);
        
        I(i) = m;
    }
}

void DogSet::graphCutAssignement(const Eigen::MatrixXd& Vin, const Eigen::MatrixXi& Fin, const double smoothnessScale)
{
    
    computeDistances(Vin, Fin);
    computeDistancesAlongMeshNormals(Vin, Fin);
    computeNormalDeviations(Vin, Fin);
    
    const int nv = Vin.rows();
    GCoptimizationGeneralGraph mrf(nv, n);
  
    Eigen::MatrixXi E;
    igl::edges(Fin, E);
    
#ifndef CUTOLD
    
    Eigen::VectorXd L(E.rows());
    
    for (int i = 0; i < E.rows(); ++i)
        L(i) = (Vin.row(E(i, 0)) - Vin.row(E(i, 1))).norm();
    
    L /= L.maxCoeff();
    
    escale = smoothnessScale * cover_smoothness;
    
    for (int i = 0; i < E.rows(); ++i)
    {
        mrf.setNeighbors(E(i, 0), E(i, 1), L(i) * escale);
    }
    
    mrf.setLabelCost((GCO_ENERGYTYPE)0 );
    //   mrf.setLabelCost((GCO_ENERGYTYPE)100000);
#else
    
    for (int i = 0; i < E.rows(); ++i)
     {
         mrf.setNeighbors(E(i, 0), E(i, 1));
     }
    
    mrf.setSmoothCost(smoothnessTerm);
    mrf.setLabelCost((GCO_ENERGYTYPE)0);
#endif
    const int seed = time(0);
    std::cout << "seed: " << seed << "\n";
    srand(seed);
    
    mrf.setDataCost(dataTerm);
   // mrf.setSmoothCost(smoothnessTerm);
    

    mrf.setVerbosity(1);
  //  mrf.setLabelOrder(true);
    
    
    std::vector<int> order(n);
    for(int i = 0; i < n; ++i) order[i] = i;
    std::random_shuffle(order.begin(), order.end());
    mrf.setLabelOrder(order.data(), order.size());
    
    
   // std::cout << "    start optimization\n";
    auto energy = mrf.expansion();

    auto en_data = mrf.giveDataEnergy();
    auto en_smooth = mrf.giveSmoothEnergy();
    auto en_label = mrf.giveLabelEnergy();
    
  //  std::cout << "    Energy: " << energy << ", Data Energy: " << en_data << ", Smoothness Energy: " << en_smooth << ", Label Energy: " << en_label << std::endl;
    
    I.resize(nv);
    mrf.whatLabel(0, I.size(), I.data());
}


//void DogSet::set_weight_normals(double& weight)
//{
//	weight_normals = weight;
//}