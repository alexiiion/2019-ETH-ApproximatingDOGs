#include "ProcessSeams.hpp"
#include <igl/edges.h>
#include <igl/adjacency_list.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <Eigen/Sparse>
#include <iostream>
#include <set>
#include <unordered_map>

#include <igl/writeOFF.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/avg_edge_length.h>
#include <igl/remove_unreferenced.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_loop.h>

#include <igl/remove_duplicate_vertices.h>
#include <igl/facet_components.h>



ProcessSeams::ProcessSeams(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_, const Eigen::VectorXi& L_)
: V(V_), F(F_), L(L_)
{
}


std::vector<std::vector<int>> ProcessSeams::getSeamVertices()
{
    int n0 = 1 + I.maxCoeff();
    std::vector<std::vector<int>> copies(n0);

    for (int i = 0; i < I.size(); ++i)
    {
        if (I(i) != -1)
            copies[I(i)].push_back(i);
    }

    //for (auto& c : copies)
    //{
    //    if (c.size() < 1)
    //        continue;
    //
    //    const int nc = c.size();
    //    std::cout << "  c.size():" << c.size() << std::endl;
    //
    //    if (nc > 2)
    //        std::cout << "  --> more than 2: " << std::endl;
    //
    //    for (int i = 0; i < nc; ++i)
    //        std::cout << "    " << c[i] << ", ";
    //    std::cout << std::endl;
    //}

    return copies;
}

std::vector<std::array<int, 2>> ProcessSeams::getPairs()
{
    //int n0 = 1 + I.maxCoeff();
    //std::vector<std::vector<int>> copies(n0);
    //
    //for(int i = 0; i < I.size(); ++i)
    //{
    //    if(I(i) != -1) 
    //        copies[I(i)].push_back(i);
    //}

    std::vector<std::vector<int>> copies = getSeamVertices();
    std::vector<std::array<int, 2>> pairs;
    

    for(auto& c : copies)
    {
        if(c.size() > 1)
        {
            const int nc = c.size();

            for(int i = 0; i < nc; ++i)
                for(int j = i + 1; j < nc; ++j)
                    pairs.push_back(std::array<int, 2>{c[i], c[j]});
        }
    }
    
    return pairs;
}


void ProcessSeams::seamEdges(Eigen::MatrixXi& E)
{
    const int nf = F.rows();
    const int nv = V.rows();
    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(F, TT);
    
    std::vector<std::vector<int>> adj;
    std::vector<std::array<int, 2>> E0;
    
    igl::adjacency_list(F, adj);
    
    for(int i = 0; i < nf; ++i)
    {
        for(int k = 0; k < 3; ++k)
        {
            int lk = L(F(i, k));
            
            if(lk <  L(F(i, (k+1) % 3)) && lk < L(F(i, (k+2) % 3)))
            {
                E0.push_back(std::array<int, 2>{F(i, (k+1) % 3), F(i, (k+2) % 3)});
            }
        }
    }
    

    E.resize(E0.size(), 2);
    for(int i = 0; i < E0.size(); ++i)
    {
        E(i, 0) = E0[i][0];
        E(i, 1) = E0[i][1];
    }
}

void vertexSubmesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& sel, Eigen::MatrixXd& Vs, Eigen::MatrixXi& Fs, Eigen::VectorXi& I)
{
    const int nv = V.rows();
    const int nv2 = sel.rows();
    
    // count new faces
    int cf = 0;
    for(int i = 0; i < F.rows(); ++i)
    {
        if(sel(F(i, 0)) != -1 && sel(F(i, 1)) != -1&& sel(F(i, 2))  != -1) ++cf;
    }
    
    Fs.resize(cf, 3);
    int cnt = 0;
    
    for(int i = 0; i < F.rows(); ++i)
    {
        if(sel(F(i, 0)) != -1 && sel(F(i, 1)) != -1&& sel(F(i, 2))  != -1)
        {
            Fs(cnt, 0) = sel(F(i, 0));
            Fs(cnt, 1) = sel(F(i, 1));
            Fs(cnt, 2) = sel(F(i, 2));
            
            ++cnt;
        }
    }
    
    Vs.resize(nv2, 3);
    I.resize(nv2);
    
    for(int i = 0; i < nv; ++i)
    {
        if(sel(i) != -1)
        {
            Vs.row(sel(i)) = V.row(i);
            I(sel(i)) = i;
        }
    }
}

void boundaryFlag(const Eigen::MatrixXi& F, const int nv, std::vector<char>& b)
{
    b.clear();
    b.resize(nv, 0);
    
    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(F, TT);
    
    for(int i = 0; i < F.rows(); ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            if(TT(i, j) == -1)
            {
                b[F(i, j)] = 1;
                b[F(i, (j + 1) % 3)] = 1;
            }
        }
    }
}

template<typename T>
class SortedLess
{
public:
    bool operator()(const std::array<T, 2>& v0, const std::array<T, 2>& v1) const
    {
        if(v0[0] < v1[0]) return true;
        else if(v1[0] > v0[0]) return false;
        else return v0[1] < v1[0];
    }
};


void removeNAN(Eigen::SparseMatrix<double>& A)
{
    for(int i = 0; i < A.nonZeros(); ++i)
    {
        double& v = A.valuePtr()[i];
        
        if(!std::isfinite(v) || std::abs(v) > 1e16)
            v = .0;
    }
}




void ProcessSeams::splitMeshes(const double smoothFactor)
{
    std::vector<int> labels(L.data(), L.data() + L.size());
    std::sort(labels.begin(), labels.end());
    labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
    
    const int nc = labels.size();
    int nv = V.rows();
    const int nf = F.rows();

    std::vector<std::vector<int>> adj;
    igl::adjacency_list(F, adj);
    std::vector<Eigen::MatrixXi> Fs;
    
    // get faces containing at least one vertex of label i and vertices with label < i
    for(int j = 0; j < nc; ++j)
    {
        const int l = labels[j];
        
        std::vector<std::array<int, 3>> Fl;
        
        for(int i = 0; i < nf; ++i)
        {
            int lb[3]{L(F(i, 0)), L(F(i, 1)), L(F(i, 2))};
            
            bool foundL = false;
            bool greaterL = false;
            
            for(int k = 0; k < 3; ++k)
            {
                if(lb[k] == l) foundL = true;
                if(lb[k] > l) greaterL = true;
            }
            
            if(foundL && !greaterL)
            {
                Fl.push_back(std::array<int, 3>{F(i, 0), F(i, 1), F(i, 2)});
            }
        }
        
        Fs.push_back(Eigen::MatrixXi(Fl.size(), 3));
        
        for(int i = 0; i < Fl.size(); ++i)
        {
            Fs.back()(i, 0) = Fl[i][0];
            Fs.back()(i, 1) = Fl[i][1];
            Fs.back()(i, 2) = Fl[i][2];
        }
    }
    
    std::vector<Eigen::MatrixXd> Vs(nc);
    std::vector<Eigen::VectorXi> Is(nc);
        
    for(int i = 0; i < nc; ++i)
    {
        Eigen::VectorXi tmp;
        Eigen::MatrixXd NV;
        Eigen::MatrixXi NF;
        
        igl::remove_unreferenced(V, Fs[i], Vs[i], NF, tmp, Is[i]);
        Fs[i].swap(NF);
    }
    

    int nvMerged = 0;
    int nfMerged = 0;
    
   
    for(int i = 0; i < nc; ++i)
    {
        nvMerged += Vs[i].rows();
        nfMerged += Fs[i].rows();
    }
    
    // merge
    int offsetV = 0;
    int offsetF = 0;
    I.resize(nvMerged);
    V.resize(nvMerged, 3);
    F.resize(nfMerged, 3);
    L.resize(nvMerged);
        
    for(int l = 0; l < nc; ++l)
    {
        const int nvl = Vs[l].rows();
        
        I.block(offsetV, 0, nvl, 1) = Is[l];
        V.block(offsetV, 0, nvl, 3) = Vs[l];
        L.block(offsetV, 0, nvl, 1).setConstant(labels[l]);
        
        for(int i = 0; i < Fs[l].rows(); ++i)
        {
            for(int k = 0; k < 3; ++k)
                F(offsetF + i, k) = offsetV + Fs[l](i, k);
        }
        
        offsetF += Fs[l].rows();
        offsetV += nvl;
    }
}

void ProcessSeams::stitchMesh(const Eigen::MatrixXd& current_V, Eigen::MatrixXd& out_V, Eigen::MatrixXi& out_F)
{
    Eigen::MatrixXd V_copy = current_V;
    std::vector<std::vector<int>> seam_vertices = getSeamVertices();

    for (auto& s : seam_vertices)
    {
        const int nv = s.size();
        if (nv < 1)
            continue;

        //compute average location of the seam vertices
        Eigen::RowVector3d average(0, 0, 0);
        for (int i = 0; i < nv; ++i)
            average += V_copy.row(s[i]);

        average *= 1. / nv;

        //snap seam vertices to the average location
        for (int i = 0; i < nv; ++i)
            V_copy.row(s[i]) = average;
    }

    Eigen::MatrixXd newV; Eigen::MatrixXi newVI, newVJ, newF;
    igl::remove_duplicate_vertices(V_copy, this->F, 1e-7, newV, newVI, newVJ, newF);

    // count connected components of the stitched mesh
    Eigen::MatrixXi C;
    igl::facet_components(newF, C);
    std::cout << "Output has " << newV.rows() << " vertices and " << newF.rows() << " faces" << " with " << C.maxCoeff() + 1 << " connected components" << std::endl;

    out_V = newV;
    out_F = newF;
}