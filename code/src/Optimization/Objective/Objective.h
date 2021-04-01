#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>

class Objective 
{
public:
  virtual Objective* clone() const = 0;

  virtual double obj(const Eigen::VectorXd& x) const = 0;
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const = 0;

  int get_hessian_IJV_size() 
  {
	  return IJV.size();
  }
  
  const std::vector<Eigen::Triplet<double>>& update_and_get_hessian_ijv(const Eigen::VectorXd& x) 
  {
    updateHessianIJV(x); 
	return IJV; 
  }
  
  // For objectives that depeends on a state/reference shape (which might also be updated at various times)
  virtual void set_ref(const Eigen::VectorXd& x0) {};

  void check_grad(const Eigen::VectorXd& x) const;

  std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M) 
  {
    std::vector<Eigen::Triplet<double>> v;
  
	for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    
	return v;
  }

  // builds the hessian from an IJV
  virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x) 
  {
    //Eigen::SparseMatrix<double> H(x.rows(),x.rows());
    updateHessianIJV(x);
  
	if (!is_H_cached) {
      cachedH =  Eigen::SparseMatrix<double>(x.rows(),x.rows());
      igl::sparse_cached_precompute(IJV, cached_ijv_data, cachedH);
      is_H_cached = true;
    } else {
      igl::sparse_cached(IJV, cached_ijv_data, cachedH);
    }
    
	return cachedH;
  }


 protected:
   std::vector<Eigen::Triplet<double> > IJV; 
   Eigen::VectorXi cached_ijv_data;
   
   Eigen::SparseMatrix<double> cachedH;
   bool is_H_cached = false;


 private:
  // return 0 sparse matrix if not implemented
  virtual void updateHessianIJV(const Eigen::VectorXd& x) { /*empty on purpose */ } 

  void finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy = 1) const;
};