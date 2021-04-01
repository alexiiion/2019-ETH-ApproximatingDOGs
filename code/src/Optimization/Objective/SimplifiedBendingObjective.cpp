#include "SimplifiedBendingObjective.h"

#include "DataModel.h"
#include "Logger.h"

SimplifiedBendingObjective::SimplifiedBendingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x)  : quadTop(quadTop) {
	// Number of hessian triplets
	IJV.resize(51*quadTop.stars.rows()/5+27*quadTop.bnd3.rows()/4);
	// Change the weights to 1 to get to the old bending energy (this is for unregular grid, needed for intersecting creases)
	int vnum = x.rows()/3;
	init_edge_lengths.resize(4*quadTop.stars.rows()/5 + 2*quadTop.bnd3.rows()/4); 
	int cnt = 0;
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double ex_f_l = sqrt( pow(p0_x-pxf_x,2) + pow(p0_y-pxf_y,2) + pow(p0_z-pxf_z,2) );
		double ex_b_l = sqrt( pow(p0_x-pxb_x,2) + pow(p0_y-pxb_y,2) + pow(p0_z-pxb_z,2) );
		double ey_f_l = sqrt( pow(p0_x-pyf_x,2) + pow(p0_y-pyf_y,2) + pow(p0_z-pyf_z,2) );
		double ey_b_l = sqrt( pow(p0_x-pyb_x,2) + pow(p0_y-pyb_y,2) + pow(p0_z-pyb_z,2) );
		
		// Normalize weights to 1 and choose the ratio as the opposite of the ratio (to have linear precision)
		init_edge_lengths[cnt++] = ex_b_l/(pow(ex_f_l,2)+pow(ex_b_l,2));
		init_edge_lengths[cnt++] = ex_f_l/(pow(ex_f_l,2)+pow(ex_b_l,2));
		init_edge_lengths[cnt++] = ey_b_l/(pow(ey_b_l,2)+pow(ey_f_l,2));
		init_edge_lengths[cnt++] = ey_f_l/(pow(ey_b_l,2)+pow(ey_f_l,2));
	}

	 for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	    double ex_f_l = sqrt( pow(p0_x-pxf_x,2) + pow(p0_y-pxf_y,2) + pow(p0_z-pxf_z,2) );
		double ex_b_l = sqrt( pow(p0_x-pxb_x,2) + pow(p0_y-pxb_y,2) + pow(p0_z-pxb_z,2) );

		init_edge_lengths[cnt++] = ex_b_l/(pow(ex_f_l,2)+pow(ex_b_l,2));
		init_edge_lengths[cnt++] = ex_f_l/(pow(ex_f_l,2)+pow(ex_b_l,2));
	}
}

double SimplifiedBendingObjective::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  int vnum = x.rows()/3;
  
  int cnt = 0;  
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t2 = len_ex_b*p0_x+len_ex_f*p0_x-len_ex_b*pxb_x-len_ex_f*pxf_x;
		double t3 = len_ex_b*p0_y+len_ex_f*p0_y-len_ex_b*pxb_y-len_ex_f*pxf_y;
		double t4 = len_ex_b*p0_z+len_ex_f*p0_z-len_ex_b*pxb_z-len_ex_f*pxf_z;
		double t5 = len_ey_b*p0_x+len_ey_f*p0_x-len_ey_b*pyb_x-len_ey_f*pyf_x;
		double t6 = len_ey_b*p0_y+len_ey_f*p0_y-len_ey_b*pyb_y-len_ey_f*pyf_y;
		double t7 = len_ey_b*p0_z+len_ey_f*p0_z-len_ey_b*pyb_z-len_ey_f*pyf_z;
		e += t2*t2+t3*t3+t4*t4+t5*t5+t6*t6+t7*t7;


  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

    double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];
    
	double t2 = len_ex_b*p0_x+len_ex_f*p0_x-len_ex_b*pxb_x-len_ex_f*pxf_x;
	double t3 = len_ex_b*p0_y+len_ex_f*p0_y-len_ex_b*pxb_y-len_ex_f*pxf_y;
	double t4 = len_ex_b*p0_z+len_ex_f*p0_z-len_ex_b*pxb_z-len_ex_f*pxf_z;
	e += t2*t2+t3*t3+t4*t4;

  }

  DataModel::output_bending_energy = e;
  write_log(DataModel::log_level_optimization) << "    SimplifiedBendingObjective: " << DataModel::output_bending_energy << std::endl;

  return e;
}

Eigen::VectorXd SimplifiedBendingObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;

  int cnt = 0;
  #pragma clang loop vectorize(enable)
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
        //local_grad.setZero();

        //int p_0_i = i*s+j, p_xf_i = i*s+j+1,p_xb_i = i*s+j-1,p_yf_i = (i+1)*s+j, p_yb_i = (i-1)*s+j;
        int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
        const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
        const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  const double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t2 = len_ex_b+len_ex_f;
		double t3 = len_ey_b+len_ey_f;
		double t4 = len_ex_b*p0_x;
		double t5 = len_ex_f*p0_x;
		double t13 = len_ex_b*pxb_x;
		double t14 = len_ex_f*pxf_x;
		double t6 = t4+t5-t13-t14;
		double t7 = len_ex_b*p0_y;
		double t8 = len_ex_f*p0_y;
		double t15 = len_ex_b*pxb_y;
		double t16 = len_ex_f*pxf_y;
		double t9 = t7+t8-t15-t16;
		double t10 = len_ex_b*p0_z;
		double t11 = len_ex_f*p0_z;
		double t17 = len_ex_b*pxb_z;
		double t18 = len_ex_f*pxf_z;
		double t12 = t10+t11-t17-t18;
		double t19 = len_ey_b*p0_x;
		double t20 = len_ey_f*p0_x;
		double t28 = len_ey_b*pyb_x;
		double t29 = len_ey_f*pyf_x;
		double t21 = t19+t20-t28-t29;
		double t22 = len_ey_b*p0_y;
		double t23 = len_ey_f*p0_y;
		double t30 = len_ey_b*pyb_y;
		double t31 = len_ey_f*pyf_y;
		double t24 = t22+t23-t30-t31;
		double t25 = len_ey_b*p0_z;
		double t26 = len_ey_f*p0_z;
		double t32 = len_ey_b*pyb_z;
		double t33 = len_ey_f*pyf_z;
		double t27 = t25+t26-t32-t33;

		grad(p_0_i) += t2*t6*2.0+t3*t21*2.0;
		grad(p_0_i+v_num) += t2*t9*2.0+t3*t24*2.0;
		grad(p_0_i+2*v_num) += t2*t12*2.0+t3*t27*2.0;
		grad(p_xb_i+0) +=len_ex_b*t6*-2.0;
		grad(p_xb_i+v_num) += len_ex_b*t9*-2.0;
		grad(p_xb_i+2*v_num) += len_ex_b*t12*-2.0;
		grad(p_xf_i) += len_ex_f*t6*-2.0;
		grad(p_xf_i+v_num) += len_ex_f*t9*-2.0;
		grad(p_xf_i+2*v_num) += len_ex_f*t12*-2.0;
		grad(p_yb_i) += len_ey_b*t21*-2.0;
		grad(p_yb_i+v_num) += len_ey_b*t24*-2.0;
		grad(p_yb_i+2*v_num) += len_ey_b*t27*-2.0;
		grad(p_yf_i) += len_ey_f*t21*-2.0;
		grad(p_yf_i+v_num) += len_ey_f*t24*-2.0;
		grad(p_yf_i+2*v_num) += len_ey_f*t27*-2.0;
  }


  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];

	double t2 = len_ex_b+len_ex_f;
	double t3 = len_ex_b*p0_x;
	double t4 = len_ex_f*p0_x;
	double t12 = len_ex_b*pxb_x;
	double t13 = len_ex_f*pxf_x;
	double t5 = t3+t4-t12-t13;
	double t6 = len_ex_b*p0_y;
	double t7 = len_ex_f*p0_y;
	double t14 = len_ex_b*pxb_y;
	double t15 = len_ex_f*pxf_y;
	double t8 = t6+t7-t14-t15;
	double t9 = len_ex_b*p0_z;
	double t10 = len_ex_f*p0_z;
	double t16 = len_ex_b*pxb_z;
	double t17 = len_ex_f*pxf_z;
	double t11 = t9+t10-t16-t17;

	grad(p_0_i) +=t2*t5*2.0;
	grad(p_0_i+v_num) += t2*t8*2.0;
	grad(p_0_i+2*v_num) += t2*t11*2.0;
	grad(p_xb_i+0) += len_ex_b*t5*-2.0;
	grad(p_xb_i+v_num) += len_ex_b*t8*-2.0;
	grad(p_xb_i+2*v_num) += len_ex_b*t11*-2.0;
	grad(p_xf_i) += len_ex_f*t5*-2.0;
	grad(p_xf_i+v_num) += len_ex_f*t8*-2.0;
	grad(p_xf_i+2*v_num) += len_ex_f*t11*-2.0;

  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  //cout << "new grad.norm() = " << grad.norm() << endl; //exit(1);
  return grad;
}


void SimplifiedBendingObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  // Number of ijv values
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0; int cnt = 0;
  #pragma clang loop vectorize(enable)
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
        //local_grad.setZero();

        //int p_0_i = i*s+j, p_xf_i = i*s+j+1,p_xb_i = i*s+j-1,p_yf_i = (i+1)*s+j, p_yb_i = (i-1)*s+j;
        const int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
        const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
        const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  const double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t2 = len_ex_b+len_ex_f;
		double t3 = len_ey_b+len_ey_f;
		double t4 = t2*t2;
		double t5 = t4*2.0;
		double t6 = t3*t3;
		double t7 = t6*2.0;
		double t8 = t5+t7;
		double t9 = len_ex_b*len_ex_b;
		double t10 = t9*2.0;
		double t11 = len_ex_b*len_ex_f*2.0;
		double t12 = len_ex_f*len_ex_f;
		double t13 = t12*2.0;
		double t14 = len_ey_b*len_ey_b;
		double t15 = t14*2.0;
		double t16 = len_ey_b*len_ey_f*2.0;
		double t17 = len_ey_f*len_ey_f;
		double t18 = t17*2.0;

		
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i, t8);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xb_i, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_yb_i, len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_yf_i, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t8);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_yb_i+vnum, len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_yf_i+vnum, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t8);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yb_i+2*vnum,len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yf_i+2*vnum, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_0_i, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t10);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t10);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t10);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t13);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t13);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t11);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t13);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_0_i, len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_yb_i, t15);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_yf_i, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_0_i+vnum, len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yb_i+vnum, t15);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yf_i+vnum, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_0_i+2*vnum, len_ey_b*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yb_i+2*vnum, t15);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yf_i+2*vnum, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_0_i, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_yb_i, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_yf_i, t18);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_0_i+vnum, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yb_i+vnum, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yf_i+vnum, t18);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_0_i+2*vnum, len_ey_f*t3*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yb_i+2*vnum, t16);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yf_i+2*vnum, t18);

   }

   for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
		const int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];

		double t2 = len_ex_b+len_ex_f;
		double t3 = t2*t2;
		double t4 = t3*2.0;
		double t5 = len_ex_b*len_ex_b;
		double t6 = t5*2.0;
		double t7 = len_ex_b*len_ex_f*2.0;
		double t8 = len_ex_f*len_ex_f;
		double t9 = t8*2.0;

		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i, t4);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xb_i, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t4);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t4);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_0_i, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t6);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t6);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, len_ex_b*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t6);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t9);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t9);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, len_ex_f*t2*-2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t7);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t9);
	}
}