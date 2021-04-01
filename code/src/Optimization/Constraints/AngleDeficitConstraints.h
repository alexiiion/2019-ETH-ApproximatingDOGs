#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

#include "DataModel.h"
#include "Utils.h"
#include "Logger.h"

class AngleDeficitConstraints : public Constraints {
public:
	AngleDeficitConstraints(const std::vector<int>& inner_v, const std::vector<std::vector<int>>& inner_v_nbs, const Eigen::VectorXd& mass_matrix) 
		: inner_v(inner_v), inner_v_nbs(inner_v_nbs) 
	{
		const_n = inner_v.size();

		int ijv_size = 0;
		for (auto nb_list : inner_v_nbs) 
			ijv_size += 9 * nb_list.size();
		IJV.resize(ijv_size);
		
		sqrtM.resize(mass_matrix.rows());
		for (int i = 0; i < mass_matrix.rows(); i++) 
			sqrtM(i) = sqrt(mass_matrix(i));


		if(sqrtM.hasNaN())
			write_log(1) << linebreak << "ERROR in <" << __FUNCTION__ << "> AngleDeficitConstraints constructor:: sqrtM.hasNaN() = true " << sqrtM.hasNaN() << std::endl << std::endl;
		//write_log(0) << "AngleDeficitConstraints constructor:: sqrtM: " << linebreak << sqrtM << std::endl;
	};

	virtual AngleDeficitConstraints* clone() const { return new AngleDeficitConstraints(*this); }

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const 
	{
		int vnum = x.rows() / 3;
		Eigen::VectorXd constraints_vals(const_n);
		for (int i = 0; i < inner_v.size(); i++) 
		{
			double angle_sum = 0;

			int p0_i = inner_v[i];
			const double p0_x(x(p0_i + 0)); const double p0_y(x(p0_i + 1 * vnum)); const double p0_z(x(p0_i + 2 * vnum));
			
			int number_neighbors = inner_v_nbs[i].size();
			for (int j = -1; j < number_neighbors-1; j++) 
			{
				int current_index = j == -1 ? number_neighbors - 1 : j;
				int next_index = j+1;

				int p_xf_i = inner_v_nbs[i][current_index], p_yf_i = inner_v_nbs[i][next_index];
				//int p_xf_i = inner_v_nbs[i][j], p_yf_i = inner_v_nbs[i][j + 1];

				const double pxf_x(x(p_xf_i + 0)); const double pxf_y(x(p_xf_i + 1 * vnum)); const double pxf_z(x(p_xf_i + 2 * vnum));
				const double pyf_x(x(p_yf_i + 0)); const double pyf_y(x(p_yf_i + 1 * vnum)); const double pyf_z(x(p_yf_i + 2 * vnum));

				double t2 = p0_x - pxf_x;
				double t3 = p0_y - pxf_y;
				double t4 = p0_z - pxf_z;
				double t5 = p0_x - pyf_x;
				double t6 = p0_y - pyf_y;
				double t7 = p0_z - pyf_z;
				//double cur_angle = acos((t2*t5 + t3 * t6 + t4 * t7)*1.0 / sqrt(t2*t2 + t3 * t3 + t4 * t4)*1.0 / sqrt(t5*t5 + t6 * t6 + t7 * t7));
				double cos_cur_angle = clip((t2*t5 + t3*t6 + t4*t7)*1.0 / sqrt(t2*t2 + t3*t3 + t4*t4)*1.0 / sqrt(t5*t5 + t6*t6 + t7*t7), -1, 1);
				double cur_angle = acos(cos_cur_angle);
				
				if(std::isnan(cur_angle))
					write_log(2) << linebreak << "WARNING in <" << __FUNCTION__ << "> current angle is NaN!" << std::endl << std::endl;
				
				angle_sum += cur_angle;
			}

			write_log(DataModel::log_level_optimization + 1) << "      angle_sum = " << angle_sum << '\n';
			constraints_vals(i) = sqrtM(i)*(angle_sum - 2 * M_PI);
		}

		if(constraints_vals.hasNaN())
			write_log(DataModel::log_level_optimization) << "    AngleDeficitConstraints norm: " << constraints_vals.norm() << std::endl;
			

		write_log(DataModel::log_level_optimization) << "    AngleDeficitConstraints norm: " << constraints_vals.norm() << std::endl;
		return constraints_vals;
	}

	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int vnum = x.rows() / 3;
		int ijv_count = 0;
		for (int i = 0; i < inner_v.size(); i++) {
			int p0_i = inner_v[i];
			const double p0_x(x(p0_i + 0)); const double p0_y(x(p0_i + 1 * vnum)); const double p0_z(x(p0_i + 2 * vnum));

			int number_neighbors = inner_v_nbs[i].size();
			for (int j = -1; j < number_neighbors - 1; j++)
			{
				int current_index = j == -1 ? number_neighbors - 1 : j;
				int next_index = j + 1;

				int p_xf_i = inner_v_nbs[i][current_index], p_yf_i = inner_v_nbs[i][next_index];
				//int p_xf_i = inner_v_nbs[i][j], p_yf_i = inner_v_nbs[i][j + 1];

				const double pxf_x(x(p_xf_i + 0)); const double pxf_y(x(p_xf_i + 1 * vnum)); const double pxf_z(x(p_xf_i + 2 * vnum));
				const double pyf_x(x(p_yf_i + 0)); const double pyf_y(x(p_yf_i + 1 * vnum)); const double pyf_z(x(p_yf_i + 2 * vnum));

				double t3 = p0_x - pxf_x;
				double t4 = p0_y - pxf_y;
				double t5 = p0_z - pxf_z;
				double t6 = p0_x - pyf_x;
				double t7 = p0_y - pyf_y;
				double t8 = p0_z - pyf_z;
				double t17 = t3 * t6;
				double t18 = t4 * t7;
				double t19 = t5 * t8;
				double t2 = t17 + t18 + t19;
				double t9 = t3 * t3;
				double t10 = t4 * t4;
				double t11 = t5 * t5;
				double t12 = t9 + t10 + t11;
				double t13 = t6 * t6;
				double t14 = t7 * t7;
				double t15 = t8 * t8;
				double t16 = t13 + t14 + t15;
				double t20 = 1.0 / sqrt(t16);
				double t21 = p0_x * 2.0;
				double t22 = 1.0 / sqrt(t12);
				double t23 = t2 * t2;
				double t24 = 1.0 / t12;
				double t25 = 1.0 / t16;
				double t31 = t23 * t24*t25;
				double t26 = -t31 + 1.0;
				double t27 = 1.0 / sqrt(t26);
				double t28 = 1.0 / pow(t12, 3.0 / 2.0);
				double t29 = p0_y * 2.0;
				double t30 = 1.0 / pow(t16, 3.0 / 2.0);
				double t32 = p0_z * 2.0;
				double t33 = pxf_x * 2.0;
				double t34 = t21 - t33;
				double t35 = t2 * t20*t28*t34*(1.0 / 2.0);
				double t36 = pxf_y * 2.0;
				double t37 = t29 - t36;
				double t38 = t2 * t20*t28*t37*(1.0 / 2.0);
				double t39 = pxf_z * 2.0;
				double t40 = t32 - t39;
				double t41 = t2 * t20*t28*t40*(1.0 / 2.0);
				double t42 = pyf_x * 2.0;
				double t43 = t21 - t42;
				double t44 = t2 * t22*t30*t43*(1.0 / 2.0);
				double t45 = pyf_y * 2.0;
				double t46 = t29 - t45;
				double t47 = t2 * t22*t30*t46*(1.0 / 2.0);
				double t48 = pyf_z * 2.0;
				double t49 = t32 - t48;
				double t50 = t2 * t22*t30*t49*(1.0 / 2.0);

				IJV[ijv_count++] = Eigen::Triplet<double>(i, p0_i + 0, sqrtM(i)*(t27*(t35 + t44 + t20 * t22*(pxf_x + pyf_x - t21))));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p0_i + 1 * vnum, sqrtM(i)*(t27*(t38 + t47 + t20 * t22*(pxf_y + pyf_y - t29))));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p0_i + 2 * vnum, sqrtM(i)*(t27*(t41 + t50 + t20 * t22*(pxf_z + pyf_z - t32))));

				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_xf_i + 0, sqrtM(i)*(-t27 * (t35 - t6 * t20*t22)));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_xf_i + 1 * vnum, sqrtM(i)*(-t27 * (t38 - t7 * t20*t22)));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_xf_i + 2 * vnum, sqrtM(i)*(-t27 * (t41 - t8 * t20*t22)));

				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_yf_i + 0, sqrtM(i)*(-t27 * (t44 - t3 * t20*t22)));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_yf_i + 1 * vnum, sqrtM(i)*(-t27 * (t47 - t4 * t20*t22)));
				IJV[ijv_count++] = Eigen::Triplet<double>(i, p_yf_i + 2 * vnum, sqrtM(i)*(-t27 * (t50 - t5 * t20*t22)));
			}
		}
	}

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Empty on purpose (gauss newton hessian approximation)
	};

private:
	std::vector<int> inner_v;
	std::vector<std::vector<int>> inner_v_nbs;
	Eigen::VectorXd sqrtM;
};
