
#include "../include/Laser.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Runge_Kutta.cuh"

#include <cmath>









__device__ double CalculationE1(const particle& first, const particle& second)
{
	//坐标平方和
	const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	//const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	//const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);


	return  0.5 * px_py_pz_squre_sum_first - 2.0 / sqrt(loc_squre_sum_first) +
		1.0 / loc_squre_sum_first * Q_constant * Q_constant / 4.0 / A_hardness *
		exp(A_hardness * (1.0 - pow((loc_squre_sum_first * px_py_pz_squre_sum_first /
			                            Q_constant * Q_constant), 2))) +
		1.0 / sqrt(distance_squre) / 2.0;
}

__device__ double CalculationE2(const particle& first, const particle& second)
{
	//坐标平方和
	//const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	//const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);


	return  0.5 * px_py_pz_squre_sum_second - 2.0 / sqrt(loc_squre_sum_second) +
		1.0 / loc_squre_sum_second * Q_constant * Q_constant / 4.0 / A_hardness *
		exp(A_hardness * (1.0 - pow((loc_squre_sum_second * px_py_pz_squre_sum_second /
			                            Q_constant * Q_constant), 2))) +
		1.0 / sqrt(distance_squre) / 2.0;
}









__device__ double compute_qq_single(const size_t& now_step)
{
	double t1 = 0.5 * DX * (now_step + 1);
	return  pow((sin(Omega1 / 2.0 / (2 * N1_const + N2_const)*t1)), 2);

}

__device__ double compute_e_for_check(const size_t& now_step, const double& e1_single, const double& e2_single)
{
	return  sqrt(pow(e1_single, 2) + pow(e2_single, 2));
}

__device__ double compute_e1_single(const size_t& now_step, const double& qq_now_single, const double& EE0)
{
	double tao = 0.0;
	double t1 = 0.5 * DX * (now_step + 1);
	return  (EE0 / (1.0 + TP_const)) * qq_now_single * sin(Omega1 * t1 + tao) -
		(EE0*TP_const / (1.0 + TP_const)) * qq_now_single * sin(Omega2 * t1 + 2 * tao);
}

__device__ double compute_e2_single(const size_t& now_step, const double& qq_now_single, const double& EE0)
{
	double tao = 0.0;
	double t1 = 0.5 * DX * (now_step + 1);
	return  (EE0 / (1.0 + TP_const)) * qq_now_single * cos(Omega1 * t1 + tao) +
		(EE0*TP_const / (1.0 + TP_const)) * qq_now_single * cos(Omega2 * t1 + 2 * tao);

}

__host__  double compute_ee0_by_index(const int index)
{
	double EE0 = 2.742*pow(10, 3)*sqrt(pow(10.0, (12.0 + double(index)*0.2)));
	EE0 = EE0 / (5.1421*(pow(10.0, 11.0)));
	return  EE0;
}