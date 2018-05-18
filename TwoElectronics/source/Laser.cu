
#include "../include/Laser.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Runge_Kutta.cuh"

#include <cmath>



__device__ double CalculationE1(const nucleus& first, const nucleus& second)
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

__device__ double CalculationE2(const nucleus& first, const nucleus& second)
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

