
#include "../include/Laser.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Runge_Kutta.cuh"

#include <cmath>





__global__ void pre_second_step_E_forcheck(const double* E1, const double* E2, double* E_check)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		E_check[idx] = sqrt(pow(E1[idx], 2) + pow(E2[idx], 2));
	}
}








__global__ void pre_second_step_e1(const double* QQ, const double EE0, double* E1)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		double t1 = 0.5 * DX * idx;
		/*curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;*/
		double tao = 0.0;
		E1[idx] = (EE0 / (1.0 + TP_const)) * QQ[idx] * sin(Omega1 * t1 + tao) -
			(EE0*TP_const / (1.0 + TP_const)) * QQ[idx] * sin(Omega2 * t1 + 2 * tao);

	}

}


__global__ void pre_second_step_e2(const double* QQ, const double EE0, double* E2)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		double t1 = 0.5 * DX * idx;
		/*	curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;*/
		double tao = 0.0;

		E2[idx] = (EE0 / (1.0 + TP_const)) * QQ[idx] * cos(Omega1 * t1 + tao) +
			(EE0*TP_const / (1.0 + TP_const)) * QQ[idx] * cos(Omega2 * t1 + 2 * tao);

	}

}
















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

