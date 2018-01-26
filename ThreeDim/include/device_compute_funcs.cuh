#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH

#include <curand.h>
#include <curand_kernel.h>
#include <cmath>
#include <vector_types.h>

#include "../include/sci_const.h"
//#include "../include/device_compute_funcs.cuh"
#include "../include/nucleus.hpp"
//#include <crt/host_defines.h>

//计算总动能
__device__ double E_kall(const nucleus& first, const nucleus& second);

//Px Py Pz
__device__ void px_py_pz_distribution(nucleus& first, nucleus& second, double ekall, int i);

//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);

//第一个核，三个坐标的二阶导
__device__ nucleus fx_first_nucleus(const nucleus& first, const nucleus& second);

//第二个核，三个坐标的二阶导
__device__ nucleus fx_second_nucleus(const nucleus& first, const nucleus& second);

//用于双核粒子的随机数化
__device__ void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec);
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec);

#endif //DEVICE_COMPUTE_FUNCS_CUH

__device__ double E_kall(const nucleus& first, const nucleus& second)
{
	return E_total - (-1.0 / sqrt(pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing - 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + elec_elec*elec_elec)) -
					(-1.0 / sqrt(pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + elec_elec*elec_elec))
					- (1.0 / sqrt(nucleus_distance(first, second) + elec_nucl*elec_nucl)) -
					(-1.0 / sqrt(pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + elec_elec*elec_elec)) -
					(-1.0 / sqrt(pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + elec_elec*elec_elec));
}

__device__ void px_py_pz_distribution(nucleus& first, nucleus& second,double ekall,int i)
{



	//double ekall = E_kall(first, second);
	curandStatePhilox4_32_10_t s;
	//unsigned long long seed = threadIdx.x;
	// seed a random number generator 
	curand_init(i, 0, 0, &s);

	double random1 = curand_uniform_double(&s);
	double random2 = curand_uniform_double(&s);
	double random3 = curand_uniform_double(&s);
	double random4 = curand_uniform_double(&s);
	double random5 = curand_uniform_double(&s);

	double theta1 = random1 * PI;
	double theta2 = random2 * PI;
	double phi1 = random3 * 2 * PI;
	double phi2 = random4 * 2 * PI;

	first.px = sqrt(2.0*ekall*random5)*cos(theta1)*sin(phi1);
	first.py = sqrt(2.0*ekall*random5)*sin(theta1)*sin(phi1);
	first.pz = sqrt(2.0*ekall*random5)*cos(phi1);

	second.px = sqrt(2.0*ekall*(1 - random5))*cos(theta2)*sin(phi2);
	second.py = sqrt(2.0*ekall*(1 - random5))*sin(theta2)*sin(phi2);
	second.pz = sqrt(2.0*ekall*(1 - random5))*cos(phi2);
}


__device__  double nucleus_distance(const nucleus& first, const nucleus& second)
{
	return (pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2));
}




__device__  nucleus fx_first_nucleus(const nucleus& first, const nucleus& second)
{
	nucleus fx_first;

	fx_first.x = (first.x - second.x)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (first.x - nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- (first.x + nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.y = (first.y - second.y)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- first.y
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- first.y
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 3) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.z = (first.z - second.z)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (first.z - nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- (first.z + nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.px = fx_first.py = fx_first.pz = 0;
	return fx_first;
}

__device__  nucleus fx_second_nucleus(const nucleus& first, const nucleus& second)
{
	nucleus fx_second;
	fx_second.x = (second.x - first.x)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (second.x - nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- (second.x + nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3));

	fx_second.y = (second.y - first.y)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- second.y
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- second.y
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 3) +
			second.y*second.y + elec_elec*elec_elec), 3));

	fx_second.z = (first.z - second.z)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (second.z - nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- (second.z + nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3));
	fx_second.px = fx_second.py = fx_second.pz = 0;
	return fx_second;
}

