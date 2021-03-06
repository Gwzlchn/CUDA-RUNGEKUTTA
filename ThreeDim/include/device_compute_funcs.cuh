﻿#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH

#include <curand.h>
#include <curand_kernel.h>
#include <cmath>


#include "sm_20_atomic_functions.h"
#include "../include/sci_const.h"
#include "../include/nucleus.hpp"


struct derivative
{
	double px;
	double py;
	double pz;
	double fx;
	double fy;
	double fz;
};

//初始分配函数
//计算总动能
__device__ double E_kall(const nucleus& first, const nucleus& second);

//Px Py Pz分配
__device__ void px_py_pz_distribution(nucleus& first, nucleus& second, double ekall, int i);



//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);

//第一个核，三个坐标的二阶导
__device__ double3 fx_first_nucleus(const nucleus& first, const nucleus& second);

//第二个核，三个坐标的二阶导
__device__ double3 fx_second_nucleus(const nucleus& first, const nucleus& second);

//龙哥库塔方法
__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second);
__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second,
	const double& e_laser_t1, const double& e_laser_t2,
	const double& e_laser_t3, const double& e_laser_t4);



__device__ derivative fisrt_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second);
__device__ derivative second_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second);
__device__ nucleus first_and_second_k_add_dx_raw(const derivative& k_one_to_four, const nucleus& raw_nucleus);
__device__ nucleus first_and_second_k_add_dx_div(const derivative& k_one_to_four, const nucleus& raw_nucleus);
__device__ void k_one_to_four_add(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4,
	nucleus& raw_nucleus);


//第二步runge-kutta
__device__ derivative fisrt_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e_laser);

__device__ derivative second_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e_laser);

//E1与E2
__device__ double CalculationE1	(const nucleus& first, const nucleus& second);
__device__ double CalculationE2	(const nucleus& first, const nucleus& second);
//__device__ void count_ee1_and_ee2(const nucleus& first, const nucleus& second, unsigned long long* e_laser);
#endif //DEVICE_COMPUTE_FUNCS_CUH








__device__ double E_kall(const nucleus& first, const nucleus& second)
{
	return E_total - (-1.0 / sqrt(pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y +elec_nucl * elec_nucl)) -
					(-1.0 / sqrt(pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + elec_nucl * elec_nucl))
					- (1.0 / sqrt(nucleus_distance(first, second) + elec_elec*elec_elec)) -
					(-1.0 / sqrt(pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl)) -
					(-1.0 / sqrt(pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + elec_nucl * elec_nucl));
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

	double theta1 = random1 * 2 * PI;
	double theta2 = random2 * 2 * PI;
	double phi1 = random3 * PI;
	double phi2 = random4 * PI;

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




__device__  double3 fx_first_nucleus(const nucleus& first, const nucleus& second)
{
	double3 fx_fy_fz_first;

	fx_fy_fz_first.x = (first.x - second.x)
					/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
				- (first.x - nuclear_spacing / 2.0 * sin(PI*rotation))
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3))
				- (first.x + nuclear_spacing / 2.0 * sin(PI*rotation))
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3));

	fx_fy_fz_first.y = (first.y - second.y)
					/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
				- first.y
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3))
				- first.y
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3));

	fx_fy_fz_first.z = (first.z - second.z)
					/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
				- (first.z - nuclear_spacing / 2.0 * cos(PI*rotation))
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3))
				- (first.z + nuclear_spacing / 2.0 * cos(PI*rotation))
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + elec_nucl * elec_nucl), 3));

	
	return fx_fy_fz_first;
}

__device__ double3 fx_second_nucleus(const nucleus& first, const nucleus& second)
{
	double3 fx_fy_fz_second;
	
	fx_fy_fz_second.x = (second.x - first.x)
							/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
						- (second.x - nuclear_spacing / 2.0 * sin(PI*rotation))
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3))
						- (second.x + nuclear_spacing / 2.0 * sin(PI*rotation))
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3));

	fx_fy_fz_second.y = (second.y - first.y)
							/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
						- second.y
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3))
						- second.y
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3));

	fx_fy_fz_second.z = ( second.z- first.z )
							/ sqrt(pow((nucleus_distance(first, second) + elec_elec * elec_elec), 3))
						- (second.z - nuclear_spacing / 2.0 * cos(PI*rotation))
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3))
						- (second.z + nuclear_spacing / 2.0 * cos(PI*rotation))
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + elec_nucl*elec_nucl), 3));
	

	return fx_fy_fz_second;
}




//在pair_k里面 x为二阶导！px为原始px；
////K1=(/g1(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1),&
//&f1(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1), &
//&g2(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1), &
//&f2(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1), &
//&g3(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1), &
//&f3(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1), &
//&g4(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1), &
//&f4(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1), &
//&g5(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1), &
//&f5(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1) - EE, &
//&g6(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1), &
//&f6(y11, y12, y13, y14, y15, y16, y17, y18, y111, y112, y113, y114, t1, a, a1) - EE / )
__device__ derivative fisrt_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second)
{
	double3 first_fx = fx_first_nucleus(first, second);
	derivative first_px_fx;
	first_px_fx.px = first.px;
	first_px_fx.py = first.py;
	first_px_fx.pz = first.pz;
	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y;
	first_px_fx.fz = first_fx.z;

	return first_px_fx;
	
}


__device__ derivative second_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second)
{
	double3 second_fx = fx_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second.px;
	second_px_fx.py = second.py;
	second_px_fx.pz = second.pz;
	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y;
	second_px_fx.fz = second_fx.z;
	return second_px_fx;

}


__device__ void k_one_to_four_add(const derivative& K1,const derivative& K2,const derivative& K3,const derivative& K4,
				nucleus& raw_nucleus)
{
	raw_nucleus.x = raw_nucleus.x + DX * (K1.px + 2.0*K2.px + 2.0*K3.px + K4.px) / 6.0;
	raw_nucleus.y = raw_nucleus.y + DX * (K1.py + 2.0*K2.py + 2.0*K3.py + K4.py) / 6.0;
	raw_nucleus.z = raw_nucleus.z + DX * (K1.pz + 2.0*K2.pz + 2.0*K3.pz + K4.pz) / 6.0;
	raw_nucleus.px = raw_nucleus.px + DX * (K1.fx + 2.0*K2.fx + 2.0*K3.fx + K4.fx) / 6.0;
	raw_nucleus.py = raw_nucleus.py + DX * (K1.fy + 2.0*K2.fy + 2.0*K3.fy + K4.fy) / 6.0;
	raw_nucleus.pz = raw_nucleus.pz + DX * (K1.fz + 2.0*K2.fz + 2.0*K3.fz + K4.fz) / 6.0;



	return;
}


 nucleus first_and_second_k_add_dx_raw(const derivative& k_one_to_four,const nucleus& raw_nucleus)
{
	double now_dx = DX;

	nucleus k_add;
	k_add.x = raw_nucleus.x + now_dx * k_one_to_four.px ;
	k_add.y = raw_nucleus.y + now_dx * k_one_to_four.py ;
	k_add.z = raw_nucleus.z + now_dx * k_one_to_four.pz ;
	k_add.px = raw_nucleus.px + now_dx * k_one_to_four.fx ;
	k_add.py = raw_nucleus.py + now_dx * k_one_to_four.fy ;
	k_add.pz = raw_nucleus.pz + now_dx * k_one_to_four.fz ;

	return k_add;
}

__device__ nucleus first_and_second_k_add_dx_div(const derivative& k_one_to_four, const nucleus& raw_nucleus)
{
	double now_dx = DX/2.0;

	nucleus k_add;
	k_add.x = raw_nucleus.x + now_dx * k_one_to_four.px;
	k_add.y = raw_nucleus.y + now_dx * k_one_to_four.py;
	k_add.z = raw_nucleus.z + now_dx * k_one_to_four.pz;
	k_add.px = raw_nucleus.px + now_dx * k_one_to_four.fx;
	k_add.py = raw_nucleus.py + now_dx * k_one_to_four.fy;
	k_add.pz = raw_nucleus.pz + now_dx * k_one_to_four.fz;

	return k_add;
}
   



__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second)
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	const derivative second_k1 = second_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	const nucleus first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	const nucleus second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	const derivative second_k2 = second_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	const nucleus first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	const nucleus second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	const derivative second_k3 = second_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	const nucleus first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	const nucleus second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);
	const derivative second_k4 = second_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);

	k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);

	return;
}



__device__ derivative fisrt_k_one_to_four_second_step
(const nucleus& first, const nucleus& second,const double& e_laser)
{
	const double3 first_fx = fx_first_nucleus(first, second);
	derivative first_px_fx;
	first_px_fx.px = first.px;
	first_px_fx.py = first.py;
	first_px_fx.pz = first.pz;
	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y;
	first_px_fx.fz = first_fx.z + e_laser;

	return first_px_fx;

}


__device__ derivative second_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e_laser)
{
	const double3 second_fx = fx_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second.px;
	second_px_fx.py = second.py;
	second_px_fx.pz = second.pz;
	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y;
	second_px_fx.fz = second_fx.z + e_laser;
	return second_px_fx;

}







__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second,
	const double& e_laser_t1,const double& e_laser_t2, 
	const double& e_laser_t3, const double& e_laser_t4 )
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four_second_step(step_one_first, step_one_second,e_laser_t1);
	const derivative second_k1 = second_k_one_to_four_second_step(step_one_first, step_one_second, e_laser_t1);
	const nucleus first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	const nucleus second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four_second_step(first_k1_add, second_k1_add, e_laser_t2);
	const derivative second_k2 = second_k_one_to_four_second_step(first_k1_add, second_k1_add, e_laser_t2);
	const nucleus first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	const nucleus second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four_second_step(first_k2_add, second_k2_add,e_laser_t3);
	const derivative second_k3 = second_k_one_to_four_second_step(first_k2_add, second_k2_add, e_laser_t3);
	const nucleus first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	const nucleus second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four_second_step(first_k3_add, second_k3_add, e_laser_t4);
	const derivative second_k4 = second_k_one_to_four_second_step(first_k3_add, second_k3_add, e_laser_t4);

	k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);
}





__device__ double CalculationE1(const nucleus& first, const nucleus& second)
{
	return  0.5*(pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2)) +
			(-1.0 / sqrt(pow(first.z + (nuclear_spacing / 2.0)*cos(PI*rotation), 2) +
						pow(first.x + (nuclear_spacing / 2.0)*sin(PI*rotation), 2) + 
						pow(first.y, 2) + pow(elec_nucl, 2))) +
			(-1.0 / sqrt(pow(first.z - (nuclear_spacing / 2.0)*cos(PI*rotation), 2) +
						pow(first.x - (nuclear_spacing / 2.0)*sin(PI*rotation), 2) +
						pow(first.y, 2) + pow(elec_nucl, 2))) +
			0.5*(1.0 / sqrt(pow(first.x - second.x, 2) + pow(first.y - second.y, 2) +
						pow(first.z - second.z, 2) + pow(elec_elec, 2)));
}

__device__ double CalculationE2(const nucleus& first, const nucleus& second)
{
	return 0.5*(pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2)) +
			(-1 / sqrt(pow(second.z + (nuclear_spacing / 2)*cos(PI*rotation), 2) +
						pow(second.x + (nuclear_spacing / 2)*sin(PI*rotation), 2) +
						pow(second.y, 2) + pow(elec_nucl, 2))) +
			(-1 / sqrt(pow(second.z - (nuclear_spacing / 2)*cos(PI*rotation), 2) +
						pow(second.x - (nuclear_spacing / 2)*sin(PI*rotation), 2) + 
						pow(second.y, 2) + pow(elec_nucl, 2))) +
			0.5*(1 / sqrt(pow(first.x - second.x, 2) + pow(first.y - second.y, 2) + 
						pow(first.z - second.z, 2) + pow(elec_elec, 2)));
}


