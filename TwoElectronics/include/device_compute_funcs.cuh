#ifndef DEVICE_COMPUTE_FUNCS_CUH
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

//--------------------------初始分配函数
//初始矩阵，取最小值

//计算总动能
__device__ double E_kall(const nucleus& first, const nucleus& second);


//获得六个随机数，如果有1 则置0 cuda 范围是（0，1] 但fortran随机数范围是[0, 1)
__device__ void get_six_random(double2& two_random, double4& four_random,const int& seed);

// x y z Px Py Pz分配 seed 种子，min_r 对应 rr min_p 对应 pp
__device__ void distribution(nucleus& first, nucleus& second, 
							 const int& seed,const double& min_r,const double& min_p);



//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);


//第一个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_first_nucleus(const nucleus& first, const nucleus& second);

//第二个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_second_nucleus(const nucleus& first, const nucleus& second);

//第一个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_first_nucleus(const nucleus& first, const nucleus& second);

//第二个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_second_nucleus(const nucleus& first, const nucleus& second);

//龙哥库塔方法
__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second);
__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second,
	const double4 e1_laser_now,const double4 e2_laser_now);


//第一个粒子 K1~K4 第一步循环
__device__ derivative fisrt_k_one_to_four_fisrt_step
(const nucleus& first, const nucleus& second);
//第二个粒子 K1~K4 第一步循环
__device__ derivative second_k_one_to_four_fisrt_step
(const nucleus& first, const nucleus& second);

//计算完 K3 更新下一步参数用，K3 不除以2
__device__ nucleus first_and_second_k_add_dx_raw
(const derivative& k_one_to_four, const nucleus& raw_to_add);
//计算完 K1或 K2 更新下一步参数用，K1 或 K2 除以2
__device__ nucleus first_and_second_k_add_dx_div
(const derivative& k_one_to_four, const nucleus& raw_to_add);
//计算完 K1~K4 完整的相加 
__device__ void k_one_to_four_add
(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4, nucleus& raw_to_add);


//第二步runge-kutta
__device__ derivative fisrt_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e1_laser, const double& e2_laser);


__device__ derivative second_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e1_laser, const double& e2_laser);

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

__device__ void get_six_random(double2& two_random, double4& four_random,const int& seed)
{
	curandStatePhilox4_32_10_t s;
	curand_init(seed, 0, 0, &s);
	two_random = curand_uniform2_double(&s);
	four_random = make_double4(curand_uniform_double(&s), curand_uniform_double(&s), 
				curand_uniform_double(&s), curand_uniform_double(&s) );
	if (two_random.x == 1)
		two_random.x = 0;
	if (two_random.y == 1)
		two_random.y = 0;
	
	if (four_random.x == 1)
		four_random.x = 0;
	if (four_random.y == 1)
		four_random.y = 0;
	if (four_random.z == 1)
		four_random.z = 0;
	if (four_random.w == 1)
		four_random.w= 0;
	

}


__device__ void distribution(nucleus& first, nucleus& second, 
					const int& seed, const double& min_r, const double& min_p)
{

	double2 two_rand;
	double4 four_rand;
	get_six_random(two_rand, four_rand,seed);

	double theta1 = two_rand.x * 2.0 * PI;
	double phi1 = two_rand.y * PI;
	

	first.x = min_r * sin(phi1) * cos(theta1);
	first.y = min_r * sin(phi1) * sin(theta1);
	first.z = min_r * cos(phi1);

	second.x = -first.x;
	second.y = -first.y;
	second.z = -first.z;
	
	
	double phi2 = four_rand.x * PI;
	double phi3 = four_rand.y * PI;
	double theta2 = four_rand.z * 2.0 * PI;
	double theta3 = four_rand.w * 2.0 * PI;
	
	first.px = min_p * cos(theta2)*sin(phi2);
	first.py = min_p * sin(theta2)*sin(phi2);
	first.pz = min_p * cos(phi2);

	second.px = min_p * cos(theta3)*sin(phi3);
	second.py = min_p * sin(theta3)*sin(phi3);
	second.pz = min_p * cos(phi3);
}


__device__  double nucleus_distance(const nucleus& first, const nucleus& second)
{
	return (pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2));
}


//第一个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_first_nucleus(const nucleus& first, const nucleus& second)
{
	double Q_squre = pow(Q_constant, 2);
	//坐标平方和
	double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	
	//一阶导平方和
	double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);



	//第一个核 一阶导三个计算公式 对应 g1 g3 g5
	double3 gx_gy_gz;
	gx_gy_gz.x = first.px * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first/ Q_squre, 2))));

	gx_gy_gz.y = first.py * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));
	
	gx_gy_gz.z = first.pz * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));
	
	return gx_gy_gz;
}

//第二个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_second_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);
	
	//坐标平方和
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);

	//一阶导平方和
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);



	//第二个核 一阶导三个计算公式 对应 g2 g4 g6
	double3 gx_gy_gz;
	gx_gy_gz.x = second.px * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	gx_gy_gz.y = second.py * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	gx_gy_gz.z = second.pz * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	return gx_gy_gz;
}

//第一个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_first_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);

	//坐标平方和
	const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);
	//两核距离平方的1.5 次方,对应 sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
	const double distance_1_5_power = pow(distance_squre, 1.5);

	//坐标平方和的1.5 次方 ，对应 sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0)
	const double loc_1_5_power_first = pow(loc_squre_sum_first, 1.5);
	const double loc_1_5_power_second = pow(loc_squre_sum_second, 1.5);

	//临时变量1 ,第一个粒子（一阶导平方和的平方 / Q方）
	const double temp1 = pow(px_py_pz_squre_sum_first, 2) / Q_squre;
	//临时变量2，第一个粒子（坐标平方和 * 一阶导平方和 / Q方）的平方
	const double temp2 = pow((loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre), 2);


	double3 fx_fy_fz;
	fx_fy_fz.x = first.x * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.x - second.x) / distance_1_5_power;

	fx_fy_fz.y = first.y * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.y - second.y) / distance_1_5_power;

	fx_fy_fz.z = first.z * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.z - second.z) / distance_1_5_power;

	return  fx_fy_fz;
}


//第二个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_second_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);

	//坐标平方和
	const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);
	//两核距离平方的1.5 次方,对应 sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
	const double distance_1_5_power = pow(distance_squre, 1.5);

	//坐标平方和的1.5 次方 ，对应 sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0)
	const double loc_1_5_power_first = pow(loc_squre_sum_first, 1.5);
	const double loc_1_5_power_second = pow(loc_squre_sum_second, 1.5);

	//临时变量1 ,第二个粒子（一阶导平方和的平方 / Q方）
	//对应(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0
	const double temp1 = pow(px_py_pz_squre_sum_second, 2) / Q_squre;
	//临时变量2，第二个粒子（坐标平方和 * 一阶导平方和 / Q方）的平方
	const double temp2 = pow((loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre), 2);


	double3 fx_fy_fz;
	fx_fy_fz.x = second.x * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.x - second.x) / distance_1_5_power;

	fx_fy_fz.y = second.y * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.y - second.y) / distance_1_5_power;

	fx_fy_fz.z = second.z * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.z - second.z) / distance_1_5_power;

	return  fx_fy_fz;

}










//第一个粒子 K1~K4 第一步循环
__device__ derivative fisrt_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second)
{
	//二阶导 三个数
	const double3 first_fx = fx_fy_fz_first_nucleus(first, second);
	//一阶导 三个数
	const double3 first_gx = gx_gy_gz_first_nucleus(first, second);
	
	derivative first_px_fx;
	first_px_fx.px = first_gx.x;
	first_px_fx.py = first_gx.y;
	first_px_fx.pz = first_gx.z;

	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y;
	first_px_fx.fz = first_fx.z;

	return first_px_fx;
	
}


__device__ derivative second_k_one_to_four_fisrt_step(const nucleus& first, const nucleus& second)
{
	//二阶导 三个数
	double3 second_fx = fx_fy_fz_second_nucleus(first, second);
	//一阶导 三个数
	double3 second_gx = gx_gy_gz_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second_gx.x;
	second_px_fx.py = second_gx.y;
	second_px_fx.pz = second_gx.z;
	
	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y;
	second_px_fx.fz = second_fx.z;
	return second_px_fx;

}



__device__ nucleus first_and_second_k_add_dx_raw(const derivative& k_one_to_four, const nucleus& raw_to_add)
{
	double now_dx = DX;

	nucleus k_add;
	k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	return k_add;
}

__device__ nucleus first_and_second_k_add_dx_div(const derivative& k_one_to_four, const nucleus& raw_to_add)
{
	double now_dx = DX / 2.0;

	nucleus k_add;
	k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	return k_add;
}






__device__ void k_one_to_four_add(const derivative& K1,const derivative& K2,const derivative& K3,const derivative& K4,
				nucleus& raw_to_add)
{
	raw_to_add.x = raw_to_add.x + DX * (K1.px + 2.0*K2.px + 2.0*K3.px + K4.px) / 6.0;
	raw_to_add.y = raw_to_add.y + DX * (K1.py + 2.0*K2.py + 2.0*K3.py + K4.py) / 6.0;
	raw_to_add.z = raw_to_add.z + DX * (K1.pz + 2.0*K2.pz + 2.0*K3.pz + K4.pz) / 6.0;
	raw_to_add.px = raw_to_add.px + DX * (K1.fx + 2.0*K2.fx + 2.0*K3.fx + K4.fx) / 6.0;
	raw_to_add.py = raw_to_add.py + DX * (K1.fy + 2.0*K2.fy + 2.0*K3.fy + K4.fy) / 6.0;
	raw_to_add.pz = raw_to_add.pz + DX * (K1.fz + 2.0*K2.fz + 2.0*K3.fz + K4.fz) / 6.0;



	return;
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





//第一个粒子 K1~K4 第二步循环
__device__ derivative fisrt_k_one_to_four_second_step
(const nucleus& first, const nucleus& second,const double& e1_laser,const double& e2_laser)
{
	const double3 first_fx =  fx_fy_fz_first_nucleus(first, second);
	const double3 first_gx  = gx_gy_gz_first_nucleus(first, second);
	derivative first_px_fx;
	first_px_fx.px = first_gx.x;
	first_px_fx.py = first_gx.y;
	first_px_fx.pz = first_gx.z;
	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y - e2_laser;
	first_px_fx.fz = first_fx.z - e1_laser;

	return first_px_fx;

}

//第二个粒子 K1~K4 第二步循环
__device__ derivative second_k_one_to_four_second_step
(const nucleus& first, const nucleus& second, const double& e1_laser, const double& e2_laser)
{
	const double3 second_fx = fx_fy_fz_second_nucleus(first, second);
	const double3 second_gx = gx_gy_gz_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second.px;
	second_px_fx.py = second.py;
	second_px_fx.pz = second.pz;
	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y - e2_laser;
	second_px_fx.fz = second_fx.z - e1_laser;
	return second_px_fx;

}







__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second,
	const double4 e1_laser_now, const double4 e2_laser_now)
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four_second_step(step_one_first, step_one_second,e1_laser_now.x,e2_laser_now.x);
	const derivative second_k1 = second_k_one_to_four_second_step(step_one_first, step_one_second, e1_laser_now.x, e2_laser_now.x);
	const nucleus first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	const nucleus second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	const derivative second_k2 = second_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	const nucleus first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	const nucleus second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	const derivative second_k3 = second_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	const nucleus first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	const nucleus second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);
	const derivative second_k4 = second_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);

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


