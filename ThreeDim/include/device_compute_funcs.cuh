#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH

#include <curand.h>
#include <curand_kernel.h>
#include <cmath>


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
__device__ void update_step_two(nucleus& step_two_first, nucleus& step_two_second);


__device__ derivative fisrt_k_one_to_four(const nucleus& first, const nucleus& second);
__device__ derivative second_k_one_to_four(const nucleus& first, const nucleus& second);
__device__ nucleus first_and_second_k_add(const derivative& k_one_to_four, const nucleus& raw_nucleus, int choose);
__device__ void k_one_to_four_add(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4,
	nucleus& raw_nucleus);

#endif //DEVICE_COMPUTE_FUNCS_CUH








__device__ double E_kall(const nucleus& first, const nucleus& second)
{
	const double power_e_n = elec_nucl * elec_nucl;
	return E_total - (-1.0 / sqrt(pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + power_e_n)) -
					(-1.0 / sqrt(pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + power_e_n))
					- (1.0 / sqrt(nucleus_distance(first, second) + elec_nucl*elec_nucl)) -
					(-1.0 / sqrt(pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + power_e_n)) -
					(-1.0 / sqrt(pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								second.y*second.y + power_e_n));
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
	const double pow_elec_elec = elec_elec * elec_elec; // A1*A1
	const double pow_elec_nucl = elec_nucl * elec_nucl; // A*A

	fx_fy_fz_first.x = (first.x - second.x)
					/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
				- (first.x - nuclear_spacing / 2.0 * sin(PI*rotation))
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3))
				- (first.x + nuclear_spacing / 2.0 * sin(PI*rotation))
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3));

	fx_fy_fz_first.y = (first.y - second.y)
					/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
				- first.y
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3))
				- first.y
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3));

	fx_fy_fz_first.z = (first.z - second.z)
					/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
				- (first.z - nuclear_spacing / 2.0 * cos(PI*rotation))
					/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3))
				- (first.z + nuclear_spacing / 2.0 * cos(PI*rotation))
					/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
								pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
								first.y*first.y + pow_elec_nucl), 3));

	
	return fx_fy_fz_first;
}

__device__ double3 fx_second_nucleus(const nucleus& first, const nucleus& second)
{
	double3 fx_fy_fz_second;
	const double pow_elec_elec = elec_elec * elec_elec; // A1*A1
	const double pow_elec_nucl = elec_nucl * elec_nucl; // A*A
	fx_fy_fz_second.x = (second.x - first.x)
							/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
						- (second.x - nuclear_spacing / 2.0 * sin(PI*rotation))
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3))
						- (second.x + nuclear_spacing / 2.0 * sin(PI*rotation))
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3));

	fx_fy_fz_second.y = (second.y - first.y)
							/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
						- second.y
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3))
						- second.y
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3));

	fx_fy_fz_second.z = ( second.z- first.z )
							/ sqrt(pow((nucleus_distance(first, second) + pow_elec_elec), 3))
						- (second.z - nuclear_spacing / 2.0 * cos(PI*rotation))
							/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3))
						- (second.z + nuclear_spacing / 2.0 * cos(PI*rotation))
							/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
										pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
										second.y*second.y + pow_elec_nucl), 3));
	

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
__device__ derivative fisrt_k_one_to_four(const nucleus& first, const nucleus& second)
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


__device__ derivative second_k_one_to_four(const nucleus& first, const nucleus& second)
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


__device__ nucleus first_and_second_k_add(const derivative& k_one_to_four,const nucleus& raw_nucleus,int choose)
{
	double now_dx = DX;
	if (choose == 2)
		now_dx = DX / 2.0;

	nucleus k_add;
	k_add.x = raw_nucleus.x + now_dx * k_one_to_four.px ;
	k_add.y = raw_nucleus.y + now_dx * k_one_to_four.py ;
	k_add.z = raw_nucleus.z + now_dx * k_one_to_four.pz ;
	k_add.px = raw_nucleus.px + now_dx * k_one_to_four.fx ;
	k_add.py = raw_nucleus.py + now_dx * k_one_to_four.fy ;
	k_add.pz = raw_nucleus.pz + now_dx * k_one_to_four.fz ;

	return k_add;
}




__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second)
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four(step_one_first, step_one_second);
	const derivative second_k1 = second_k_one_to_four(step_one_first, step_one_second);
	const nucleus first_k1_add = first_and_second_k_add(first_k1, step_one_first, 2);
	const nucleus second_k1_add = first_and_second_k_add(second_k1, step_one_second, 2);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four(first_k1_add, second_k1_add);
	const derivative second_k2 = second_k_one_to_four(first_k1_add, second_k1_add);
	const nucleus first_k2_add = first_and_second_k_add(first_k2, step_one_first, 2);
	const nucleus second_k2_add = first_and_second_k_add(second_k2, step_one_second, 2);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four(first_k2_add, second_k2_add);
	const derivative second_k3 = second_k_one_to_four(first_k2_add, second_k2_add);
	const nucleus first_k3_add = first_and_second_k_add(first_k3, step_one_first, 1);
	const nucleus second_k3_add = first_and_second_k_add(second_k3, step_one_second, 1);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four(first_k3_add, second_k3_add);
	const derivative second_k4 = second_k_one_to_four(first_k3_add, second_k3_add);

	k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);

	return;
}

__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second)
{

}

