#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH

#include <curand.h>
#include <curand_kernel.h>
#include <cmath>


#include "../include/sci_const.h"
#include "../include/nucleus.hpp"


//初始分配函数
//计算总动能
__device__ double E_kall(const nucleus& first, const nucleus& second);

//Px Py Pz分配
__device__ void px_py_pz_distribution(nucleus& first, nucleus& second, double ekall, int i);



//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);

//第一个核，三个坐标的二阶导
__device__ nucleus fx_first_nucleus(const nucleus& first, const nucleus& second);

//第二个核，三个坐标的二阶导
__device__ nucleus fx_second_nucleus(const nucleus& first, const nucleus& second);

//龙哥库塔方法
__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second);
__device__ void update_step_two(nucleus& step_two_first, nucleus& step_two_second);



__device__ nuclei pair_add_dx(const nucleus& first, const nucleus& second, nuclei k_one_to_four, int choose);



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



//choose = 2 ：balabalba+dx/2 
//=1: something+ dx 
//代替这个
//y21 = x1(1, j) + h*K1(1) / 2
//y22 = px1(1, j) + h*K1(2) / 2
//y23 = x2(1, j) + h*K1(3) / 2
//y24 = px2(1, j) + h*K1(4) / 2
//y25 = y1(1, j) + h*K1(5) / 2
//y26 = Py1(1, j) + h*K1(6) / 2
//y27 = y2(1, j) + h*K1(7) / 2
//y28 = Py2(1, j) + h*K1(8) / 2
//y211 = z1(1, j) + h*K1(9) / 2
//y212 = pz1(1, j) + h*K1(10) / 2
//y213 = z2(1, j) + h*K1(11) / 2
//y214 = pz2(1, j) + h*K1(12) / 2
__device__ nuclei pair_add_dx(const nucleus& first, const nucleus& second,nuclei k_one_to_four,int choose)
{
	double h;
	if (choose == 1)
		h = dx;
	else if (choose == 2)
		h = dx / 2.0;
	
	nuclei add_dx;
	add_dx.first.x = first.x + h*k_one_to_four.first.x;
	add_dx.first.y = first.y + h*k_one_to_four.first.y;
	add_dx.first.z = first.z + h*k_one_to_four.first.z;
	add_dx.first.px = first.px + h*k_one_to_four.first.px;
	add_dx.first.py = first.py + h*k_one_to_four.first.py;
	add_dx.first.pz = first.pz + h*k_one_to_four.first.pz;

	add_dx.second.x = second.x + h*k_one_to_four.second.x;
	add_dx.second.y = second.y + h*k_one_to_four.second.y;
	add_dx.second.z = second.z + h*k_one_to_four.second.z;
	add_dx.second.px = second.px + h*k_one_to_four.second.px;
	add_dx.second.py = second.py + h*k_one_to_four.second.py;
	add_dx.second.pz = second.pz + h*k_one_to_four.second.pz;

	return add_dx;

}

__device__ nuclei pair_k(const nucleus& first, const nucleus& second)
{
	nucleus fir_k1 = fx_first_nucleus(first, second);
	fir_k1.px = first.px;
	fir_k1.py = first.py;
	fir_k1.pz = first.pz;
	nucleus sec_k1 = fx_second_nucleus(first, second);
	sec_k1.px = second.px;
	sec_k1.py = second.py;
	sec_k1.pz = second.pz;

	return { fir_k1, sec_k1 };
	
}


//__device__ nuclei pair_k2(const nucleus& first, const nucleus& second)
//{
//	//nuclei k2 = pair_add_dx(first,second)
//	nucleus fir_k1 = fx_first_nucleus(first, second);
//	fir_k1.px = first.px;
//	fir_k1.py = first.py;
//	fir_k1.pz = first.pz;
//	nucleus sec_k1 = fx_second_nucleus(first, second);
//	sec_k1.px = second.px;
//	sec_k1.py = second.py;
//	sec_k1.pz = second.pz;
//	return { fir_k2, sec_k2 };
//}
__device__ void k_one_to_four_add(const nuclei& k1,const nuclei& k2,const nuclei& k3,const nuclei& k4,
				nucleus& raw_first,nucleus& raw_second)
{
	nuclei k_total;
	k_total.first.x = 1 / 6 * (k1.first.x + 2 * k2.first.x + 2 * k3.first.x + k4.first.x);
	k_total.first.y = 1 / 6 * (k1.first.y + 2 * k2.first.y + 2 * k3.first.y + k4.first.y);
	k_total.first.z = 1 / 6 * (k1.first.z + 2 * k2.first.z + 2 * k3.first.z + k4.first.z);
	k_total.first.px = 1 / 6 * (k1.first.px + 2 * k2.first.px + 2 * k3.first.px + k4.first.px);
	k_total.first.py = 1 / 6 * (k1.first.py + 2 * k2.first.py + 2 * k3.first.py + k4.first.py);
	k_total.first.pz = 1 / 6 * (k1.first.pz + 2 * k2.first.pz + 2 * k3.first.pz + k4.first.pz);
	k_total.second.x = 1 / 6 * (k1.second.x + 2 * k2.second.x + 2 * k3.second.x + k4.second.x);
	k_total.second.y = 1 / 6 * (k1.second.y + 2 * k2.second.y + 2 * k3.second.y + k4.second.y);
	k_total.second.z = 1 / 6 * (k1.second.z + 2 * k2.second.z + 2 * k3.second.z + k4.second.z);
	k_total.second.px = 1 / 6 * (k1.second.px + 2 * k2.second.px + 2 * k3.second.px + k4.second.px);
	k_total.second.py = 1 / 6 * (k1.second.py + 2 * k2.second.py + 2 * k3.second.py + k4.second.py);
	k_total.second.pz = 1 / 6 * (k1.second.pz + 2 * k2.second.pz + 2 * k3.second.pz + k4.second.pz);

	raw_first.x += k_total.first.x;
	raw_first.y += k_total.first.y;
	raw_first.z += k_total.first.z;
	raw_first.px += k_total.first.px;
	raw_first.py += k_total.first.py;
	raw_first.pz += k_total.first.pz;

	raw_second.x += k_total.second.x;
	raw_second.y += k_total.second.y;
	raw_second.z += k_total.second.z;
	raw_second.px += k_total.second.px;
	raw_second.py += k_total.second.py;
	raw_second.pz += k_total.second.pz;

	return;
}


__device__ void update_step_one(nucleus& step_one_first, nucleus& step_one_second)
{
	/*nucleus x_y_z_fir_fx1 = fx_first_nucleus(step_one_first, step_one_second);
	nucleus x_y_z_sec_fx1 = fx_second_nucleus(step_one_first, step_one_second);*/
	nuclei nuclei_k1 = pair_k(step_one_first, step_one_second);
	nuclei nuclei_k1_add = pair_add_dx(step_one_first, step_one_second, nuclei_k1, 2);
	nuclei nuclei_k2 = pair_k(nuclei_k1_add.first, nuclei_k1_add.second);
	nuclei nuclei_k2_add = pair_add_dx(step_one_first, step_one_second, nuclei_k2, 2);
	nuclei nuclei_k3 = pair_k(nuclei_k2_add.first, nuclei_k2_add.second);
	nuclei nuclei_k3_add = pair_add_dx(step_one_first, step_one_second, nuclei_k3, 1);
	nuclei nuclei_k4 = pair_k(nuclei_k3_add.first, nuclei_k3_add.second);


	k_one_to_four_add(nuclei_k1, nuclei_k2, nuclei_k3, nuclei_k4, step_one_first, step_one_second);
	return;
}

__device__ void update_step_two(nucleus& step_one_first, nucleus& step_one_second)
{

}

