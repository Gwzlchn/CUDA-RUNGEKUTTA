#include "../include/Init_First_Second.cuh"
#include "../include/Sci_Constant.h"
#include <curand.h>
#include <curand_kernel.h>





void get_min_r_min_p(int nx, int ny, double& min_r, double& min_p)
{
	double *R_Arr = (double*)malloc(nx * sizeof(double));
	double *P_Arr = (double*)malloc(ny * sizeof(double));

	for (int i = 0; i < nx; i++)
		R_Arr[i] = 0.5 + 0.01 * i;
	for (int i = 0; i < ny; i++)
		P_Arr[i] = 0.0 + 0.01*i;


	double** mat = new double*[nx];
	for (int i = 0; i<nx; i++)
		mat[i] = new double[ny];

	double Vh, Vk, Ek;
	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			Vh = pow(Q_constant, 2) / (4.0*A_hardness*pow(R_Arr[i], 2)) *
				exp(A_hardness * (1.0 - pow((R_Arr[i] * P_Arr[j] / Q_constant), 4)));
			Vk = -2.0 / R_Arr[i];
			Ek = P_Arr[j] * P_Arr[j] / 2.0;
			mat[i][j] = Vh + Vk + Ek + 1.065;
		}
	}

	int min_x_index, min_y_index;
	double min = mat[0][0];
	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (min > mat[i][j])
			{
				min = mat[i][j];
				min_x_index = i;
				min_y_index = j;
			}
		}
	}

	min_r = R_Arr[min_x_index];
	min_p = P_Arr[min_y_index];

	return;
}







__device__ void get_six_random(double2& two_random, double4& four_random, const int& seed)
{
	curandStatePhilox4_32_10_t s;
	curand_init(seed, 0, 0, &s);
	two_random = curand_uniform2_double(&s);
	four_random = make_double4(curand_uniform_double(&s), curand_uniform_double(&s),
	                           curand_uniform_double(&s), curand_uniform_double(&s));
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
		four_random.w = 0;


}


__device__ void distribution(particle& first, particle& second,
                              int seed,  double min_r, double min_p)
{

	double2 two_rand;
	double4 four_rand;
	get_six_random(two_rand, four_rand, seed);

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


__device__  double nucleus_distance(const particle& first, const particle& second)
{
	return (pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2));
}


//第一个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_first_nucleus(const particle& first, const particle& second)
{
	double Q_squre = pow(Q_constant, 2);
	//坐标平方和
	double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);

	//一阶导平方和
	double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);



	//第一个核 一阶导三个计算公式 对应 g1 g3 g5
	double3 gx_gy_gz;
	gx_gy_gz.x = first.px * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));

	gx_gy_gz.y = first.py * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));

	gx_gy_gz.z = first.pz * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));

	return gx_gy_gz;
}

//第二个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_second_nucleus(const particle& first, const particle& second)
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
__device__ double3 fx_fy_fz_first_nucleus(const particle& first, const particle& second)
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
__device__ double3 fx_fy_fz_second_nucleus(const particle& first, const particle& second)
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
