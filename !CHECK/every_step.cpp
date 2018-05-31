


#include <ctime>
#include <cstdio>
#include <cmath>

struct double4
{
	double x, y, z, w;
};
struct double3
{
	double x, y, z;
};

double4 make_double4(double x, double y, double z, double w)
{
	double4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}

//龙哥库塔用到变量
#define   one_steps  10000

#define two_steps 40000


#define PI  3.1415926535897932384626433832795	//圆周率
//双电子

#define E0  1.59 // 基态能量

#define  E_sec_ion  1.065 // 第二电离能,Ip2

#define   Q_constant    1.225

#define   A_hardness    2.0 //硬度参数

#define   Omega1    0.0584

#define   Omega2    0.117

#define   N1_const    2.0

#define   N2_const    6.0

#define   T0_const    (2 * PI / Omega1)

//下面的10000为one_steps
#define DX 1
//#define   DX    ((2.0*N1_const + N2_const) * T0_const / (two_steps))

#define  TP_const 0.8

//初始矩阵相关
#define NX_const  351
#define NY_const  301


//粒子数
#define Pairs_Total 1UL

//激光场迭代次数

#define Iter_Count  21

//最开始  ee0 输出检查时
#define  EE0_Check  (sqrt(1e15 / 3.51e16))















struct derivative
{
	double px;
	double py;
	double pz;
	double fx;
	double fy;
	double fz;
};


struct particle {
	double x;
	double y;
	double z;

	double px;
	double py;
	double pz;
};

struct  particle_pair {
	particle first;
	particle second;
};



void PrintStruct(particle_pair* ToSaveNuclei, size_t n, const char* FileName) //输出的粒子的信息
{
	//format date and time. 
	struct tm *ptr;
	time_t lt;
	char str[9];


	FILE* file;
	file = fopen(FileName, "w");


	if (!file) perror("cannot open file");
	for (size_t i = 0; i < n; i++)
	{
		fprintf(file, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
			ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
			ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
			ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz
		);
		fprintf(file, "\n");
	}


	fclose(file);

	return;
}

void PrintK1K2K3K4(double4* array,size_t size,const char* FileName)
{
	FILE* file;
	file = fopen(FileName, "w");


	if (!file) perror("cannot open file");
	for (size_t i = 0;  (i < size); i++)
	{
		

		fprintf(file, "%-.12lf\t%-.12lf\t%-.12lf\t%-.12lf\t ",
			array[i].x,array[i].y,array[i].z,array[i].w
		);
		fprintf(file, "\n");
	
	}


	fclose(file);

	return;
	
}




   particle first_and_second_k_add_dx_raw(const derivative& k_one_to_four, const particle& raw_to_add)
 {
	 double now_dx = DX;

	 particle k_add;
	 k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	 k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	 k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	 k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	 k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	 k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	 return k_add;
 }

   particle first_and_second_k_add_dx_div(const derivative& k_one_to_four, const particle& raw_to_add)
 {
	 double now_dx = DX / 2.0;

	 particle k_add;
	 k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	 k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	 k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	 k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	 k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	 k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	 return k_add;
 }











 double nucleus_distance(const particle& first, const particle& second)
 {
	 return (pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2));
 }


 //第一个核，三个坐标的一阶导
   double3 gx_gy_gz_first_nucleus(const particle& first, const particle& second)
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
   double3 gx_gy_gz_second_nucleus(const particle& first, const particle& second)
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
   double3 fx_fy_fz_first_nucleus(const particle& first, const particle& second)
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
   double3 fx_fy_fz_second_nucleus(const particle& first, const particle& second)
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




 //第一个粒子 K1~K4 第二步循环
   derivative fisrt_k_one_to_four_second_step
 (const particle& first, const particle& second, const double& e1_laser, const double& e2_laser)
 {
	 const double3 first_fx = fx_fy_fz_first_nucleus(first, second);
	 const double3 first_gx = gx_gy_gz_first_nucleus(first, second);
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
   derivative second_k_one_to_four_second_step
 (const particle& first, const particle& second, const double& e1_laser, const double& e2_laser)
 {
	 const double3 second_fx = fx_fy_fz_second_nucleus(first, second);
	 const double3 second_gx = gx_gy_gz_second_nucleus(first, second);
	 derivative second_px_fx;
	 second_px_fx.px = second_gx.x;
	 second_px_fx.py = second_gx.y;
	 second_px_fx.pz = second_gx.z;
	 second_px_fx.fx = second_fx.x;
	 second_px_fx.fy = second_fx.y - e2_laser;
	 second_px_fx.fz = second_fx.z - e1_laser;
	 return second_px_fx;

 }










   void k_one_to_four_add(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4,
	 particle& raw_to_add)
 {
	 raw_to_add.x = raw_to_add.x + DX * (K1.px + 2.0*K2.px + 2.0*K3.px + K4.px) / 6.0;
	 raw_to_add.y = raw_to_add.y + DX * (K1.py + 2.0*K2.py + 2.0*K3.py + K4.py) / 6.0;
	 raw_to_add.z = raw_to_add.z + DX * (K1.pz + 2.0*K2.pz + 2.0*K3.pz + K4.pz) / 6.0;
	 raw_to_add.px = raw_to_add.px + DX * (K1.fx + 2.0*K2.fx + 2.0*K3.fx + K4.fx) / 6.0;
	 raw_to_add.py = raw_to_add.py + DX * (K1.fy + 2.0*K2.fy + 2.0*K3.fy + K4.fy) / 6.0;
	 raw_to_add.pz = raw_to_add.pz + DX * (K1.fz + 2.0*K2.fz + 2.0*K3.fz + K4.fz) / 6.0;



	 return;
 }

   //第一个粒子 K1~K4 第一步循环
 derivative fisrt_k_one_to_four_fisrt_step(const particle& first, const particle& second)
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


 derivative second_k_one_to_four_fisrt_step(const particle& first, const particle& second)
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





   void update_step_one(particle& step_one_first, particle& step_one_second,
	   particle_pair& every_step,double4* every_step_k)
   {
	   
	   
	   	every_step_k->x = step_one_first.x;
	   	(every_step_k + 1)->x = step_one_first.px;
	  	(every_step_k + 2)->x = step_one_second.x;
	   	(every_step_k + 3)->x = step_one_second.px;
		(every_step_k + 4)->x =step_one_first.y;
		(every_step_k + 5)->x = step_one_first.py;
		(every_step_k + 6)->x = step_one_second.y;
		(every_step_k + 7)->x = step_one_second.py;
		(every_step_k + 8)->x = step_one_first.z;
		(every_step_k + 9)->x = step_one_first.pz;
		(every_step_k + 10)->x =  step_one_second.z;
		(every_step_k + 11)->x =  step_one_second.pz;
	   
	   //计算K1
	   const derivative first_k1 = fisrt_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	   const derivative second_k1 = second_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	   const particle first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	   const particle second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	   //K2
	   const derivative first_k2 = fisrt_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	   const derivative second_k2 = second_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	   const particle first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	   const particle second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	   //K3
	   const derivative first_k3 = fisrt_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	   const derivative second_k3 = second_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	   const particle first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	   const particle second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	   //K4
	   const derivative first_k4 = fisrt_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);
	   const derivative second_k4 = second_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);

	   k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	   k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);


	   every_step.first.x = step_one_first.x;
	   every_step.first.y = step_one_first.y;
	   every_step.first.z = step_one_first.z;
	   every_step.first.px = step_one_first.px;
	   every_step.first.py = step_one_first.py;
	   every_step.first.pz = step_one_first.pz;


	   every_step.second.x = step_one_second.x;
	   every_step.second.y = step_one_second.y;
	   every_step.second.z = step_one_second.z;
	   every_step.second.px = step_one_second.px;
	   every_step.second.py = step_one_second.py;
	   every_step.second.pz = step_one_second.pz;


//g1	
	
	every_step_k->y = first_k1_add.x;
	every_step_k->z = first_k2_add.x;
	every_step_k->w = first_k3_add.x;
//f1
	//(every_step_k + 1)->x = first_k1_add.fx;
	(every_step_k + 1)->y = first_k1_add.px;
	(every_step_k + 1)->z = first_k2_add.px;
	(every_step_k + 1)->w = first_k3_add.px;
//g2

	//(every_step_k + 2)->x = second_k1_add.px;
	(every_step_k + 2)->y = second_k1_add.x;
	(every_step_k + 2)->z = second_k2_add.x;
	(every_step_k + 2)->w = second_k3_add.x;
//f2

	//(every_step_k + 3)->x = second_k1_add.px;
	(every_step_k + 3)->y = second_k1_add.px;
	(every_step_k + 3)->z = second_k2_add.px;
	(every_step_k + 3)->w = second_k3_add.px;
//g3
	//(every_step_k + 4)->x = first_k1_add.y;
	(every_step_k + 4)->y = first_k1_add.y;
	(every_step_k + 4)->z = first_k2_add.y;
	(every_step_k + 4)->w = first_k3_add.y;
	
//f3
	//(every_step_k + 5)->x = first_k1_add.py;
	(every_step_k + 5)->y = first_k1_add.py;
	(every_step_k + 5)->z = first_k2_add.py;
	(every_step_k + 5)->w = first_k3_add.py;

//g4
	//(every_step_k + 6)->x = second_k1_add.y;
	(every_step_k + 6)->y = second_k1_add.y;
	(every_step_k + 6)->z = second_k2_add.y;
	(every_step_k + 6)->w = second_k3_add.y;	
	
	
	
//f4
	//(every_step_k + 7)->x = second_k1_add.py;
	(every_step_k + 7)->y = second_k1_add.py;
	(every_step_k + 7)->z = second_k2_add.py;
	(every_step_k + 7)->w = second_k3_add.py;
	
	
//g5
	//(every_step_k + 8)->x = first_k1_add.z;
	(every_step_k + 8)->y = first_k1_add.z;
	(every_step_k + 8)->z = first_k2_add.z;
	(every_step_k + 8)->w = first_k3_add.z;
	
	
//f5
	//(every_step_k + 9)->x = first_k1_add.pz;
	(every_step_k + 9)->y = first_k1_add.pz;
	(every_step_k + 9)->z = first_k2_add.pz;
	(every_step_k + 9)->w = first_k3_add.pz;
//g6
	//(every_step_k + 10)->x = second_k1_add.z;
	(every_step_k + 10)->y = second_k1_add.z;
	(every_step_k + 10)->z = second_k2_add.z;
	(every_step_k + 10)->w = second_k3_add.z;
//f6
	//(every_step_k + 11)->x = second_k1_add.pz;
	(every_step_k + 11)->y = second_k1_add.pz;
	(every_step_k + 11)->z = second_k2_add.pz;
	(every_step_k + 11)->w = second_k3_add.pz;

	
	   return;
   }





 void update_step_two_every_step(particle& step_one_first, particle& step_one_second,
	 const double4 e1_laser_now, const double4 e2_laser_now, particle_pair& every_step)
 {
	 //计算K1
	 const derivative first_k1 = fisrt_k_one_to_four_second_step(step_one_first, step_one_second, e1_laser_now.x, e2_laser_now.x);
	 const derivative second_k1 = second_k_one_to_four_second_step(step_one_first, step_one_second, e1_laser_now.x, e2_laser_now.x);
	 const particle first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	 const particle second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	 //K2
	 const derivative first_k2 = fisrt_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	 const derivative second_k2 = second_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	 const particle first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	 const particle second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	 //K3
	 const derivative first_k3 = fisrt_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	 const derivative second_k3 = second_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	 const particle first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	 const particle second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	 //K4
	 const derivative first_k4 = fisrt_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);
	 const derivative second_k4 = second_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);

	 k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	 k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);

	 every_step.first.x = step_one_first.x;
	 every_step.first.y = step_one_first.y;
	 every_step.first.z = step_one_first.z;
	 every_step.first.px = step_one_first.px;
	 every_step.first.py = step_one_first.py;
	 every_step.first.pz = step_one_first.pz;


	 every_step.second.x = step_one_second.x;
	 every_step.second.y = step_one_second.y;
	 every_step.second.z = step_one_second.z;
	 every_step.second.px = step_one_second.px;
	 every_step.second.py = step_one_second.py;
	 every_step.second.pz = step_one_second.pz;
 }


 void first_step_every_step(particle_pair init,particle_pair* every_step,double4* every_step_k)
{
	for(size_t i =0 ;i<one_steps;i++)
	{
		update_step_one(init.first, init.second, every_step[i],&every_step_k[12 * i]);
	}
}




 

 double compute_qq_single(const size_t& now_step)
 {
	 double t1 = 0.5 * DX * (now_step + 1);
	 return  pow((sin(Omega1 / 2.0 / (2 * N1_const + N2_const)*t1)), 2);

 }


 double compute_e1_single(const size_t& now_step, const double& qq_now_single, const double& EE0)
 {
	 double tao = 0.0;
	 double t1 = 0.5 * DX * (now_step + 1);
	 return  (EE0 / (1.0 + TP_const)) * qq_now_single * sin(Omega1 * t1 + tao) -
		 (EE0*TP_const / (1.0 + TP_const)) * qq_now_single * sin(Omega2 * t1 + 2 * tao);
 }

 double compute_e2_single(const size_t& now_step, const double& qq_now_single, const double& EE0)
 {
	 double tao = 0.0;
	 double t1 = 0.5 * DX * (now_step + 1);
	 return  (EE0 / (1.0 + TP_const)) * qq_now_single * cos(Omega1 * t1 + tao) +
		 (EE0*TP_const / (1.0 + TP_const)) * qq_now_single * cos(Omega2 * t1 + 2 * tao);

 }




 void pairs_second_step_on_gpu_every_step
 (particle_pair second_arr, const size_t size, double* E1_array, double* E2_array,
	 particle_pair* every_step_arr)
 {
	

	 double4 e1_laser = make_double4(0.0, 0.0, 0.0, 0.0);
	 double4 e2_laser = make_double4(0.0, 0.0, 0.0, 0.0);
	 int idx_of_laser = -1; // 相当于nn
							//double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
							//double now_t = 0.0; //当前时间，相当于t(1)


	
	
		 for (int i = 0; i < two_steps; i++)
		 {

			 if (idx_of_laser == -1)
			 {
				 e1_laser = make_double4(0.0, E1_array[0], E1_array[0], E1_array[1]);
				 e2_laser = make_double4(0.0, E2_array[0], E2_array[0], E2_array[1]);
			 }
			 else
			 {
				 e1_laser = make_double4(E1_array[idx_of_laser], E1_array[idx_of_laser + 1], E1_array[idx_of_laser + 1], E1_array[idx_of_laser + 2]);
				 e2_laser = make_double4(E2_array[idx_of_laser], E2_array[idx_of_laser + 1], E2_array[idx_of_laser + 1], E2_array[idx_of_laser + 2]);
			 }
			 idx_of_laser += 2;

			 update_step_two_every_step(second_arr.first, second_arr.second,
				 e1_laser, e2_laser, every_step_arr[i]);


		 }


	
 }


 int main()
 {
	 particle_pair init = {
		 0.6499352856,
		 -0.6063439233,
		 -0.3058286632,
		 -0.0390684108,
		 -0.1062714257,
		 -1.2950598609,
		 -0.6499352856,
		 0.6063439233,
		 0.3058286632,
		 1.1479844474,
		 0.2725119169,
		 0.5457737293
	 };

	 double* qq_arr = new double[2 * two_steps];
	 double* e1_arr = new double[2 * two_steps];
	 double* e2_arr = new double[2 * two_steps];

	 for (size_t i = 0; i<2 * two_steps; i++)
	 {
		 qq_arr[i] = compute_qq_single(i); 
	 }

	 for (size_t i = 0; i<2 * two_steps; i++)
	 {
		 e1_arr[i] = compute_e1_single(i, qq_arr[i], EE0_Check);
		 e2_arr[i] = compute_e2_single(i, qq_arr[i], EE0_Check);
	 }


	 particle_pair* every_step = new particle_pair[one_steps];
	 double4* every_step_k = new double4[12 * one_steps];

	 //pairs_second_step_on_gpu_every_step(init, 1, e1_arr, e2_arr, every_step);
	 first_step_every_step(init, every_step,every_step_k);
	PrintK1K2K3K4(every_step_k,12*one_steps,"K1K2K3K4_add.dat");
	 //PrintStruct(every_step, one_steps, "no_laser_every_step.dat");
	
 
 }
