


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
//#define DX 1
#define   DX    ((2.0*N1_const + N2_const) * T0_const / (one_steps))

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


















double h = 0.1;

double g1(double x1,double px1,double t1)
{
	return px1;
}

double f1(double x1,double px1,double t1)
{
	return -cos(t1);
}



void fill_t_arr(double* t_arr,size_t size)
{
	for(int i= 0 ;i!= size;i++){
		t_arr[i] = h*i;
	}
}

void first_step_every_step
(double* x_arr,double* px_arr,double* t_arr,size_t iter_times)
{
	for(size_t i = 0;i!=iter_times;i++){
		double t1 = t_arr[i];
		double y1[2] = {x_arr[i],px_arr[i]};
		double K1[2] = {
			g1(y1[0],y1[1],t1),
			f1(y1[0],y1[1],t1)
		};
		double y2[2] = {
			K1[0] * h /2.0 + x_arr[i],
			K1[1] * h /2.0 + px_arr[i]
		};
		double t2 = t_arr[i] + h/2.0;
		double K2[2] = {
			g1(y2[0],y2[1],t2),
			f1(y2[0],y2[1],t2)
		};
		double y3[2] = {
			K2[0] * h /2.0 + x_arr[i],
			K2[1] * h /2.0 + px_arr[i]
		};
		double t3 = t_arr[i] + h/2.0;
		double K3[2] = {
			g1(y3[0],y3[1],t2),
			f1(y3[0],y3[1],t2)
		};
		double y4[2] = {
			K3[0] * h /2.0 + x_arr[i],
			K3[1] * h /2.0 + px_arr[i]
		};
		double t4 = t_arr[i]  + h;
		double K4[2] = {
			g1(y4[0],y4[1],t2),
			f1(y4[0],y4[1],t2)
		};

		x_arr[i] = 


	}

}












 

 

 int main()
 {
	//  particle_pair init = {
	// 	 0.6499352856,
	// 	 -0.6063439233,
	// 	 -0.3058286632,
	// 	 -0.0390684108,
	// 	 -0.1062714257,
	// 	 -1.2950598609,
	// 	 -0.6499352856,
	// 	 0.6063439233,
	// 	 0.3058286632,
	// 	 1.1479844474,
	// 	 0.2725119169,
	// 	 0.5457737293
	//  };


	particle_pair init = {
		1.0,
		0,
		0,
		0,
		0,
		0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0

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
	//PrintK1K2K3K4(every_step_k,12*one_steps,"K1K2K3K4_add.dat");
	 PrintStruct(every_step, one_steps, "no_laser_every_step_dx_raw.dat");
	
 
 }
