
#ifndef SCI_CONST_H
#define SCI_CONST_H

#include <crt/host_defines.h>
#include <cmath>



//随机数有关
__device__ double stddev = 0.7;		//方差

__device__ double mean = 2.0;			//均值


//核粒子有关
__device__ double rotation = 0.5;	//转轴角度，对应之前的kk

__device__ double nuclear_spacing = 4.0;	//核间距,对应之前的R

__device__ double PI = 3.1415926535897932384626433832795;	//圆周率

__device__ double elec_nucl =1.25;	//电子与核之间参数，对应之前A



__device__ double elec_elec = 0.1;	//电子与电子之间参数。对应之前A1



//!!注意这个4.0应该是R,但windows不支持。服务器可以
__device__ double E_total = -1.0165 - 1.0 / 4.0;	//两个电子总能量 对应之前 E0



//龙哥库塔用到变量
__device__ int  one_steps = 10000;

__device__ int	two_steps = 40000;
int two_steps_in_host = 40000;


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

#define   DX    ((2.0*N1_const + N2_const) * T0_const / (40000))

#define  TP_const 0.8

//初始矩阵相关
const int NX_const = 351;
const int NY_const = 301;


#endif
