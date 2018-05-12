
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


//下面的10000为one_steps
#include<iostream>


// int main(int argc, char const *argv[])
// {
//     double dx = DX;

//     printf("%.10f\n",T0_const);
//     printf("%.15f",dx);

//     return 0;
// }
