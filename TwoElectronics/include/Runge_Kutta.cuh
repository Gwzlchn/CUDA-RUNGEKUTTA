#ifndef __RUNGEKUTTA_CUH
#define __RUNGEKUTTA_CUH


#include "device_launch_parameters.h"
#include "Struct_Defines.h"



//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);


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


#endif // RUNGEKUTTA