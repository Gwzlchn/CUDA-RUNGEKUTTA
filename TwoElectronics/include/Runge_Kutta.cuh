#ifndef __RUNGEKUTTA_CUH
#define __RUNGEKUTTA_CUH


#include "device_launch_parameters.h"
#include "Struct_Defines.h"



//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const particle& first, const particle& second);



//龙哥库塔方法
__device__ void update_step_one(particle& step_one_first, particle& step_one_second);
__device__ void update_step_two(particle& step_one_first, particle& step_one_second,
	const double4 e1_laser_now, const double4 e2_laser_now);



//第一个粒子 K1~K4 第一步循环
__device__ derivative first_k_one_to_four_first_step
(const particle& first, const particle& second);
//第二个粒子 K1~K4 第一步循环
__device__ derivative second_k_one_to_four_first_step
(const particle& first, const particle& second);

//计算完 K3 更新下一步参数用，K3 不除以2
__device__ particle first_and_second_k_add_dx_raw
(const derivative& k_one_to_four, const particle& raw_to_add);
//计算完 K1或 K2 更新下一步参数用，K1 或 K2 除以2
__device__ particle first_and_second_k_add_dx_div
(const derivative& k_one_to_four, const particle& raw_to_add);
//计算完 K1~K4 完整的相加 
__device__ void k_one_to_four_add
(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4, particle& raw_to_add);


//第二步runge-kutta
__device__ derivative first_k_one_to_four_second_step
(const particle& first, const particle& second, const double& e1_laser, const double& e2_laser);


__device__ derivative second_k_one_to_four_second_step
(const particle& first, const particle& second, const double& e1_laser, const double& e2_laser);


//每一步都保存
__device__ void update_step_two_every_step(particle& step_one_first, particle& step_one_second,
	const double4 e1_laser_now, const double4 e2_laser_now, particle_pair& every_step);

#endif // RUNGEKUTTA