#pragma once

#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH


#include "nucleus.hpp"
#include <crt/host_defines.h>



//计算总动能
__device__ double E_kall(const nucleus& first, const nucleus& second);
//Px Py Pz
__device__ void px_py_pz_distribution(nucleus& first, nucleus& second);
//两核之间距离的平方
//返回（x1-x2)^2 +（y1-y2)^2 +（z1-z2)^2
__device__  double nucleus_distance(const nucleus& first, const nucleus& second);

//第一个核，三个坐标的二阶导
__device__ nucleus fx_first_nucleus(const nucleus& first, const nucleus& second);
//第二个核，三个坐标的二阶导
__device__ nucleus fx_second_nucleus(const nucleus& first, const nucleus& second);



__device__ void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec);
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec);

#endif //DEVICE_COMPUTE_FUNCS_CUH



