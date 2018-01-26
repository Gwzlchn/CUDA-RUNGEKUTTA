#ifndef GLOBAL_FUNCS_H
#define GLOBAL_FUNCS_H

//#include <math.h>
//#include <thrust/iterator/counting_iterator.h> 
//#include <thrust/functional.h> 
//#include <thrust/transform_reduce.h> 
//#include <curand_kernel.h> 
//#include <curand.h>
//#include <cstdio>
//#include <cstring>
//
//#include "device_launch_parameters.h"
//#include "cuda_runtime.h"
//#include "../include/sci_const.h"
//
//#include "device_compute_funcs.cuh"
#include "nucleus.hpp"
////生成双精度01均匀分布随机数
////参数:	Array:双精度数组	Size:数组长度
//void UniformRandomArrayD(double* Array, const long Size);
//
////生成双精度正态分布随机数
////参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
////void NormalRandomArrayD(double* Array, const long Size);
//
////用于双核粒子的随机数化
////参数:	Array:粒子数组	Size:数组长度	Angle:偏移角(0)
//void NucleiRandomD(nuclei* Array, const long Size);

//void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName);

//用于双核粒子的随机数化

void compute_on_gpu_one(const long pairs);

void NucleiRandomD(nuclei* Array, const long Size);

void NucleiFisrtStep(nuclei* first_array, const long size);


void NucleiSecondStep(nuclei* second_array, const long size);

#endif //GLOBAL_FUNCS_H