#ifndef GLOBAL_FUNS_H
#define GLOBAL_FUNS_H
#include "nucleus.hpp"
#include <math.h>	//可以使用M_PI
#include <curand.h>

//生成双精度正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
void NormalRandomArray(double* Array, const long Size, double Mean = 0, double Stddev = 0.7);

//生成双精度双正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Nudis:核间距(2)	Stddev:方差(0.7)
void DoubleNormalRandomArray(double* Array, const long Size, double Nudis = 2, double Stddev = 0.7);

//用于粒子的随机数化
//参数:	Array:粒子数组	Size:数组长度
void NucleusRandom(nucleus* Array, const long Size);

void NormalRandom(nuclei* raw_nuclei, const long n);
void InitialNuclei(nuclei* randomed_nuclei, const long raw_count, long & left);

void FirstStep(nuclei* inited_nuclei, const long n);



#endif //GLOBAL_FUNS_H
