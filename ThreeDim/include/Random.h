#
#ifndef RANDOM_H
#define RANDOM_H


#include "nucleus.hpp"


//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
void UniformRandomArrayD(double* Array, const long Size);

//生成双精度正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
//void NormalRandomArrayD(double* Array, const long Size);

//用于双核粒子的随机数化
//参数:	Array:粒子数组	Size:数组长度	Angle:偏移角(0)
void NucleiRandomD(nuclei* Array, const long Size);

#endif //RANDOM_H