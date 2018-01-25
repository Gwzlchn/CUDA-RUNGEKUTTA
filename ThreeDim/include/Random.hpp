#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include <curand.h>

#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include "nucleus.hpp"
#include "sci_const.cuh"

//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
void UniformRandomArrayD(double* Array, const long Size);

//生成双精度正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
void NormalRandomArrayD(double* Array, const long Size, double Mean = 0, double Stddev = 0.7);

//生成双精度双正态分布随机数
//参数:	Array1:双精度数组1	Array2:双精度数组2	Array3:双精度数组3	Array2:双精度数组4	
//Size:数组长度	Nudis:半核间距(2)	Stddev:方差(0.7)
__global__ void DoubleNormalRandomArrayD(double* Array1, double* Array2, double* Array3, double* Array4,
	const long Size, double Nudis = 2, double Stddev = 0.7);

//用于双核粒子的随机数化
//参数:	Array:粒子数组	Size:数组长度	Angle:偏移角(0)
void NucleiRandomD(nuclei* Array, const long Size, double Angle = 0);

#endif //RANDOM_H



//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
void UniformRandomArrayD(double* Array, const long Size)
{
	curandGenerator_t gen;											//生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
	curandGenerateUniformDouble(gen, Array, Size);					//生成0-1均匀分布随机数，存储到缓冲器中
	curandDestroyGenerator(gen);                         			//释放指针
	return;
}

//生成双精度正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
void NormalRandomArrayD(double* Array, const long Size, double Mean, double Stddev)
{
	curandGenerator_t gen;											//生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
	curandGenerateNormalDouble(gen, Array, Size, Mean, Stddev);		//生成正态分布随机数，存储到缓冲器中
	curandDestroyGenerator(gen);                         			//释放指针
	return;
}

//生成双精度双正态分布随机数
//参数:	Array1:双精度数组1	Array2:双精度数组2	Array3:双精度数组3	Array2:双精度数组4	
//Size:数组长度	Nudis:半核间距(2)	Stddev:方差(0.7)
__global__ void DoubleNormalRandomArrayD(double* Array1, double* Array2, double* Array3, double* Array4,
	const long Size, double Nudis, double Stddev)
{
	int i = threadIdx.x;
	double temp1 = 1;
	double temp2 = 1;
	double temp3 = 1;

	Array1[i] = (Array1[i] - 0.5) * 20;
	Array3[i] = (Array3[i] - 0.5) * 20;

	temp1 = exp((-pow((Array1[i] - Nudis), Nudis)) / (Nudis * pow(Stddev, Nudis)))
		+ exp((-pow((Array1[i] + Nudis), Nudis)) / (Nudis * pow(Stddev, Nudis)));
	temp2 = exp((-pow((Array3[i] - Nudis), Nudis)) / (Nudis * pow(Stddev, Nudis)))
		+ exp((-pow((Array3[i] + Nudis), Nudis)) / (Nudis * pow(Stddev, Nudis)));

	if (Array2[i] > temp1 && Array4[i] > temp2)
	{
		Array1[i] = -99;
		Array3[i] = -99;
	}
	return;
}

//用于双核粒子的随机数化
//参数:	Array:粒子数组	Size:数组长度 Angle:偏移角
void NucleiRandomD(nuclei* Array, const long Size, double Angle)
{
	int i(0);
	int j(0);
	size_t DoubleSize = 2 * Size * sizeof(double);
	double *DTempArr1, *DTempArr2, *DTempArr3, *DTempArr4;
	cudaMalloc((void**)&DTempArr1, DoubleSize);
	cudaMalloc((void**)&DTempArr2, DoubleSize);
	cudaMalloc((void**)&DTempArr3, DoubleSize);
	cudaMalloc((void**)&DTempArr4, DoubleSize);

	while (i<Size)
	{
		UniformRandomArrayD(DTempArr1, 2 * Size);
		UniformRandomArrayD(DTempArr2, 2 * Size);
		UniformRandomArrayD(DTempArr3, 2 * Size);
		UniformRandomArrayD(DTempArr4, 2 * Size);

		int threadsPerBlock = 256;
		int threadsPerGrid = (2 * Size + threadsPerBlock - 1) / threadsPerBlock;
		DoubleNormalRandomArrayD <<<threadsPerGrid, threadsPerBlock>>>(DTempArr1, DTempArr2, DTempArr3, DTempArr4, 2 * Size);
		while (i<Size && (i + j)<2 * Size)
		{
			if (DTempArr1[i + j] == -99)
			{
				j++;
			}
			else {
				Array[i].first.x = DTempArr1[i + j] * sin(Angle*PI);
				Array[i].first.y = 0;
				Array[i].first.z = DTempArr1[i + j] * cos(Angle*PI);
				Array[i].second.x = DTempArr3[i + j] * sin(Angle*PI);
				Array[i].second.y = 0;
				Array[i].second.z = DTempArr3[i + j] * cos(Angle*PI);
				i++;
			}
		}
	}
}