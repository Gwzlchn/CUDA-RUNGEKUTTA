#include "../include/global_funcs.h"
#include <curand.h>
#include "device_launch_parameters.h"

//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
void UniformRandomArray(double* Array, const long Size)
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
void NormalRandomArray(double* Array, const long Size, double Mean , double Stddev)
{
	curandGenerator_t gen;											//生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
	curandGenerateNormalDouble(gen, Array, Size, Mean, Stddev);		//生成正态分布随机数，存储到缓冲器中
	curandDestroyGenerator(gen);                         			//释放指针
	return;
}

//生成双精度双正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Nudis:核间距(2)	Stddev:方差(0.7)
void DoubleNormalRandomArray(double* Array, const long Size, double Nudis, double Stddev)
{
	UniformRandomArray(Array, Size);
	return;
}

//用于粒子的随机数化
//参数:	Array:粒子数组	Size:数组长度
void NucleusRandom(nucleus* Array, const long Size)
{
	double* arr = new double[Size];
	NormalRandomArray(arr, Size);
	for (int i = 0; i < Size; i++)
	{

	}
}

void NormalRandomNuclei(nuclei* raw_nuclei,double* random_arr ,const long n)
{
	for(long i=0;i<(n/2);i++)
	{
		raw_nuclei[i].init_first.x = random_arr[i];
		raw_nuclei[i].init_second.x = random_arr[i + (n / 2)];
	}
}

void NormalRandom(nuclei* raw_nuclei, const long n)
{

}
void InitialNuclei(nuclei* randomed_nuclei, const long raw_count, long & left)
{

}
void FirstStep(nuclei* inited_nuclei, const long n)
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;

}