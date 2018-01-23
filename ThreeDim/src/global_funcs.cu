#include "../include/global_funcs.h"
#include <curand.h>

//生成双精度正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Variance:方差(0.7)
void NormalRandomArray(double* Array, const long Size, int Mean=0, double Variance=0.7)
{
	curandGenerator_t gen;											//生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//步骤1：指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//步骤2：随机数初始化
	curandGenerateNormalDouble(gen, Array, Size, Mean, Variance);	//步骤3：生成随机数，存储到缓冲器中
	curandDestroyGenerator(gen);                         			//释放指针
	return;
}

//生成双精度双正态分布随机数
//参数:	Array:双精度数组	Size:数组长度	Nudis:核间距(2)	Variance:方差(0.7)
void DoubleNormalRandomArray(double* Array, const long Size, int Nudis = 2, double Variance = 0.7)
{
	curandGenerator_t gen;											//生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//步骤1：指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//步骤2：随机数初始化
	curandGenerateNormalDouble(gen, Array, Size, Nudis, Variance);	//步骤3：生成随机数，存储到缓冲器中
	curandDestroyGenerator(gen);                         			//释放指针
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

}