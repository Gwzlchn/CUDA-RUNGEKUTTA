#pragma comment(lib, "cudart.lib")
#pragma comment(lib, "curand.lib")
#include "../include/global_funcs.h"
#include "../include/sci_const.h"
#include "../include/device_compute_funcs.cuh"
#include "../include/common.hpp"
#include "../include/PrintStruct.h"

#include "device_launch_parameters.h"
#include <stdio.h>
#include <curand_kernel.h>
#include <cmath>
#include <vector_types.h>
#include <cuda_runtime.h>

//生成双精度双正态分布随机数
__global__ void DoubleNormalRandomArrayD(nuclei* Array, const long Size)
{
	double A1, A2, A3, A4;
	double Ekall = -1;
	double temp1 = 1;
	double temp2 = 1;

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	
	curandState s;
	int seed = -i;
	curand_init(seed, 0, 0, &s);
	
	while (Ekall < 0)
	{
		A2 = A4 = 2;

		while (A2 > temp1 && A4 > temp2)
		{
			A1 = curand_uniform_double(&s);
			A2 = curand_uniform_double(&s);
			A3 = curand_uniform_double(&s);
			A4 = curand_uniform_double(&s);
			
			A1 = (A1 - 0.5) * 20;
			A3 = (A3 - 0.5) * 20;

			temp1 = exp((-pow((A1 - mean), 2)) / (mean * stddev * stddev))
				+ exp((-pow((A1 + mean), 2)) / (mean * stddev * stddev));
			temp2 = exp((-pow((A3 - mean), 2)) / (mean * stddev * stddev))
				+ exp((-pow((A3 + mean), 2)) / (mean * stddev * stddev));
		}
		//printf("%lf\t%lf\n", A1,A3);
	
		Array[i].first.x = A1 * sin(rotation*PI);
		Array[i].first.y = 0;
		Array[i].first.z = A1 * cos(rotation*PI);

		Array[i].second.x = A3 * sin(rotation*PI);
		Array[i].second.y = 0;
		Array[i].second.z = A3 * cos(rotation*PI);

		Ekall = E_kall(Array[i].first, Array[i].second);
		
		//printf("%lf\n", Ekall);
	}
	px_py_pz_distribution(Array[i].first, Array[i].second,Ekall,i);
	return;
}

//用于双核粒子的随机数化
void NucleiRandomD(nuclei* Array, const long Size)
{
	int dimx = 512;
	dim3 block(dimx);
	dim3 grid((Size + block.x - 1) / block.x, 1);
	DoubleNormalRandomArrayD <<< grid, block >>> (Array, Size);
}

void compute_on_gpu_one(const long pairs)
{
	long long nBytes = pairs * sizeof(nuclei);
	printf("Use %lld Bytes %lfMB\n", nBytes, nBytes / double(1024 * 1024));

	nuclei* test;
	nuclei* host;
	double Start = seconds();
	cudaMalloc((void **)(&test), nBytes);
	host = (nuclei*)malloc(nBytes);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	NucleiRandomD(test, pairs);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	float costtime;
	cudaEventElapsedTime(&costtime, start, stop);

	cudaMemcpy(host, test, nBytes, cudaMemcpyDeviceToHost);
	PrintStruct(host, pairs, "testOne.dat", costtime,0);
}

//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
//void UniformRandomArrayD(double* Array, const long Size)
//{
//	curandGenerator_t gen;											//生成随机数变量
//	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
//	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
//	curandGenerateUniformDouble(gen, Array, Size);					//生成0-1均匀分布随机数，存储到缓冲器中
//	curandDestroyGenerator(gen);                         			//释放指针
//	return;
//}
//
////生成双精度正态分布随机数
////参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
//void NormalRandomArrayD(double* Array, const long Size, double Mean, double Stddev)
//{
//	curandGenerator_t gen;											//生成随机数变量
//	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
//	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
//	curandGenerateNormalDouble(gen, Array, Size, Mean, Stddev);		//生成正态分布随机数，存储到缓冲器中
//	curandDestroyGenerator(gen);                         			//释放指针
//	return;
//}

////生成双精度双正态分布随机数
////参数:	Array1:双精度数组1	Array2:双精度数组2	Array3:双精度数组3	Array2:双精度数组4	
////Size:数组长度	Nudis:半核间距(2)	Stddev:方差(0.7)
//__global__ void DoubleNormalRandomArrayD(double* Array1, double* Array2, double* Array3, double* Array4,
//	const long Size )
//{
//	int i = threadIdx.x;
//	double temp1 = 1;
//	double temp2 = 1;
//
//	Array1[i] = (Array1[i] - 0.5) * 20;
//	Array3[i] = (Array3[i] - 0.5) * 20;
//
//	temp1 = exp((-pow((Array1[i] - nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)))
//		+ exp((-pow((Array1[i] + nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)));
//	temp2 = exp((-pow((Array3[i] - nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)))
//		+ exp((-pow((Array3[i] + nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)));
//
//	if (Array2[i] > temp1 && Array4[i] > temp2)
//	{
//		Array1[i] = -99;
//		Array3[i] = -99;
//	}
//	return;
//}
//
////线性传参
//__global__ void LinearTransmissionD(nuclei* Array, double* DTempArr1, double* DTempArr3, const long Size, int& i, int& j)
//{
//	int p, q;
//	cudaMalloc((void **)(&p), 4);
//	cudaMalloc((void **)(&q), 4);
//	cudaMemcpy(&p, &i, 4, cudaMemcpyHostToDevice);
//	cudaMemcpy(&p, &i, 4, cudaMemcpyHostToDevice);
//	while (i < Size && (i + j) < 2 * Size)
//	{
//		if (DTempArr1[i + j] == -99)
//		{
//			j++;
//		}
//		else {
//			Array[i].first.x = DTempArr1[i + j] * sin(rotation*PI);
//			Array[i].first.y = 0;
//			Array[i].first.z = DTempArr1[i + j] * cos(rotation*PI);
//			Array[i].second.x = DTempArr3[i + j] * sin(rotation*PI);
//			Array[i].second.y = 0;
//			Array[i].second.z = DTempArr3[i + j] * cos(rotation*PI);
//			i++;
//		}
//	}
//	cudaMemcpy(&i, &p, 4, cudaMemcpyDeviceToHost);
//	cudaMemcpy(&j, &q, 4, cudaMemcpyDeviceToHost);
//	return;
//}
//
////用于双核粒子的随机数化
////参数:	Array:粒子数组	Size:数组长度 Angle:偏移角
//void NucleiRandomD(nuclei* Array, const long Size)
//{
//	int i(0);
//	int j(0);
//	size_t DoubleSize = 2 * Size * sizeof(double);
//	double *DTempArr1, *DTempArr2, *DTempArr3, *DTempArr4;
//	cudaMalloc((void**)&DTempArr1, DoubleSize);
//	cudaMalloc((void**)&DTempArr2, DoubleSize);
//	cudaMalloc((void**)&DTempArr3, DoubleSize);
//	cudaMalloc((void**)&DTempArr4, DoubleSize);
//
//	while (i < Size)
//	{
//		UniformRandomArrayD(DTempArr1, 2 * Size);
//		UniformRandomArrayD(DTempArr2, 2 * Size);
//		UniformRandomArrayD(DTempArr3, 2 * Size);
//		UniformRandomArrayD(DTempArr4, 2 * Size);
//
//		int threadsPerBlock = 256;
//		int threadsPerGrid = (2 * Size + threadsPerBlock - 1) / threadsPerBlock;
//		DoubleNormalRandomArrayD <<<threadsPerGrid, threadsPerBlock >>> (DTempArr1, DTempArr2, DTempArr3, DTempArr4, 2 * Size);
//		LinearTransmissionD <<<1,1>>>(Array, DTempArr1, DTempArr3, Size, i, j);
//	}
//}