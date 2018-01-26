
#include "./include/common.hpp"
#include "./include/global_funcs.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

////生成双精度01均匀分布随机数
////参数:	Array:双精度数组	Size:数组长度
//void UniformRandomArrayD(double* Array, const long Size);
//
////生成双精度正态分布随机数
////参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
//void NormalRandomArrayD(double* Array, const long Size, double Mean = 0, double Stddev = 0.7);
//
////用于双核粒子的随机数化
////参数:	Array:粒子数组	Size:数组长度	Angle:偏移角(0)
//extern "C" void NucleiRandomD(nuclei* Array, const long Size, double Angle = 0);

int main()
{
	double rotation = 0.0;//转轴角度，对应之前的kk
	//计时
	double Start, Elaps;

	printf("Starting...\n");

	//选择设备
	int dev = 0;
	cudaDeviceProp deviceProp;
	CHECK(cudaGetDeviceProperties(&deviceProp, dev));
	printf("Using Device %d: %s\n", dev, deviceProp.name);
	CHECK(cudaSetDevice(dev));

	long pairs = 100;
	
	
	compute_on_gpu_one(pairs);


	return 0;
}



