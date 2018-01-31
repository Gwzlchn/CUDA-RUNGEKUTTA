#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "./include/common.hpp"
#include "./include/global_funcs.h"

int main()
{
	//计时
	double Start, Elaps;

	printf("Starting...\n");

	//选择设备
	int dev = 0;
	cudaDeviceProp deviceProp;
	CHECK(cudaGetDeviceProperties(&deviceProp, dev));
	printf("Using Device %d: %s\n", dev, deviceProp.name);
	CHECK(cudaSetDevice(dev));

	long pairs = 10000;
	
	compute_on_gpu_one(pairs);

	return 0;
}



