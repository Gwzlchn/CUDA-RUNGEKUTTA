#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "../include/Erorr_Check.hpp"



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


	compute_on_gpu_one(pairs, "TwoElec");

	return 0;
}



