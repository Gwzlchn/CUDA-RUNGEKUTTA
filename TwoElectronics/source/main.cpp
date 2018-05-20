#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "../include/Erorr_Check.hpp"
#include "../include/Sci_Constant.h"


void compute_on_gpu_all(size_t pairs)
{
	
}






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


	compute_on_gpu_all(Pairs_Total);

	return 0;
}



