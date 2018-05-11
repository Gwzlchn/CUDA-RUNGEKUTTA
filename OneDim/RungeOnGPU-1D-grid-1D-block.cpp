#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "kernel_funcs.h"
#include "common.hpp"

//计时用变量
double iStart;
double iElaps;
	
	
	
int main()
{
   printf("Starting...\n");

    //选择设备
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    // 数据量，nx是粒子数
    int nx = 1000000;
    int ny = 9;
    int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
    printf("Matrix size: nx %d ny %d\n", nx, ny);

	
	
	//申请主机内存空间
	double *h_gpuRef;
	h_gpuRef = (double *)malloc(nBytes);
	//申请GPU内存空间
	double *d_Result;
    CHECK(cudaMalloc((void **)&d_Result, nBytes));
	
    //以随机数填充初始化，并启动核函数完成fx,px初始化
    iStart = seconds();
    InitialMatrix(d_Result,nx,ny);
    iElaps = seconds() - iStart;
    printf("initialize matrix elapsed %f sec\n", iElaps);

	
    
	//在一个GPU上启动核函数,并将值储存在外文件
    
	ComputeOnGPU1(d_Result,nx,ny,h_gpuRef);
    
    //printf("RungeOnGPU1  elapsed %f sec\n",iElaps);


	// 释放GPU内存，释放主机内机内存
    CHECK(cudaFree(d_Result));
    free(h_gpuRef);

    // reset device
    CHECK(cudaDeviceReset());


    return (0);
}
