#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include<iostream>
#include<fstream>
#include "GPUFunctions.h"
#include "common.h"
#include "HostFunctions.hpp"

//计时用变量
double iStart;
double iElaps;
	
	
using namespace std;
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
    int nx = 100000;
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
    
	ifstream ReadInRandom;
	ReadInRandom.open("initOneCol.dat");
	int i=0;
	 while(!ReadInRandom.eof()) //读取数据到数组  
        {  
  
            ReadInRandom>>h_gpuRef[i*ny];  
            //file>>tempChar[i];  
            i++;  
        }  
	
	
	ReadInRandom.close();
	
	StoreData(h_gpuRef,nx,ny,"test.dat");
	return 0;
	
}