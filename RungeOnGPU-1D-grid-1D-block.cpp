#include "common.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include "HostFunctions.hpp"
#include "GPUFunctions.h"


//#include<cutil_math.h>



//在GPU端用生成正态分布数组
	
double iStart;
double iElaps;
	
	
	
int main()
{
   printf("Starting...\n");

    // set up device
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    // set up data size of matrix
    int nx = 1 << 5;
    int ny = 1 << 3;
	
	
	
	

    int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
    printf("Matrix size: nx %d ny %d\n", nx, ny);
	
	double *h_gpuRef;
	h_gpuRef = (double *)malloc(nBytes);
	//memset(h_gpuRef, 0, nBytes);
	
	
	
	double *d_Result;
    CHECK(cudaMalloc((void **)&d_Result, nBytes));

    // initialize data at host side
    iStart = seconds();
    NormalRandom(d_Result, nx);
    iElaps = seconds() - iStart;
    printf("initialize matrix elapsed %f sec\n", iElaps);

   


    // invoke kernel at host side
	
	
	
	
	int dimx = 256;
    dim3 block(dimx, 1);
    dim3 grid((nx + block.x - 1) / block.x, 1);
    

    iStart = seconds();
	ComputeOnGPU1(d_Result,nx,ny,grid,block);
    //RungeOnGPU1D<<<grid, block>>>(d_Result, nx, ny);
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    printf("sumMatrixOnGPU1D  elapsed %f sec\n",iElaps);

    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(h_gpuRef, d_Result, nBytes, cudaMemcpyDeviceToHost));

    // check device results
    //checkResult(hostRef, gpuRef, nxy);

	
	
	
	
	

	//Store DATA	 
	iStart = seconds();
	StoreData(h_gpuRef,nx,ny,"gpu.dat");
	//StoreData(h_Random,1,ny,"h_Random.dat");
	iElaps = seconds() - iStart;
    printf("STORE THE DATA elapsed %lf sec\n",iElaps);
    
	
	
	// free device global memory
    CHECK(cudaFree(d_Result));
    

  
   
    
    free(h_gpuRef);

    // reset device
    CHECK(cudaDeviceReset());


    return (0);
}
