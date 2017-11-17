#include"GPUFunctions.h"
#include "device_launch_parameters.h"
#include "HostFunctions.hpp"
#include "common.h"
#include<cuda_runtime.h>
#include<curand.h>
#include<math.h>


__device__ double fx(double x)
{

	return -1.0/(pow(sqrt(pow(x,2.0)+pow(A,2.0)),3.0));
}


__device__ double Ekall(double x)
{

	return E0+1.0/(sqrt(pow(x,2.0)+pow(A,2.0)));
}

__device__ double Px(double x)
{
	return sqrt(2*Ekall(x));
}

//数据初始化应该单独用一个kernel函数，计算fx px的初值
//待完成。mark一下
//1118wzl已完成

__global__ void InitialKernel(double* Result,int nx,int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
    if (ix < nx ){
        for (int iy = 0; iy <ny; iy++)
        {
            int idx = iy * nx + ix;
            if((idx>=1*nx)&&(idx<2*nx)){
				if(Ekall(Result[idx-nx])>=0.0)
					Result[idx] = Px(double(Result[idx-nx]));
				else Result[idx] = 0.0;
			}
				
			if((idx>=2*nx)&&(idx<3*nx)){
				if(Result[idx-1*nx]>0.0)
					Result[idx] = fx(double(Result[idx-2*nx]));
				else Result[idx] = 0.0;
			}
				
			if((idx>=3*nx)&&(idx<4*nx)){
				if(Result[idx-2*nx]>0.0)
					Result[idx] = Result[idx-3*nx];
				else Result[idx] = 0.0;
			}
			if((idx>=4*nx)&&(idx<5*nx))
				Result[idx] = Result[idx-3*nx];
			if((idx>=5*nx)&&(idx<6*nx))
				Result[idx] = Result[idx-3*nx];
		}
	}

}


void NormalRandom(double *ip, const int size){
    

	curandGenerator_t gen;                                  //生成随机数变量
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);//步骤1：指定算法
    curandSetPseudoRandomGeneratorSeed(gen, 11ULL);         //步骤2：随机数初始化
    curandGenerateNormalDouble(gen, ip, size, 0, 1);        //步骤3：生成随机数，存储到缓冲器中（第1个数字为均值，第二个为方差）
    curandDestroyGenerator(gen);                         //释放指针
	return;
	
	
}


void  InitialMatrix(double* d_Result,int nx,int ny){
	NormalRandom(d_Result,nx);
	//分配grid,block大小
	int dimx = 256;
    dim3 block(dimx, 1);
    dim3 grid((nx + block.x - 1) / block.x, 1);
	InitialKernel<<<grid,block>>>(d_Result,nx,ny);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaGetLastError());
	int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
	double *h_gpuRef;
	h_gpuRef = (double *)malloc(nBytes);
	CHECK(cudaMemcpy(h_gpuRef, d_Result, nBytes, cudaMemcpyDeviceToHost));
	//保存数据
	double iStart = seconds();
	StoreData(h_gpuRef,nx,ny,"init.dat");
	//StoreData(h_Random,1,ny,"h_Random.dat");
	double iElaps = seconds() - iStart;
    printf("STORE THE DATA elapsed %lf sec\n",iElaps);
	
	
}






__device__ double updateXi(double xi,double pxi)
{
	double K1=pxi;
	double K2=xi+K1/2.0*DX;
	double K3=xi+K2/2.0*DX;
	double K4=xi+K3*DX;
	
	return xi+DX*(K1+2*K2+2*K3+K4)/6.0;
}


__device__ double updatePxi(double pxi,double fxi)
{
	double K1=fxi;
	double K2=pxi+K1/2.0*DX;
	double K3=pxi+K2/2.0*DX;
	double K4=pxi+K3*DX;
	
	return pxi+DX*(K1+2*K2+2*K3+K4)/6.0;
}





__global__ void ComputeKernel(double* Result,int nx,int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;

    if (ix < nx ){
		
		return;
		
		}
}







 void ComputeOnGPU1(double* Result,int nx,int ny,double* h_gpuRef){
	
	
		//分配grid,block大小
	int dimx = 256;
    dim3 block(dimx, 1);
    dim3 grid((nx + block.x - 1) / block.x, 1);
	ComputeKernel<<<grid,block>>>(Result,nx,ny);
	 CHECK(cudaDeviceSynchronize());
	    //如果核函数错误，返回信息
    CHECK(cudaGetLastError());
	 // GPU数据拷贝回主机
	int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
	CHECK(cudaMemcpy(h_gpuRef, Result, nBytes, cudaMemcpyDeviceToHost));
	//保存数据
	double iStart = seconds();
	StoreData(h_gpuRef,nx,ny,"gpu.dat");
	//StoreData(h_Random,1,ny,"h_Random.dat");
	double iElaps = seconds() - iStart;
    printf("STORE THE DATA elapsed %lf sec\n",iElaps);
	return;
}









