#include"GPUFunctions.h"
#include "device_launch_parameters.h"
#include "HostFunctions.hpp"
#include "common.h"
#include<cuda_runtime.h>
#include<curand.h>
#include<math.h>


__device__ double fx(double x)
{

	return -(x/(pow(sqrt(pow(x,2.0)+pow(A,2.0)),3.0)));
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
	//第一列已经是随机数了
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
    if (ix < nx ){
        for (int iy = 0; iy <ny; iy++)
        {
            int idx = iy * nx + ix;
			//第二列为第一列各自的px初值，如果出现根号下小于零的情况，直接赋值0，计算部分判断简单些（nan判定很烦……）
            if((idx>=1*nx)&&(idx<2*nx)){
				if(Ekall(Result[idx-nx])>=0.0)
					Result[idx] = Px(double(Result[idx-nx]));
				else Result[idx] = 0.0;
			}
			//第三列为第一列各自的fx初值，出现小于零情况同理。
			if((idx>=2*nx)&&(idx<3*nx)){
				if(Result[idx-1*nx]>0.0)
					Result[idx] = fx(double(Result[idx-2*nx]));
				else Result[idx] = 0.0;
			}
			
			//第四五六列为前三列的复制，为了compute函数准备
			if((idx>=3*nx)&&(idx<4*nx)){
				if(Result[idx-2*nx]>0.0)
					Result[idx] = Result[idx-3*nx];
				else Result[idx] = 0.0;
			}
			if((idx>=4*nx)&&(idx<5*nx))
				Result[idx] = Result[idx-3*nx];
			
		}
	}

}


void NormalRandom(double *ip, const int size){
    

	curandGenerator_t gen;                                  //生成随机数变量
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);//步骤1：指定算法
    curandSetPseudoRandomGeneratorSeed(gen, 11ULL);         //步骤2：随机数初始化
    curandGenerateNormalDouble(gen, ip, size, 0, 0.7);        //步骤3：生成随机数，存储到缓冲器中（第1个数字为均值，第二个为方差）
    curandDestroyGenerator(gen);                         	//释放指针
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
	
	
	
	//保存数据仅仅为了测试用，写好compute部分以后肯定不用保存这个数据了……
	int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
	double *h_gpuRef;
	h_gpuRef = (double *)malloc(nBytes);
	CHECK(cudaMemcpy(h_gpuRef, d_Result, nBytes, cudaMemcpyDeviceToHost));
	//保存数据
	double iStart = seconds();
	StoreData(h_gpuRef,nx,ny,"init.dat");
	double iElaps = seconds() - iStart;
    printf("STORE THE InitialKernel DATA elapsed %lf sec\n",iElaps);
	
	
}






__device__ void updateXi(double& xi,double& pxi)
{
	//const double DX=0.027;
	double K1  = pxi;
	double K11 = fx(xi);
	
	double K2  = pxi + K11/2.0*DX;
	double K22 = fx(xi + K1/2.0*DX);
	
	double K3  = pxi + K22/2.0*DX;
	double K33 = fx(xi + K2/2.0*DX);
	
	double K4  = pxi+K33*DX;
	double K44 = fx(xi+K3*DX);
	
	xi  = xi  + DX * (K1  + 2*K2  + 2*K3  + K4)/6.0;
	pxi = pxi + DX * (K11 + K22*2 + K33*2 + K44)/6.0;
	return;
}










__global__ void ComputeKernel(double* Result,int nx,int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int idxOfXi  = 3 * nx + ix;
	unsigned int idxOfPxi = 4 * nx + ix;
	
	unsigned int Steps = TOSTOP/DX;
	
	
    if(ix<nx && Result[idxOfXi]!=0.0){
		for(int i=0;i<Steps;i++){
			updateXi(Result[idxOfXi],Result[idxOfPxi]);
			
		}
	}
}







 void ComputeOnGPU1(double* Result,int nx,int ny,double* h_gpuRef){
	
	
	//分配grid,block大小
	int dimx = 512;
    dim3 block(dimx);
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
    printf("STORE THE ComputeKernel DATA elapsed %lf sec\n",iElaps);
	return;
}









