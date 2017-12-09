#include"kernel_funcs.h"
#include"device_funcs.cuh"
#include<curand.h>
#include "common.hpp"
#include"host_funcs.hpp"

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


__global__ void ComputeKernel(double* Result,int nx,int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int idxOfXi  = 3 * nx + ix;
	unsigned int idxOfPxi = 4 * nx + ix;
	unsigned int idxOfXiTwo  = 5 * nx + ix;
	unsigned int idxOfPxiTwo = 6 * nx + ix;
	unsigned int idxOfTemp = 7 * nx + ix;


    if(ix<nx && Result[idxOfXi]!=0.0){
		for(int i=0;i<STEPSFIRST;i++){
			updateXi(Result[idxOfXi],Result[idxOfPxi]);
		}
		Result[idxOfXiTwo] = Result[idxOfXi];
		Result[idxOfPxiTwo] = Result[idxOfPxi];
		
		for(int i=0;i<STEPSSECOND;i++){
			updateXiAtStepTwo(Result[idxOfXiTwo],Result[idxOfPxiTwo],i*DX);
		}
		
		double TempE=0.5 * (pow(Result[idxOfPxiTwo],2.0)) - (1.0 / sqrt( pow(Result[idxOfXiTwo],2.0)+ pow(A,2.0)));
		if( TempE <= 0.0)
			Result[idxOfTemp]=-999;
	}
}



int CountZeros(double* h_Result,int nx)
{

	unsigned int idxOfXi  = nx ;
	unsigned int idxOfTemp = 7 * nx ;
	int count=0;
	for(int i=0;i<nx;i++){
		if(h_Result[idxOfXi+i] == 0.0f) count++;
		//if(h_Result[idxOfTemp+i] == -999) nonZeros++;
	}
	
	
	
	return count;
}

int CountTooBig(double* h_Result,int nx)
{

	unsigned int idxOfXi  = nx ;
	unsigned int idxOfTemp = 7 * nx ;
	int count=0;
	for(int i=0;i<nx;i++){
		//if(h_Result[idxOfXi+i] == 0.0) count++;
		if(h_Result[idxOfTemp+i] == -999) count++;
	}
	
	
	
	return count;
}









 void ComputeOnGPU1(double* Result,int nx,int ny,double* h_gpuRef){
	
	
	//分配grid,block大小
	int dimx = 512;
    dim3 block(dimx);
    dim3 grid((nx + block.x - 1) / block.x, 1);
	
	double iStart = seconds();
	ComputeKernel<<<grid,block>>>(Result,nx,ny);
	 CHECK(cudaDeviceSynchronize());
	//如果核函数错误，返回信息
    CHECK(cudaGetLastError());
	double iElaps = seconds() - iStart;
	printf("RungeOnGPU  elapsed %f sec\n",iElaps);
	
	// GPU数据拷贝回主机
	int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
	CHECK(cudaMemcpy(h_gpuRef, Result, nBytes, cudaMemcpyDeviceToHost));
	
	int zeros=0,nonzeros=0;
	zeros = CountZeros(h_gpuRef,nx);
	nonzeros = CountTooBig(h_gpuRef,nx);
	printf("The Number of Zeros is %d,\t The Number of NonZeros is %d \n",zeros,nonzeros);
	double per = (nx - zeros - nonzeros)/(nx - zeros);
	printf("Percentage is %lf  \n",per);
	
	//保存数据
	iStart = seconds();
	StoreData(h_gpuRef,nx,ny,"gpuStepTwo1202.dat");
	//StoreData(h_Random,1,ny,"h_Random.dat");
	iElaps = seconds() - iStart;
    printf("STORE THE ComputeKernel DATA elapsed %lf sec\n",iElaps);
	return;
}

