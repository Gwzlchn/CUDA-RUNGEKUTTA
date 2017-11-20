#ifndef GPUFUNCTIONS_H
#define GPUFUNCTIONS_H
#include<cuda_runtime.h>



//只有从main函数里调用的才需要在此处声明?

//科学常数应该定义在device.h里面，这里混淆了

__device__ const double PI=3.14159265358979323846; 
__device__ const double A=1;
__device__ const double E0=-0.5;
__device__ const double omega=0.057;
//步长DX 终点TIME
__device__ const double T0 = 2*PI/omega;
__device__ const int 	STEPS=40000;
__device__ const double DX=10*T0/STEPS;
__device__ const int STEPSFIRST=10000;
__device__ const int STEPSSECOND=40000;



void  InitialMatrix(double* d_Result,int nx,int ny);
void NormalRandom(double *ip, const int size);
void ComputeOnGPU1(double* Result,int nx,int ny,double* h_gpuRef);


#endif //GPUFUNCTIONS_H





