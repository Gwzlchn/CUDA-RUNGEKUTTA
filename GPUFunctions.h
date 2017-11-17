#ifndef GPUFUNCTIONS_H
#define GPUFUNCTIONS_H
#include<cuda_runtime.h>




__device__ const double PI=3.14159265358979323846; 
__device__ const double A=1;
__device__ const double E0=0.5;

//步长DX 终点TIME
__device__ const double DX=0.027;
__device__ const int TOSTOP=1000;

//void InitialRandom(double *ip, const int size);

void NormalRandom(double *ip, const int size);



void ComputeOnGPU1(double* Result,int nx,int ny,dim3 grid,dim3 block);









#endif //GPUFUNCTIONS_H





