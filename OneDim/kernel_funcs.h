#ifndef KERNEL_FUNCS_H
#define KERNEL_FUNCS_H

void  InitialMatrix(double* d_Result,int nx,int ny);
void NormalRandom(double *ip, const int size);
void ComputeOnGPU1(double* Result,int nx,int ny,double* h_gpuRef);
void CountZeros(double* h_Result,int nx,int& Zeros, int& nonZeros);

#endif