#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI    3.14159265358979323846 
#define rd (rand()/(RAND_MAX+1.0))  //此处是关于随机生成正态分布的定义

__device__ const double A=1;
__device__ const double E0=0.5;

//步长DX 终点TIME
__device__ const double DX=0.0027;
__device__ const int TOSTOP=10000;

//#include<cutil_math.h>

/*
 * This example demonstrates a simple vector sum on the GPU and on the host.
 * sumArraysOnGPU splits the work of the vector sum across CUDA threads on the
 * GPU. A 1D thread block and 1D grid are used. sumArraysOnHost sequentially
 * iterates through vector elements on the host.
 */

 void StoreData(double *Matrix, const int NX,const int NY,const char name[])
{
   
	
	FILE* fp;
	fp = fopen(name, "w");
    if (!fp)
    {
        perror("cannot open file");
	}
	
    for (int i = 0; i < NX; i++)
    {
		
		for(int j=0;j<NY;j++){
			fprintf(fp,"%-.10lf\t\t",*(Matrix+j*NX+i));
			
		}
		fprintf(fp,"\n");
       
    }

	return;
}

//区间[min,max]上的均匀分布
double rand_m(double min, double max)
{
    return min+(max-min)*rand()/(RAND_MAX+1.0);
}

//求均值为miu，方差为sigma的正态分布函数在x处的函数值
double normal(double x, double miu,double sigma)
{
    return 1.0/sqrt(2*PI)/sigma*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}

//按照矩形区域在函数值曲线上下位置分布情况得到正态分布函数x值
double rand_normal_distribution(double miu,double sigma, double min ,double max)
{
    double x,y,dScope;
    do{
        x=rand_m(min,max);
        y=normal(x,miu,sigma);
        dScope=rand_m(0.0,normal(miu,miu,sigma));
    }while(dScope>y);
    return x;
}

void initialData(double *ip, const int size)
{
    int i;

    for(i = 0; i < size; i++)
    {
        ip[i] = rand_normal_distribution(0.0,1.0,-10.0,10.0);
    }
	
    return;
}

__device__ double f(double x,double y)
{
	return cos(x);
}

/*__device__ double rk4(double dx, double x, double y)
{
	double	k1 = dx * f(x, y),
		k2 = dx * f(x + dx / 2, y + k1 / 2),
		k3 = dx * f(x + dx / 2, y + k2 / 2),
		k4 = dx * f(x + dx, y + k3);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}*/



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

__device__ double updateXi(double xi,double dx)
{
	double tempPx=Px(xi);
	double K1=tempPx,
		K2=xi+(tempPx/2.0)*dx,
		K3=xi+(((tempPx/2.0)*dx)+xi)/2.0*dx,
		K4=xi+(((((tempPx/2.0)*dx)+xi)/2.0*dx)+xi)*dx;
	
	return xi+dx*(K1+2*K2+2*K3+K4)/6.0;
}

__device__ double updatePxi(double xi,double dx)
{
	double tempPx=Px(xi),tempFx=fx(xi);
	double K1=tempFx,
		K2=tempPx+(tempFx/2.0)*dx,
		K3=tempPx+(((tempFx/2.0)*dx)+tempPx)/2.0*dx,
		K4=tempPx+(((((tempFx/2.0)*dx)+xi)/2.0*dx)+tempPx)*dx;
	
	return tempPx+dx*(K1+2*K2+2*K3+K4)/6.0;
}


// grid 1D block 1D
__global__ void RungeOnGPU1D(double *MatA,double *Result, int nx, int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;

    if (ix < nx ){
		
            //int idx = iy * nx + ix;
             //Result[idx] = MatA[idx];
			 //Result[1*nx] = MatA[nx];
        
        for (int iy = 0; iy <ny; iy++)
        {
            int idx = iy * nx + ix;
            if((idx<=1*nx))
				Result[idx] = MatA[idx];
			if((idx>=1*nx)&&(idx<2*nx))
				Result[idx] = Px(double(MatA[idx-nx]));
			if((idx>=2*nx)&&(idx<3*nx))
				Result[idx] = fx(double(MatA[idx-2*nx]));
			if((idx>=3*nx)&&(idx<4*nx)){
				Result[idx] = Result[idx-3*nx];
				int i,n=1+(TOSTOP)/DX;
				for(i=1;i<n;i++)
					Result[idx]=updateXi(Result[idx],DX);
			}
			if((idx>=4*nx)&&(idx<5*nx)){
				Result[idx] = Result[idx-3*nx];
				int i,n=1+(TOSTOP)/DX;
				for(i=1;i<n;i++)
					Result[idx]=updatePxi(Result[idx],DX);
			}
			
				
				/*const double dx=0.00001;
				int i,n=1+(2*PI)/dx;
				double temp;
				for (i = 1; i < n; i++){
					temp=rk4(dx, Result[idx-2*nx] + dx * (i - 1), Result[idx]);
					Result[idx] = temp;*/
					
		}
	}

}
		
	
	

	



int main()
{
    FILE* TIME_USED;
	TIME_USED = fopen("TimeData", "w");
    if (!TIME_USED)
    {
        perror("cannot open file");
	}
	
	
	printf("Starting...\n");
	fprintf(TIME_USED,"Starting...\n");

    // set up device
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using Device %d: %s\n", dev, deviceProp.name);
	fprintf(TIME_USED,"Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    // set up data size of matrix
    int nx = 1 << 10;
    int ny = 1 << 3;
	
	
	
	

    int nxy = nx * ny;
    int nBytes = nxy * sizeof(double);
    printf("Matrix size: nx %d ny %d\n", nx, ny);
	fprintf(TIME_USED,"Matrix size: nx %d ny %d\n", nx, ny);

    // malloc host memory
    double *h_Random, *gpuRef;
   
    
    gpuRef = (double *)malloc(nBytes);
    h_Random = (double *)malloc(nBytes);

    // initialize data at host side
    double iStart = seconds();
    initialData(h_Random, nxy);
    double iElaps = seconds() - iStart;
    printf("initialize matrix elapsed %f sec\n", iElaps);
	fprintf(TIME_USED,"initialize matrix elapsed %f sec\n", iElaps);

   
    memset(gpuRef, 0, nBytes);

   
	
	

    // malloc device global memory
    double *d_Random, *d_Result;
    CHECK(cudaMalloc((void **)&d_Random, nBytes));
    
    CHECK(cudaMalloc((void **)&d_Result, nBytes));

    // transfer data from host to device
    CHECK(cudaMemcpy(d_Random, h_Random, nBytes, cudaMemcpyHostToDevice));
    

    // invoke kernel at host side
    int dimx = 256;
    dim3 block(dimx, 1);
    dim3 grid((nx + block.x - 1) / block.x, 1);

    iStart = seconds();
    RungeOnGPU1D<<<grid, block>>>(d_Random, d_Result, nx, ny);
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    printf("sumMatrixOnGPU1D <<<(%d,%d), (%d,%d)>>> elapsed %f sec\n", 
			grid.x,grid.y,block.x, block.y, iElaps);
	fprintf(TIME_USED,"sumMatrixOnGPU1D <<<(%d,%d), (%d,%d)>>> elapsed %f sec\n", 
			grid.x,grid.y,block.x, block.y, iElaps);
    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_Result, nBytes, cudaMemcpyDeviceToHost));

    // check device results
    //checkResult(hostRef, gpuRef, nxy);

	
	
	
	
	

	//Store DATA	 
	iStart = seconds();
	StoreData(gpuRef,nx,ny,"gpu.dat");
	StoreData(h_Random,nx,1,"h_Random.dat");
	iElaps = seconds() - iStart;
    printf("STORE THE DATA elapsed %lf sec\n",iElaps);
	fprintf(TIME_USED,"STORE THE DATA elapsed %lf sec\n",iElaps);
    
	
	
	// free device global memory
    CHECK(cudaFree(d_Random));
    

    // free host memory
    free(h_Random);
   
    
    free(gpuRef);

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
