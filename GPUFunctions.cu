#include"GPUFunctions.h"
#include "device_launch_parameters.h"
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






void NormalRandom(double *ip, const int size){
    

	curandGenerator_t gen;                                  //生成随机数变量
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);//步骤1：指定算法
    curandSetPseudoRandomGeneratorSeed(gen, 11ULL);         //步骤2：随机数初始化
    curandGenerateNormalDouble(gen, ip, size, 0, 2);        //步骤3：生成随机数，存储到缓冲器中（第1个数字为均值，第二个为方差）
    curandDestroyGenerator(gen);                         //释放指针
	return;
	
	
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





__global__ void ComputeKernel(double* Result,int nx,int ny)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;

    if (ix < nx ){
		
            //int idx = iy * nx + ix;
             //Result[idx] = MatA[idx];
			 //Result[1*nx] = MatA[nx];
        
        for (int iy = 0; iy <ny; iy++)
        {
            int idx = iy * nx + ix;
            
			if((idx>=1*nx)&&(idx<2*nx))
				Result[idx] = Px(double(Result[idx-nx]));
			if((idx>=2*nx)&&(idx<3*nx))
				Result[idx] = fx(double(Result[idx-2*nx]));
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






 void ComputeOnGPU1(double* Result,int nx,int ny,dim3 grid,dim3 block){
	
	
	
	ComputeKernel<<<grid,block>>>(Result,nx,ny);
	return;
}









