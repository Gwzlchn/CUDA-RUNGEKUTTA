#pragma comment(lib, "cudart.lib")
#pragma comment(lib, "curand.lib")

#include "../include/global_funcs.h"
#include "../include/sci_const.h"
#include "../include/device_compute_funcs.cuh"
#include "../include/common.hpp"
#include "../include/PrintStruct.h"

#include "device_launch_parameters.h"
#include <stdio.h>
#include <curand_kernel.h>
#include <cmath>
#include <vector_types.h>
#include <cuda_runtime.h>

//生成双精度01均匀分布随机数
//参数:	Array:双精度数组	Size:数组长度
//void UniformRandomArrayD(double* Array, const long Size)
//{
//	curandGenerator_t gen;											//生成随机数变量
//	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
//	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
//	curandGenerateUniformDouble(gen, Array, Size);					//生成0-1均匀分布随机数，存储到缓冲器中
//	curandDestroyGenerator(gen);                         			//释放指针
//	return;
//}
//
////生成双精度正态分布随机数
////参数:	Array:双精度数组	Size:数组长度	Mean:均值(0)	Stddev:方差(0.7)
//void NormalRandomArrayD(double* Array, const long Size, double Mean, double Stddev)
//{
//	curandGenerator_t gen;											//生成随机数变量
//	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);		//指定算法
//	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);					//随机数初始化
//	curandGenerateNormalDouble(gen, Array, Size, Mean, Stddev);		//生成正态分布随机数，存储到缓冲器中
//	curandDestroyGenerator(gen);                         			//释放指针
//	return;
//}


//准备矩阵，计算一次即可
//


__global__ void pre_step_init(nuclei* Array, const long& size,
				const double& min_r,const double& min_p)
{
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size)
	{
		distribution(Array[idx].first, Array[idx].second, idx, min_r, min_p);
	}
	return;
}

__global__ void first_step_on_gpu(nuclei* first_arr, const long size)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("%p\n", &first_arr);
	if(idx<size)
	{
		//printf("%d\n", idx);
		for (int i = 0; i < one_steps; i++)
			update_step_one(first_arr[idx].first, first_arr[idx].second);
	}
	
	
}


__global__ void pre_second_step_qq(double* QQ)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2 * two_steps)
	{

		double t1 = 0.5 * DX * idx;
		QQ[idx] = pow((sin(Omega1 / 2.0 / (2 * N1_const + N2_const))*t1), 2);
	    
		
	}
	
}


__global__ void pre_second_step_E_forcheck(const double* E1,const double* E2,double* E_check)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2 * two_steps)
	{

		double t1 = 0.5 * DX * idx;
		E_check[idx] = sqrt(pow(E1[idx], 2) + pow(E2[idx] ,2));
	    
		
	}
}








__global__ void pre_second_step_e1(const double* QQ,const double EE0,double* E1)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2*two_steps)
	{
		double t1 = 0.5 * DX * idx;
		curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;

		E1[idx] = (EE0 / (1.0 + TP_const)) * QQ[idx] * sin(Omega1 * t1 + tao) -
			(EE0*TP_const / (1.0 + TP_const)) * QQ[idx] * sin(Omega2 * t1 + 2 * tao);

	}

}


__global__ void pre_second_step_e2(const double* QQ, const double EE0, double* E2)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		double t1 = 0.5 * DX * idx;
		curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;

		E2[idx] = (EE0 / (1.0 + TP_const)) * QQ[idx] * sin(Omega1 * t1 + tao) +
			(EE0*TP_const / (1.0 + TP_const)) * QQ[idx] * sin(Omega2 * t1 + 2 * tao);

	}

}









__global__ void second_step_on_gpu(nuclei* second_arr , const long size, double* E1,double* E2)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	double e1_laser_t1=0.0, e1_laser_t2=0.0, e1_laser_t3=0.0, e1_laser_t4=0.0;
	double e2_laser_t1 = 0.0, e2_laser_t2 = 0.0, e2_laser_t3 = 0.0, e2_laser_t4 = 0.0;
	int idx_of_ds=-1; // 相当于nn
	double t1=0.0, t2=0.0, t3=0.0, t4=0.0;
	double now_t=0.0; //当前时间，相当于t(1)
	if (idx<size)
	{
		for (int i = 0; i < two_steps; i++)
		{
			//第一个激光场强度
			t1 = now_t;
			if (t1 == 0)
				e1_laser_t1 = 0.0;
			else
			{
				idx_of_ds = (2.0 * t1) / DX - 1;
				e1_laser_t1 = E1[idx_of_ds];
				e2_laser_t1 = E2[idx_of_ds];
			}
			//第二个激光场强度
			t2 = now_t + DX / 2.0;
			idx_of_ds = 2.0 * t2 / DX- 1;
			e1_laser_t2 = E1[idx_of_ds];
			e2_laser_t2 = E2[idx_of_ds];
			//第三个激光场强度
			t3 = now_t + DX / 2.0;
			idx_of_ds = 2 * t3 / DX- 1;
			e1_laser_t3 = E1[idx_of_ds];
			e2_laser_t3 = E2[idx_of_ds];
			//第四个激光场强度
			t4 = now_t + DX;
			idx_of_ds = 2.0 * t4 / DX - 1;
			e1_laser_t4 = E1[idx_of_ds];
			e2_laser_t4 = E2[idx_of_ds];
			double4 e1_laser = make_double4(e1_laser_t1, e2_laser_t2, e1_laser_t3, e1_laser_t4);
			double4 e2_laser = make_double4(e2_laser_t1, e2_laser_t2, e2_laser_t3, e2_laser_t4);
			update_step_two(second_arr[idx].first, second_arr[idx].second,
							e1_laser,e2_laser);
			now_t = now_t + DX;
			/*if(idx_of_ds == -1 )
				update_step_two(second_arr[idx].first, second_arr[idx].second,
									0.0,DS[0],DS[0],DS[1]);
			else
			{
				update_step_two(second_arr[idx].first, second_arr[idx].second,
					DS[idx_of_ds], DS[idx_of_ds + 1], DS[idx_of_ds + 1], DS[idx_of_ds + 2]);
			}
			idx_of_ds += 2;*/

		}
		
			
	}
}


/*
 * double ee1 = CalculationE1(second_arr[idx].first, second_arr[idx].second);
		double ee2 = CalculationE2(second_arr[idx].first, second_arr[idx].second);
		if (ee1>0 && ee2>0)
		{
			
			unsigned long long temp_idx = atomicAdd(ee1_ee2_count, 1);
			nuclei temp;
			temp.first = second_arr[idx].first;
			temp.second = second_arr[idx].second;
			second_arr_fliter[temp_idx-1] = temp;
		}
 */


void get_min_r_min_p(int nx, int ny, double& min_r, double& min_p)
{
	double *R_Arr = (double*)malloc(nx * sizeof(double));
	double *P_Arr = (double*)malloc(ny * sizeof(double));

	for (int i = 0; i < nx; i++)
		R_Arr[i] = 0.5 + 0.01 * i;
	for (int i = 0; i < ny; i++)
		P_Arr[i] = 0.0 + 0.01*i;


	double** mat = new double*[nx];
	for (int i = 0; i<nx; i++)
		mat[i] = new double[ny];

	double Vh, Vk, Ek;
	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			Vh = pow(Q_constant, 2) / (4.0*A_hardness*pow(R_Arr[i], 2)) *
				exp(A_hardness * (1.0 - pow((R_Arr[i] * P_Arr[j] / Q_constant), 4)));
			Vk = -2.0 / R_Arr[i];
			Ek = P_Arr[j] * P_Arr[j] / 2.0;
			mat[i][j] = Vh + Vk + Ek + 1.065;
		}
	}

	int min_x_index, min_y_index;
	double min = mat[0][0];
	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (min > mat[i][j])
			{
				min = mat[i][j];
				min_x_index = i;
				min_y_index = j;
			}
		}
	}

	min_r = R_Arr[min_x_index];
	min_p = P_Arr[min_y_index];

	return;
}





//用于双核粒子的随机数化  初始化
void NucleiPreRandom(nuclei* Array, const long size)
{
	//计算最小 r p;
	double min_r, min_p;
	get_min_r_min_p(NX_const, NY_const, min_r, min_p);
	
	
	int dimx = 256;
	dim3 block(dimx);
	dim3 grid((size + block.x - 1) / block.x, 1);
	pre_step_init <<< grid, block >>> (Array, size,min_r,min_p);
	CHECK(cudaDeviceSynchronize());
}


void NucleiFisrtStep(nuclei* first_array, const long size)
{
	int dimx = 32;
	dim3 block(dimx);
	dim3 grid((size + block.x - 1) / block.x, 1);
	first_step_on_gpu <<< grid, block >>> (first_array, size);
	
	
	CHECK( cudaGetLastError() );
	CHECK(cudaDeviceSynchronize());
	
}



void NucleiSecondStepPreQQ(double* QQ)
{
	//准备矢量势
	int pre_dimx = 512;
	dim3 pre_block(pre_dimx);
	dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
	pre_second_step_qq <<< pre_grid, pre_block >>> (QQ);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());

}

void NucleiSecondStepPreECheck(const double* QQ,const double EE0, double* E_check)
{
	//申请当前激光场 E1 E2 空间
	double *host_e1, *host_e2;
	double *gpu_e1, *gpu_e2;
	long bytes_of_e_laser = sizeof(double) * 2 * two_steps_in_host;
	CHECK(cudaMalloc((void **)(&gpu_e1), bytes_of_e_laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), bytes_of_e_laser));
	host_e1 = (double*)malloc(bytes_of_e_laser);
	host_e2 = (double*)malloc(bytes_of_e_laser);

	//准备矢量势
	int pre_dimx = 512;
	dim3 pre_block(pre_dimx);
	dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
	pre_second_step_e1 << < pre_grid, pre_block >> > (QQ, EE0, gpu_e1);
	pre_second_step_e2 << < pre_grid, pre_block >> > (QQ, EE0, gpu_e2);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaGetLastError());

	
	
	

	pre_second_step_E_forcheck <<< pre_grid, pre_block >>> (gpu_e1,gpu_e2,E_check);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}










void NucleiSecondStepOneLaser(nuclei* second_array , const long size,double* QQ,double EE0)
{
	//申请当前激光场 E1 E2 空间
	double *host_e1,*host_e2;
	double *gpu_e1, *gpu_e2;
	long bytes_of_e_laser = sizeof(double) * 2 * two_steps_in_host;
	CHECK(cudaMalloc((void **)(&gpu_e1), bytes_of_e_laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), bytes_of_e_laser));
	host_e1 = (double*)malloc(bytes_of_e_laser);
	host_e2 = (double*)malloc(bytes_of_e_laser);

	//准备矢量势
	int pre_dimx = 512;
	dim3 pre_block(pre_dimx);
	dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
	pre_second_step_e1 <<< pre_grid, pre_block >>> (QQ, EE0, gpu_e1);
	pre_second_step_e2 <<< pre_grid, pre_block >>> (QQ, EE0, gpu_e2);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaGetLastError());

	//计算第二步 一个激光场
 	int dimx = 32;
	dim3 block(dimx);
	dim3 grid((size + block.x - 1) / block.x, 1);
	second_step_on_gpu <<< grid, block >>> (second_array,size,gpu_e1,gpu_e2);
	
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}



void compute_on_gpu_one(const long pairs,const char* file_name)
{
	long long nBytes = pairs * sizeof(nuclei);
	printf("Use %lld Bytes %lfMB\n", nBytes, nBytes / double(1024 * 1024));
	nuclei *gpu_init,*gpu_first,*gpu_second,*gpu_second_fliter;
	nuclei *host_init,*host_first,*host_second,*host_second_fliter;
	host_init = (nuclei*)malloc(nBytes);
	host_first = (nuclei*)malloc(nBytes);
	host_second = (nuclei*)malloc(nBytes);
	host_second_fliter = (nuclei*)malloc(nBytes);


	//初始化！
	//申请init空间
	double start = seconds();
	CHECK(cudaMalloc((void **)(&gpu_init), nBytes));
	//计算
	NucleiPreRandom(gpu_init, pairs);

	//把值赋给第一步(也申请了第一步的空间)
	CHECK(cudaMalloc((void **)(&gpu_first), nBytes));
	CHECK(cudaMemcpy(gpu_first, gpu_init, nBytes, cudaMemcpyDeviceToDevice));
	//拷回并保存
	CHECK(cudaMemcpy(host_init, gpu_init, nBytes, cudaMemcpyDeviceToHost));
	
	PrintStruct(host_init, pairs, file_name, 0);
	//释放init空间
	CHECK(cudaFree(gpu_init));
	double elapse = seconds();
	printf("Inition compltete %lf\n", elapse - start);
	//初始化完成！


	//第一步计算
	//first空间在之前申请过了
	 start = seconds();
	//计算
	NucleiFisrtStep(gpu_first, pairs);
	

	//把值赋给第二步(也申请了第二步的空间)
	CHECK(cudaMalloc((void **)(&gpu_second), nBytes));
	CHECK(cudaMemcpy(gpu_second, gpu_first, nBytes, cudaMemcpyDeviceToDevice));
	//拷回并保存
	CHECK(cudaMemcpy(host_first, gpu_first, nBytes, cudaMemcpyDeviceToHost));
	PrintStruct(host_first, pairs, file_name, 1);
	//释放first空间
	CHECK(cudaFree(gpu_first));
	elapse = seconds();
	printf("FirstStep compltete %lf\n", elapse - start);
	//第一步完成！


	//第二步计算
	start = seconds();
	//准备导数
	double *gpu_qq;
	double *host_qq;
	long bytes_of_e_laser = sizeof(double) * 2 * two_steps_in_host;
	CHECK(cudaMalloc((void **)(&gpu_qq), bytes_of_e_laser));
	host_qq = (double*)malloc(bytes_of_e_laser);
	NucleiSecondStepPreQQ(gpu_qq);

	double *gpu_e_check, *host_e_check;
	CHECK(cudaMalloc((void **)(&gpu_e_check), bytes_of_e_laser));
	host_e_check = (double*)malloc(bytes_of_e_laser);
	double EE0 = sqrt(1e15 / 3.51e16);
	NucleiSecondStepPreECheck(gpu_qq, EE0, gpu_e_check);
	CHECK(cudaMemcpy(host_e_check, gpu_e_check, bytes_of_e_laser, cudaMemcpyDeviceToHost));
	PrintArray(host_e_check, 2 * two_steps_in_host, "e_check", 0);

	
	////电离率计数
	//unsigned long long*gpu_count,*host_count;
	//int bytes_of_u_long = sizeof(unsigned long long);
	//host_count = (unsigned long long*)malloc(bytes_of_u_long);
	//CHECK(cudaMalloc((void **)(&gpu_count), bytes_of_u_long));
	//CHECK(cudaMalloc((void **)(&gpu_second_fliter), nBytes));

	////检查E_check

	//



	////计算

	//NucleiSecondStep(gpu_second, gpu_second_fliter,pairs, gpu_aw, gpu_ds,gpu_count);

	////拷回并保存
	//CHECK(cudaMemcpy(host_second, gpu_second, nBytes, cudaMemcpyDeviceToHost));
	//CHECK(cudaMemcpy(host_second_fliter, gpu_second_fliter, nBytes, cudaMemcpyDeviceToHost));
	//CHECK(cudaMemcpy(host_aw, gpu_aw, bytes_of_aw_ds, cudaMemcpyDeviceToHost));
	//CHECK(cudaMemcpy(host_ds, gpu_ds, bytes_of_aw_ds, cudaMemcpyDeviceToHost));
	//CHECK(cudaMemcpy(host_count, gpu_count, bytes_of_u_long, cudaMemcpyDeviceToHost));
	//printf("%ld\n", *host_count);

	//PrintStruct(host_second, pairs,file_name , 2);
	//PrintArray(host_aw, 2 * two_steps_in_host, file_name, 0);
	//PrintArray(host_ds, 2 * two_steps_in_host, file_name, 1);
	//PrintStruct(host_second_fliter, *host_count, file_name, 3);
	////释放second空间
	//CHECK(cudaFree(gpu_second));
	//CHECK(cudaFree(gpu_aw));
	//CHECK(cudaFree(gpu_ds));
	//
	//elapse = seconds();
	//printf("SecondStep compltete %lf\n", elapse - start);
	// 第二步完成！

	//释放主机内存空间
	//free(host_aw);
	//free(host_ds);
	//free(host_first);
	//free(host_init);
	//free(host_init);




	return;
}