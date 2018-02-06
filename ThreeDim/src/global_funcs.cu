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

////生成双精度双正态分布随机数
////参数:	Array1:双精度数组1	Array2:双精度数组2	Array3:双精度数组3	Array2:双精度数组4	
////Size:数组长度	Nudis:半核间距(2)	Stddev:方差(0.7)
//__global__ void DoubleNormalRandomArrayD(double* Array1, double* Array2, double* Array3, double* Array4,
//	const long Size )
//{
//	int i = threadIdx.x;
//	double temp1 = 1;
//	double temp2 = 1;
//
//	Array1[i] = (Array1[i] - 0.5) * 20;
//	Array3[i] = (Array3[i] - 0.5) * 20;
//
//	temp1 = exp((-pow((Array1[i] - nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)))
//		+ exp((-pow((Array1[i] + nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)));
//	temp2 = exp((-pow((Array3[i] - nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)))
//		+ exp((-pow((Array3[i] + nuclear_spacing/2.0), nuclear_spacing/2.0)) / (nuclear_spacing/2.0 * pow(stddev, nuclear_spacing/2.0)));
//
//	if (Array2[i] > temp1 && Array4[i] > temp2)
//	{
//		Array1[i] = -99;
//		Array3[i] = -99;
//	}
//	return;
//}
//
////线性传参
//__global__ void LinearTransmissionD(nuclei* Array, double* DTempArr1, double* DTempArr3, const long Size, int& i, int& j)
//{
//	int p, q;
//	cudaMalloc((void **)(&p), 4);
//	cudaMalloc((void **)(&q), 4);
//	cudaMemcpy(&p, &i, 4, cudaMemcpyHostToDevice);
//	cudaMemcpy(&p, &i, 4, cudaMemcpyHostToDevice);
//	while (i < Size && (i + j) < 2 * Size)
//	{
//		if (DTempArr1[i + j] == -99)
//		{
//			j++;
//		}
//		else {
//			Array[i].first.x = DTempArr1[i + j] * sin(rotation*PI);
//			Array[i].first.y = 0;
//			Array[i].first.z = DTempArr1[i + j] * cos(rotation*PI);
//			Array[i].second.x = DTempArr3[i + j] * sin(rotation*PI);
//			Array[i].second.y = 0;
//			Array[i].second.z = DTempArr3[i + j] * cos(rotation*PI);
//			i++;
//		}
//	}
//	cudaMemcpy(&i, &p, 4, cudaMemcpyDeviceToHost);
//	cudaMemcpy(&j, &q, 4, cudaMemcpyDeviceToHost);
//	return;
//}
//
////用于双核粒子的随机数化
////参数:	Array:粒子数组	Size:数组长度 Angle:偏移角
//void NucleiRandomD(nuclei* Array, const long Size)
//{
//	int i(0);
//	int j(0);
//	size_t DoubleSize = 2 * Size * sizeof(double);
//	double *DTempArr1, *DTempArr2, *DTempArr3, *DTempArr4;
//	cudaMalloc((void**)&DTempArr1, DoubleSize);
//	cudaMalloc((void**)&DTempArr2, DoubleSize);
//	cudaMalloc((void**)&DTempArr3, DoubleSize);
//	cudaMalloc((void**)&DTempArr4, DoubleSize);
//
//	while (i < Size)
//	{
//		UniformRandomArrayD(DTempArr1, 2 * Size);
//		UniformRandomArrayD(DTempArr2, 2 * Size);
//		UniformRandomArrayD(DTempArr3, 2 * Size);
//		UniformRandomArrayD(DTempArr4, 2 * Size);
//
//		int threadsPerBlock = 256;
//		int threadsPerGrid = (2 * Size + threadsPerBlock - 1) / threadsPerBlock;
//		DoubleNormalRandomArrayD <<<threadsPerGrid, threadsPerBlock >>> (DTempArr1, DTempArr2, DTempArr3, DTempArr4, 2 * Size);
//		LinearTransmissionD <<<1,1>>>(Array, DTempArr1, DTempArr3, Size, i, j);
//	}
//}

//生成双精度双正态分布随机数double3

__global__ void DoubleNormalRandomArrayD(nuclei* Array, const long Size)
{
	
	double A1, A2, A3, A4;
	double Ekall = -1;
	double temp1 = 1;
	double temp2 = 1;

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < Size)
	{
		curandState s;
		int seed = -i;
		curand_init(seed, 0, 0, &s);

		while (Ekall < 0)
		{
			A2 = A4 = 2;

			while (A2 > temp1 && A4 > temp2)
			{
				A1 = curand_uniform_double(&s);
				A2 = curand_uniform_double(&s);
				A3 = curand_uniform_double(&s);
				A4 = curand_uniform_double(&s);

				A1 = (A1 - 0.5) * 20;
				A3 = (A3 - 0.5) * 20;

				temp1 = exp((-pow((A1 - mean), 2)) / (mean * stddev * stddev))
					+ exp((-pow((A1 + mean), 2)) / (mean * stddev * stddev));
				temp2 = exp((-pow((A3 - mean), 2)) / (mean * stddev * stddev))
					+ exp((-pow((A3 + mean), 2)) / (mean * stddev * stddev));
			}
			//printf("%lf\t%lf\n", A1,A3);

			Array[i].first.x = A1 * sin(rotation*PI);
			Array[i].first.y = 0;
			Array[i].first.z = A1 * cos(rotation*PI);

			Array[i].second.x = A3 * sin(rotation*PI);
			Array[i].second.y = 0;
			Array[i].second.z = A3 * cos(rotation*PI);

			Ekall = E_kall(Array[i].first, Array[i].second);

			//printf("%lf\n", Ekall);
		}
		px_py_pz_distribution(Array[i].first, Array[i].second, Ekall, i);
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


__global__ void pre_second_step_aw(double* AW)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2 * two_steps)
	{

		double field_strength = sqrt((2.8e15) / (3.51e16)); // 场强，对应之前ee0

		double t0 = 2 * PI / omega;
		double t1 = 0.5 * DX * idx;
		AW[idx] = (field_strength / omega) * (pow(sin((PI * t1) / (10 * t0)), 2)) * cos(omega * t1);
	}
	
}

__global__ void pre_second_step_ds(double* AW,double* DS)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2*two_steps)
	{
		if (idx == 0)
			DS[idx] = (AW[1] - AW[0]) / (0.5*DX);
		if (idx == (2 * two_steps - 1))
			DS[idx] = (AW[idx] - AW[idx - 1]) / (0.5*DX);
		else
		{
			DS[idx] = (AW[idx + 1] - AW[idx - 1]) / 2.0 /(0.5* DX);
		}
	}

}


__global__ void second_step_on_gpu(nuclei* second_arr, nuclei* second_arr_fliter , const long size,double* DS,unsigned long long* ee1_ee2_count)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	double e_laser_t1=0.0, e_laser_t2=0.0, e_laser_t3=0.0, e_laser_t4=0.0;
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
				e_laser_t1 = 0.0;
			else
			{
				idx_of_ds = (2.0 * t1) / DX - 1;
				e_laser_t1 = DS[idx_of_ds];
			}
			//第二个激光场强度
			t2 = now_t + DX / 2.0;
			idx_of_ds = 2.0 * t2 / DX- 1;
			e_laser_t2 = DS[idx_of_ds];
			//第三个激光场强度
			t3 = now_t + DX / 2.0;
			idx_of_ds = 2 * t3 / DX- 1;
			e_laser_t3 = DS[idx_of_ds];
			//第四个激光场强度
			t4 = now_t + DX;
			idx_of_ds = 2.0 * t4 / DX - 1;
			e_laser_t4 = DS[idx_of_ds];
			update_step_two(second_arr[idx].first, second_arr[idx].second,
							e_laser_t1,e_laser_t2,e_laser_t3,e_laser_t4);
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
		double ee1 = CalculationE1(second_arr[idx].first, second_arr[idx].second);
		double ee2 = CalculationE2(second_arr[idx].first, second_arr[idx].second);
		if (ee1>0 && ee2>0)
		{
			nuclei temp;
			temp.first = second_arr[idx].first;
			temp.second = second_arr[idx].second;
			second_arr_fliter[*ee1_ee2_count] = temp;
			atomicAdd(ee1_ee2_count, 1);
		}
			
	}
}

//用于双核粒子的随机数化
void NucleiRandomD(nuclei* Array, const long Size)
{
	int dimx = 512;
	dim3 block(dimx);
	dim3 grid((Size + block.x - 1) / block.x, 1);
	DoubleNormalRandomArrayD <<< grid, block >>> (Array, Size);
	CHECK(cudaDeviceSynchronize());
}


void NucleiFisrtStep(nuclei* first_array, const long size)
{
	int dimx = 32;
	dim3 block(dimx);
	dim3 grid((size + block.x - 1) / block.x, 1);
	first_step_on_gpu <<< grid, block >>> (first_array, size);
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "1st Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return;
	}
	CHECK(cudaDeviceSynchronize());
	
}




void NucleiSecondStep(nuclei* second_array, nuclei* second_array_fliter, const long size, double* aw, double* ds, unsigned long long* count)
{
	//准备矢量势
	int pre_dimx = 512;
	dim3 pre_block(pre_dimx);
	dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
	pre_second_step_aw <<< pre_grid,pre_block>>> (aw);
	CHECK(cudaDeviceSynchronize());
	pre_second_step_ds <<< pre_grid, pre_block >>> (aw, ds);
	CHECK(cudaDeviceSynchronize());
	
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "2nd Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return;
	}

	//计算第二步
	int dimx = 32;
	dim3 block(dimx);
	dim3 grid((size + block.x - 1) / block.x, 1);
	second_step_on_gpu <<< grid, block >>> (second_array, size, ds,count);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "2nd Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return;
	}
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
	NucleiRandomD(gpu_init, pairs);

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
	double *gpu_aw, *gpu_ds;
	double *host_aw, *host_ds;
	long bytes_of_aw_ds = sizeof(double) * 2 * two_steps_in_host;
	CHECK(cudaMalloc((void **)(&gpu_aw), bytes_of_aw_ds));
	CHECK(cudaMalloc((void **)(&gpu_ds), bytes_of_aw_ds));
	host_aw = (double*)malloc(bytes_of_aw_ds);
	host_ds = (double*)malloc(bytes_of_aw_ds);

	
	//电离率计数
	unsigned long long*gpu_count,*host_count;
	int bytes_of_u_long = sizeof(unsigned long long);
	host_count = (unsigned long long*)malloc(bytes_of_u_long);
	CHECK(cudaMalloc((void **)(&gpu_count), bytes_of_u_long));
	CHECK(cudaMalloc((void **)(&gpu_second_fliter), nBytes));

	//计算

	NucleiSecondStep(gpu_second, gpu_second_fliter,pairs, gpu_aw, gpu_ds,gpu_count);

	//拷回并保存
	CHECK(cudaMemcpy(host_second, gpu_second, nBytes, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(host_second_fliter, gpu_second_fliter, nBytes, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(host_aw, gpu_aw, bytes_of_aw_ds, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(host_ds, gpu_ds, bytes_of_aw_ds, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(host_count, gpu_count, bytes_of_u_long, cudaMemcpyDeviceToHost));
	printf("%ld\n", *host_count);

	PrintStruct(host_second, pairs,file_name , 2);
	PrintArray(host_aw, 2 * two_steps_in_host, file_name, 0);
	PrintArray(host_ds, 2 * two_steps_in_host, file_name, 1);
	PrintStruct(host_second_fliter, *host_count, file_name, 3);
	//释放second空间
	CHECK(cudaFree(gpu_second));
	CHECK(cudaFree(gpu_aw));
	CHECK(cudaFree(gpu_ds));
	
	elapse = seconds();
	printf("SecondStep compltete %lf\n", elapse - start);
	// 第二步完成！

	//释放主机内存空间
	//free(host_aw);
	//free(host_ds);
	//free(host_first);
	//free(host_init);
	//free(host_init);




	return;
}