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
#include <thr/xthrcommon.h>

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


__global__ void pre_step_init(nuclei* Array, const long size,
				const double min_r,const double min_p)
{
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size)
	{

		double2 two_rand;
		double4 four_rand;
		get_six_random(two_rand, four_rand, idx);

		double theta1 = two_rand.x * 2.0 * PI;
		double phi1 = two_rand.y * PI;


		Array[idx].first.x = min_r * sin(phi1) * cos(theta1);
		Array[idx].first.y = min_r * sin(phi1) * sin(theta1);
		Array[idx].first.z = min_r * cos(phi1);

		Array[idx].second.x = -Array[idx].first.x;
		Array[idx].second.y = -Array[idx].first.y;
		Array[idx].second.z = -Array[idx].first.z;


		double phi2 = four_rand.x * PI;
		double phi3 = four_rand.y * PI;
		double theta2 = four_rand.z * 2.0 * PI;
		double theta3 = four_rand.w * 2.0 * PI;

		Array[idx].first.px = min_p * cos(theta2)*sin(phi2);
		Array[idx].first.py = min_p * sin(theta2)*sin(phi2);
		Array[idx].first.pz = min_p * cos(phi2);

		Array[idx].second.px = min_p * cos(theta3)*sin(phi3);
		Array[idx].second.py = min_p * sin(theta3)*sin(phi3);
		Array[idx].second.pz = min_p * cos(phi3);
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

		double t1 = 0.5 * DX * (idx+1);
		QQ[idx] = pow((sin(Omega1 / 2.0 / (2 * N1_const + N2_const)*t1)), 2);
		//QQ[idx] = t1;
	}
	
}


__global__ void pre_second_step_E_forcheck(const double* E1,const double* E2,double* E_check)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2 * two_steps)
	{
		E_check[idx] = sqrt(pow(E1[idx], 2) + pow(E2[idx] ,2));
	}
}








__global__ void pre_second_step_e1(const double* QQ,const double EE0,double* E1)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < 2*two_steps)
	{
		double t1 = 0.5 * DX * idx;
		/*curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;*/
		double tao = 0.0;
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
	/*	curandStatePhilox4_32_10_t s;
		curand_init(idx, 0, 0, &s);
		double random = curand_uniform_double(&s);
		double tao = 2.0 * random * PI;*/
		double tao = 0.0;

		E2[idx] = (EE0 / (1.0 + TP_const)) * QQ[idx] * cos(Omega1 * t1 + tao) +
			(EE0*TP_const / (1.0 + TP_const)) * QQ[idx] * cos(Omega2 * t1 + 2 * tao);

	}

}









__global__ void second_step_on_gpu(nuclei* second_arr , const long size, double* E1,double* E2)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	double e1_laser_t1 = 0.0, e1_laser_t2 = 0.0, e1_laser_t3 = 0.0, e1_laser_t4 = 0.0;
	double e2_laser_t1 = 0.0, e2_laser_t2 = 0.0, e2_laser_t3 = 0.0, e2_laser_t4 = 0.0;
	int idx_of_ds=-1; // 相当于nn
	double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
	double now_t = 0.0; //当前时间，相当于t(1)
	if (idx<size)
	{
		for (int i = 0; i < two_steps; i++)
		{
			//第一个激光场强度
			t1 = now_t;
			if (t1 == 0)
			{
				e1_laser_t1 = 0.0;
				e2_laser_t1 = 0.0;
			}
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
			double4 e1_laser = make_double4(e1_laser_t1, e1_laser_t2, e1_laser_t3, e1_laser_t4);
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





__global__ void second_step_on_gpu_fliter
(const nuclei* second_arr, nuclei* second__arr_filter,
	const long size, unsigned long long* count_z, unsigned long long* count_zz)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < size)
	{
		double ee1 = CalculationE1(second_arr[idx].first, second_arr[idx].second);
		double ee2 = CalculationE2(second_arr[idx].first, second_arr[idx].second);

		if(ee1*ee2 < 0 )
		{
			atomicAdd(count_z, 1);
		}
		if( (ee1 > 0) && (ee2 > 0))
		{
			unsigned long long temp_idx = atomicAdd(count_zz, 1);
			/*nuclei temp;
			temp.first = second_arr[idx].first;
			temp.second = second_arr[idx].second;
			second__arr_filter[temp_idx - 1] = temp;*/
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

	
	PrintArray(host_e1, 2 * two_steps_in_host, "E1_check", 1);
	PrintArray(host_e2, 2 * two_steps_in_host, "E2_check", 1);
	

	pre_second_step_E_forcheck <<< pre_grid, pre_block >>> (gpu_e1,gpu_e2,E_check);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}




void NucleiSecondStepWholeLaserNoStream(nuclei* first_array, const long size, double* QQ)
{
	int n_streams = 1;


	unsigned long long *host_count_z_arr, *host_count_zz_arr;
	unsigned long long *gpu_count_z_arr, *gpu_count_zz_arr;
	const int size_ull = sizeof(unsigned long long);

	CHECK(cudaMalloc((void**)&gpu_count_z_arr, n_streams * size_ull));
	CHECK(cudaMalloc((void**)&gpu_count_zz_arr, n_streams * size_ull));

	//在CPU上分配页锁定内存  
	CHECK(cudaHostAlloc((void**)&host_count_z_arr, n_streams * size_ull, cudaHostAllocDefault));
	CHECK(cudaHostAlloc((void**)&host_count_zz_arr, n_streams * size_ull, cudaHostAllocDefault));
	for (int stream_index = 0; stream_index < n_streams; stream_index++)
	{
		double *gpu_e1, *gpu_e2;
		long bytes_of_e_laser = sizeof(double) * 2 * two_steps_in_host;
		CHECK(cudaMalloc((void **)(&gpu_e1), bytes_of_e_laser));
		CHECK(cudaMalloc((void **)(&gpu_e2), bytes_of_e_laser));

		double EE0 = 2.742*pow(10, 3)*sqrt(pow(10.0, (12.0 + double(stream_index)*0.2)));
		EE0 = EE0 / (5.1421*(pow(10.0, 11.0)));

		int pre_dimx = 512;
		dim3 pre_block(pre_dimx);
		dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
		pre_second_step_e1 << < pre_grid, pre_block, 0, 0 >> > (QQ, EE0, gpu_e1);
		pre_second_step_e2 << < pre_grid, pre_block, 0, 0 >> > (QQ, EE0, gpu_e2);


		//计算第二步 一个激光场
		int dimx = 32;
		dim3 block(dimx);
		dim3 grid((size + block.x - 1) / block.x, 1);

		nuclei* gpu_second_arr_once, *gpu_second_filter_once;

		long long nBytes = size * sizeof(nuclei);
		CHECK(cudaMalloc((void **)(&gpu_second_arr_once), nBytes));
		CHECK(cudaMalloc((void **)(&gpu_second_filter_once), nBytes));

		CHECK(cudaMemcpy(gpu_second_arr_once, first_array, nBytes, cudaMemcpyDeviceToDevice));

		second_step_on_gpu << < grid, block, 0, 0 >> > (gpu_second_arr_once, size, gpu_e1, gpu_e2);

		second_step_on_gpu_fliter << < grid, block, 0, 0 >> > (gpu_second_arr_once,
			gpu_second_filter_once, size, gpu_count_z_arr + stream_index, gpu_count_zz_arr + stream_index);

		CHECK(cudaMemcpy(host_count_z_arr + stream_index * size_ull,
			gpu_count_z_arr + size_ull * (stream_index),
			size_ull, cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(host_count_zz_arr + size_ull * (stream_index),
			gpu_count_zz_arr + size_ull * (stream_index),
			size_ull, cudaMemcpyDeviceToHost));




		//printf("第一列z,第二列zz");
		printf("%.10f\t", EE0);
		printf("z: %lld \t", host_count_z_arr[stream_index]);
		printf("zz: %lld \n", host_count_zz_arr[stream_index]);
		CHECK(cudaGetLastError());
	}
	CHECK(cudaGetLastError());
	//CHECK(cudaDeviceSynchronize());
}








void NucleiSecondStepWholeLaser(nuclei* first_array, const long size, double* QQ)
{
	int n_streams = 1;
	cudaStream_t *streams = (cudaStream_t *)malloc(n_streams * sizeof(cudaStream_t));

	for (int i = 0; i < n_streams; i++)
	{
		CHECK(cudaStreamCreate(&(streams[i])));
	}

	unsigned long long *host_count_z_arr, *host_count_zz_arr;
	unsigned long long *gpu_count_z_arr, *gpu_count_zz_arr;
	const int size_ull = sizeof(unsigned long long);
	
	CHECK(cudaMalloc((void**)&gpu_count_z_arr, n_streams * size_ull));
	CHECK(cudaMalloc((void**)&gpu_count_zz_arr, n_streams * size_ull));

	//在CPU上分配页锁定内存  
	CHECK(cudaHostAlloc((void**)&host_count_z_arr, n_streams * size_ull, cudaHostAllocDefault));
	CHECK(cudaHostAlloc((void**)&host_count_zz_arr, n_streams * size_ull, cudaHostAllocDefault));
	for(int stream_index = 0 ; stream_index < n_streams ; stream_index++)
	{
		double *gpu_e1, *gpu_e2;
		long bytes_of_e_laser = sizeof(double) * 2 * two_steps_in_host;
		CHECK(cudaMalloc((void **)(&gpu_e1), bytes_of_e_laser));
		CHECK(cudaMalloc((void **)(&gpu_e2), bytes_of_e_laser));
		
		double EE0 = 2.742*pow(10, 3)*sqrt(pow(10.0, (12.0 + double(stream_index)*0.2)));
		EE0 = EE0 / (5.1421*(pow(10.0, 11.0)));
				
		int pre_dimx = 512;
		dim3 pre_block(pre_dimx);
		dim3 pre_grid((2 * two_steps_in_host + pre_block.x - 1) / pre_block.x, 1);
		pre_second_step_e1 <<< pre_grid, pre_block ,0, streams[stream_index]>>> (QQ, EE0, gpu_e1);
		pre_second_step_e2 <<< pre_grid, pre_block , 0, streams[stream_index] >>> (QQ, EE0, gpu_e2);
		

		//计算第二步 一个激光场
		int dimx = 32;
		dim3 block(dimx);
		dim3 grid((size + block.x - 1) / block.x, 1);

		nuclei* gpu_second_arr_once,*gpu_second_filter_once;
		
		long long nBytes = size * sizeof(nuclei);
		CHECK(cudaMalloc((void **)(&gpu_second_arr_once), nBytes));
		CHECK(cudaMalloc((void **)(&gpu_second_filter_once), nBytes));

		CHECK(cudaMemcpy(gpu_second_arr_once, first_array, nBytes, cudaMemcpyDeviceToDevice));
		
		second_step_on_gpu <<< grid, block, 0, streams[stream_index] >>> (gpu_second_arr_once, size, gpu_e1, gpu_e2);

		second_step_on_gpu_fliter <<< grid, block, 0, streams[stream_index] >>> (gpu_second_arr_once, 
			gpu_second_filter_once, size, gpu_count_z_arr + stream_index , gpu_count_zz_arr + stream_index);
		
		cudaMemcpyAsync(host_count_z_arr + stream_index * size_ull,
			gpu_count_z_arr + size_ull * (stream_index),
			size_ull, cudaMemcpyDeviceToHost, streams[stream_index]);
		cudaMemcpyAsync(host_count_zz_arr + size_ull * (stream_index),
			gpu_count_zz_arr + size_ull * (stream_index),
			size_ull, cudaMemcpyDeviceToHost, streams[stream_index]);
	
		//printf("第一列z,第二列zz");
		printf("%.10f\t", EE0);
		printf("z: %lld \t", host_count_z_arr[stream_index]);
		printf("zz: %lld \n", host_count_zz_arr[stream_index]);
		CHECK(cudaGetLastError());
	}
	CHECK(cudaGetLastError());
	//CHECK(cudaDeviceSynchronize());
	for (int i = 0; i < n_streams; i++)
	{
		CHECK(cudaStreamDestroy(streams[i]));
		
	}
	



	free(streams);
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
	CHECK(cudaMemcpy(host_qq, gpu_qq, bytes_of_e_laser, cudaMemcpyDeviceToHost));
	//验证QQ 通过！
	//PrintArray(host_qq, 2 * two_steps_in_host, "QQ_Check", 0);


	////测试E通过！
	//double *gpu_e_check, *host_e_check;
	//CHECK(cudaMalloc((void **)(&gpu_e_check), bytes_of_e_laser));
	//host_e_check = (double*)malloc(bytes_of_e_laser);
	//double EE0 = sqrt(1e15 / 3.51e16);
	//NucleiSecondStepPreECheck(gpu_qq, EE0, gpu_e_check);
	//CHECK(cudaMemcpy(host_e_check, gpu_e_check, bytes_of_e_laser, cudaMemcpyDeviceToHost));
	//PrintArray(host_e_check, 2 * two_steps_in_host, "e_check", 0);
 
	NucleiSecondStepWholeLaser(gpu_second, pairs, gpu_qq);




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