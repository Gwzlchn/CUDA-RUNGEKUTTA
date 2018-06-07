#include "../include/Call_GPU.cuh"
#include "../include/Erorr_Check.hpp"
#include "../include/Init_First_Second.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Compute_On_GPU.cuh"
#include "../include/Laser.cuh"
#include "../include/PrintStruct.h"
#include <cstdio>
#include <cuda_runtime_api.h>

//typedef unsigned long long size_t;

dim3 get_pre_block(int dimx)
{
	return dim3(dimx);
}

dim3 get_compute_block(int dimx)
{
	return dim3(dimx);
}

dim3 get_grid(size_t size, const dim3& block)
{
	return dim3((size + block.x - 1) / block.x, 1);
}

void SaveArraysWhichOnGPU(double* gpu_array, size_t size, const char* file_name)
{
	double* host_array = (double*)malloc(sizeof(double) * size);
	CHECK(cudaMemcpy(host_array, gpu_array, sizeof(double) * size , cudaMemcpyDeviceToHost));
	PrintArray(host_array, size, file_name);

}

void SaveLaserArraysWhichOnGPU(double* e1_array, double* e2_array, double* e_check_array, size_t size,
	const char* file_name)
{
	size_t bytes_of_arr = sizeof(double) * size ;

	double* e1_host_array = (double*)malloc(bytes_of_arr);
	double* e2_host_array = (double*)malloc(bytes_of_arr);
	double* e_check_host_array = (double*)malloc(bytes_of_arr);

	CHECK(cudaMemcpy(e1_host_array, e1_array, bytes_of_arr, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(e2_host_array, e2_array, bytes_of_arr, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(e_check_host_array, e_check_array, bytes_of_arr, cudaMemcpyDeviceToHost));

	PrintLaserArrays(e1_host_array,e2_host_array,e_check_host_array, size, file_name);




}

void SavePairsWhichOnGPU(particle_pair* gpu_array, size_t size, const char* file_name)
{
	particle_pair* host_pairs = (particle_pair*)malloc(size * sizeof(particle_pair));
	CHECK(cudaMemcpy(host_pairs, gpu_array, Bytes_Of_Pairs, cudaMemcpyDeviceToHost));
	PrintStruct(host_pairs, size, file_name);
}


void Pairs_Init_Call_GPU(particle_pair * pair_array_gpu, const size_t size)
{
	//计算最小 r p;
	double min_r, min_p;
	get_min_r_min_p(NX_const, NY_const, min_r, min_p);


	
	dim3 block = get_pre_block();;
	dim3 grid = get_grid(size,block);
	pairs_init <<< grid, block >>> (pair_array_gpu, size, min_r, min_p);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());

}

void Pairs_First_Step_Call_GPU(particle_pair * pair_array_gpu, const size_t size)
{


	dim3 block = get_compute_block();;
	dim3 grid = get_grid(size, block);
	pairs_first_step_on_gpu <<< grid, block >>> (pair_array_gpu, size);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());

}

void Prepare_Laser_QQ_array(double* qq_array_gpu)
{

	dim3 block = get_pre_block();;
	dim3 grid = get_grid((2 * two_steps), block);
	pre_second_step_qq << < grid, block >> > (qq_array_gpu);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}

void Prepare_Laser_E1_array(double* qq_array_gpu,double* e1_array_gpu)
{
	dim3 block = get_pre_block();;
	dim3 grid = get_grid((2 * two_steps), block);

	pre_second_step_e1_arr << < grid, block >> > (qq_array_gpu,EE0_Check,e1_array_gpu);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());

}

void Prepare_Laser_E2_array(double* qq_array_gpu,double * e2_array_gpu)
{

	dim3 block = get_pre_block();;
	dim3 grid = get_grid((2 * two_steps), block);

	pre_second_step_e2_arr <<< grid, block >>> (qq_array_gpu, EE0_Check, e2_array_gpu);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}

void Prepare_Laser_E_Check_array(double* e1_array_gpu, double* e2_array_gpu, double* e_check_array_gpu)
{

	dim3 block = get_pre_block();;
	dim3 grid = get_grid((2 * two_steps), block);
	pre_second_step_E_arr_check << < grid, block >> > (e1_array_gpu,e2_array_gpu,e_check_array_gpu);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
}

void Pairs_Second_Step_Once_Call_GPU
(particle_pair * pair_array_first_step_gpu, double* qq_array_gpu, const size_t size, const int index,
 unsigned long long& count_z_once, unsigned long long& count_zz_once)
{
	

	double EE0 = compute_ee0_by_index(index);
	Pairs_Second_Step_Once_Use_E0_Call_GPU(pair_array_first_step_gpu, qq_array_gpu, size, EE0,
		count_z_once, count_zz_once);




}








void Pairs_Second_Step_Filter_Call_GPU
(particle_pair * pair_array_sec_step_gpu, particle_pair * pair_array_filtered_gpu,
 size_t size, unsigned long long& count_z, unsigned long long& count_zz)
{
	count_z = 0;
	count_zz = 0;
	unsigned long long  *gpu_count_z, *gpu_count_zz;
	CHECK(cudaMalloc((void**)&gpu_count_z, size_ull));
	CHECK(cudaMalloc((void**)&gpu_count_zz, size_ull));

	dim3 com_block = get_compute_block();
	dim3 com_grid = get_grid(size, com_block);
	pairs_second_step_on_gpu_fliter << < com_grid, com_block, 0, 0 >> > (pair_array_sec_step_gpu,
	                                                                     pair_array_filtered_gpu, size, gpu_count_z , gpu_count_zz );

	
	
	
	CHECK(cudaMemcpy(&count_z , gpu_count_z ,
		size_ull, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&count_zz,gpu_count_zz,
		size_ull, cudaMemcpyDeviceToHost));

}

void Pairs_Second_Step_Whole_Call_GPU(particle_pair* pair_array_gpu, const size_t size, const int iter_times)
{

	double* qq_array_gpu;
	CHECK(cudaMalloc((void**)&qq_array_gpu, Bytes_Of_Array_Laser));
	Prepare_Laser_QQ_array(qq_array_gpu);

	//保存每次迭代的z,zz
	unsigned long long* z_count_arr = new unsigned long long[iter_times];
	unsigned long long* zz_count_arr = new unsigned long long[iter_times];
	//保存每次迭代的ee0
	double* ee0_arr = new double[iter_times];
	for(int i = 0;i<iter_times;i++)
	{
		ee0_arr[i]= compute_ee0_by_index(i);
		unsigned long long z_once, zz_once;
		Pairs_Second_Step_Once_Call_GPU(pair_array_gpu, qq_array_gpu, size, i,
		                                z_once,zz_once);
		z_count_arr[i] = z_once;
		zz_count_arr[i] = zz_once;
	}

	Print_Count_Array(ee0_arr, z_count_arr, zz_count_arr, iter_times, ion_rate_file_name.c_str());

	CHECK(cudaFree(qq_array_gpu));

	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());


}



void Pairs_Second_Step_Once_Use_E0_Call_GPU(
	particle_pair * pair_array_first_step_gpu, double* qq_array_gpu, const size_t size, double EE0,
	unsigned long long& count_z_once, unsigned long long& count_zz_once)
{
	double *gpu_e1, *gpu_e2;
	CHECK(cudaMalloc((void **)(&gpu_e1), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), Bytes_Of_Array_Laser));

	//double EE0 = compute_ee0_by_index(index);
	//double EE0 = EE0_Check;
	dim3 pre_block = get_pre_block();
	dim3 pre_grid = get_grid((2 * two_steps), pre_block);
	pre_second_step_e1_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e1);
	pre_second_step_e2_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e2);

	//计算第二步循环
	particle_pair* second_array_gpu;
	dim3 com_block = get_compute_block();
	dim3 com_grid = get_grid(size, com_block);
	CHECK(cudaMalloc((void **)(&second_array_gpu), Bytes_Of_Pairs));
	CHECK(cudaMemcpy(second_array_gpu, pair_array_first_step_gpu, Bytes_Of_Pairs, cudaMemcpyDeviceToDevice));
	pairs_second_step_on_gpu << <com_grid, com_block >> > (second_array_gpu, size, gpu_e1, gpu_e2);


	//第二步循环后过滤
	particle_pair* second_array_filter_gpu;
	CHECK(cudaMalloc((void **)(&second_array_filter_gpu), Bytes_Of_Pairs));
	unsigned long long count_z, count_zz;
	Pairs_Second_Step_Filter_Call_GPU(second_array_gpu, second_array_filter_gpu, size, count_z, count_zz);
	count_z_once = count_z;
	count_zz_once = count_zz;

	SavePairsWhichOnGPU(second_array_filter_gpu, count_zz_once, second_step_file_name.c_str());

	//SavePairsWhichOnGPU(second_array_gpu,size,"OneStep.dat");

	CHECK(cudaFree(second_array_gpu));
	CHECK(cudaFree(second_array_filter_gpu));

	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());


}





void Pairs_Second_Step_Once(particle_pair* pair_array_gpu, const size_t size)
{

	double* qq_array_gpu;
	CHECK(cudaMalloc((void**)&qq_array_gpu, Bytes_Of_Array_Laser));
	Prepare_Laser_QQ_array(qq_array_gpu);

	
		
	unsigned long long z_once, zz_once;
	Pairs_Second_Step_Once_Use_E0_Call_GPU(pair_array_gpu, qq_array_gpu, size, EE0_now,
		                                z_once,zz_once);
	
	


	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());


}



void every_step(int pairs)
{

	//申请显存空间
	particle_pair *pairs_array_single_step_gpu;
	CHECK(cudaMalloc((void **)(&pairs_array_single_step_gpu), sizeof(particle_pair)));


	//计时
	double start = seconds();
	//计算
	Pairs_Init_Call_GPU(pairs_array_single_step_gpu, pairs);
	//保存
	SavePairsWhichOnGPU(pairs_array_single_step_gpu, pairs, init_file_name.c_str());
	//初始化完成
	double elapse = seconds();
	printf("Inition compltete %lf\n", elapse - start);




	double *gpu_e1, *gpu_e2, *qq_array_gpu;
	CHECK(cudaMalloc((void **)(&gpu_e1), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&qq_array_gpu), Bytes_Of_Array_Laser));

	//double EE0 = compute_ee0_by_index(index);
	double EE0 = EE0_Check;
	dim3 pre_block = get_pre_block();
	dim3 pre_grid = get_grid((2 * two_steps), pre_block);
	pre_second_step_qq << < pre_grid, pre_block >> > (qq_array_gpu);
	pre_second_step_e1_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e1);
	pre_second_step_e2_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e2);

	particle_pair *pairs_array_every_step_gpu;
	CHECK(cudaMalloc((void **)(&pairs_array_every_step_gpu), sizeof(particle_pair) * two_steps));
	start = seconds();
	//计算 后保存电离率
	//Pairs_Second_Step_Whole_Call_GPU(pairs_array_single_step_gpu, pairs, Iter_Count);
	dim3 block = get_compute_block();
	dim3 grid = get_grid(pairs, block);
	pairs_second_step_on_gpu_every_step << < 1,1 >> > (pairs_array_single_step_gpu, pairs, gpu_e1, gpu_e2, pairs_array_every_step_gpu);
	SavePairsWhichOnGPU(pairs_array_every_step_gpu, two_steps, "every_step.dat");
}
