#include "../include/Call_GPU.cuh"
#include "../include/Erorr_Check.hpp"
#include "../include/Init_First_Second.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Compute_On_GPU.cuh"
#include "../include/Laser.cuh"
#include <cstdio>
#include <cuda_runtime_api.h>

dim3 get_pre_block(int dimx)
{
	return dim3(dimx);
}

dim3 get_compute_block(int dimx)
{
	return dim3(dimx);
}

dim3 get_grid(long size, const dim3& block)
{
	return dim3((size + block.x - 1) / block.x, 1);
}





void Pairs_Init_Call_GPU(particle_pair * pair_array_gpu, const long size)
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

void Pairs_First_Steo_Call_GPU(particle_pair * pair_array_gpu, const long size)
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

void Pairs_Second_Step_Once_Call_GPU(particle_pair * pair_array_gpu,double* qq_array_gpu, const long size,const int index)
{
	double *gpu_e1, *gpu_e2;
	CHECK(cudaMalloc((void **)(&gpu_e1), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), Bytes_Of_Array_Laser));

	double EE0 = compute_ee0_by_index(index);
	dim3 pre_block = get_pre_block();
	
	dim3 pre_grid = get_grid((2*two_steps),pre_block);
	pre_second_step_e1_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e1);
	pre_second_step_e2_arr << < pre_grid, pre_block, 0, 0 >> > (qq_array_gpu, EE0, gpu_e2);




}




