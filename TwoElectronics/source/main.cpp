#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "../include/Erorr_Check.hpp"
#include "../include/Sci_Constant.h"
#include "../include/Call_GPU.cuh"
#include "../include/PrintStruct.h"
#include  "../include/Compute_On_GPU.cuh"
#include  <cmath>




void every_step(int pairs = 1)
{

	//申请显存空间
	particle_pair *pairs_array_single_step_gpu;
	CHECK(cudaMalloc((void **)(&pairs_array_single_step_gpu), sizeof(particle_pair) ));

	
	//计时
	double start = seconds();
	//计算
	Pairs_Init_Call_GPU(pairs_array_single_step_gpu, pairs);
	//保存
	SavePairsWhichOnGPU(pairs_array_single_step_gpu, pairs, init_file_name.c_str());
	//初始化完成
	double elapse = seconds();
	printf("Inition compltete %lf\n", elapse - start);




	double *gpu_e1, *gpu_e2,*qq_array_gpu;
	CHECK(cudaMalloc((void **)(&gpu_e1), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&gpu_e2), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&qq_array_gpu), Bytes_Of_Array_Laser));

	//double EE0 = compute_ee0_by_index(index);
	double EE0 = EE0_Check;
	dim3 pre_block = get_pre_block();
	dim3 pre_grid = get_grid((2 * two_steps), pre_block);
	pre_second_step_qq <<< pre_grid, pre_block >>> (qq_array_gpu);
	pre_second_step_e1_arr <<< pre_grid, pre_block, 0, 0 >>> (qq_array_gpu, EE0, gpu_e1);
	pre_second_step_e2_arr <<< pre_grid, pre_block, 0, 0 >>> (qq_array_gpu, EE0, gpu_e2);

	particle_pair *pairs_array_every_step_gpu;
	CHECK(cudaMalloc((void **)(&pairs_array_every_step_gpu), sizeof(particle_pair) * two_steps ));
	 start = seconds();
	//计算 后保存电离率
	//Pairs_Second_Step_Whole_Call_GPU(pairs_array_single_step_gpu, pairs, Iter_Count);
	dim3 block = get_compute_block();;
	dim3 grid = get_grid(1, block);
	pairs_second_step_on_gpu_every_step <<< grid, block >>> (pairs_array_single_step_gpu, 1, gpu_e1, gpu_e2, pairs_array_every_step_gpu);
	SavePairsWhichOnGPU(pairs_array_every_step_gpu, two_steps, "every_step.dat");
}




void check_laser_array_on_gpu()
{
	double start = seconds();
	//申请显存空间
	double* qq_arr_gpu, *e1_arr_gpu, *e2_arr_gpu, *e_check_arr_gpu;
	CHECK(cudaMalloc((void **)(&qq_arr_gpu), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&e1_arr_gpu), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&e2_arr_gpu), Bytes_Of_Array_Laser));
	CHECK(cudaMalloc((void **)(&e_check_arr_gpu), Bytes_Of_Array_Laser));
	
	size_t array_size = 2 * two_steps;
	
	//-------------------计算QQ------------------------
	Prepare_Laser_QQ_array(qq_arr_gpu);
	//保存
	SaveArraysWhichOnGPU(qq_arr_gpu, array_size, qq_array_file_name.c_str());

	//---------------计算E1 E2 E_Check-------------------
	Prepare_Laser_E1_array(qq_arr_gpu, e1_arr_gpu);
	Prepare_Laser_E2_array(qq_arr_gpu, e2_arr_gpu);
	Prepare_Laser_E_Check_array(e1_arr_gpu, e2_arr_gpu, e_check_arr_gpu);
	SaveLaserArraysWhichOnGPU(e1_arr_gpu, e2_arr_gpu, e_check_arr_gpu, array_size, laser_array_file_name.c_str());



	CHECK(cudaFree(qq_arr_gpu));
	CHECK(cudaFree(e1_arr_gpu));
	CHECK(cudaFree(e2_arr_gpu));
	CHECK(cudaFree(e_check_arr_gpu));

	double elapse = seconds();
	printf("Check compltete %lf\n", elapse - start);


}



void compute_on_gpu_all(size_t pairs)
{
	//申请显存空间
	particle_pair *pairs_array_single_step_gpu;
	CHECK(cudaMalloc((void **)(&pairs_array_single_step_gpu), Bytes_Of_Pairs));
	

	//-------------------------初始化---------------------

	//计时
	double start = seconds();
	//计算
	Pairs_Init_Call_GPU(pairs_array_single_step_gpu, pairs);
	//保存
	SavePairsWhichOnGPU(pairs_array_single_step_gpu, pairs, init_file_name.c_str());
	//初始化完成
	double elapse = seconds();
	printf("Inition compltete %lf\n", elapse - start);

	//----------------------------第一步--------------------------
	//计时
	//start = seconds();
	//计算
	//Pairs_First_Step_Call_GPU(pairs_array_single_step_gpu, pairs);
	//保存
	//SavePairsWhichOnGPU(pairs_array_single_step_gpu, pairs, first_file_name.c_str());
	//第一步完成
	//elapse = seconds();
	//printf("FirstStep compltete %lf\n", elapse - start);

	//---------------------------第二步----------------------------
	//计时
	start = seconds();
	//计算 后保存电离率
	//Pairs_Second_Step_Whole_Call_GPU(pairs_array_single_step_gpu, pairs, Iter_Count);
	Pairs_Second_Step_Once(pairs_array_single_step_gpu, pairs);
	//第二步完成
	elapse = seconds();
	printf("SecondStep compltete %lf\n", elapse - start);


	//整理工作 释放内存
	CHECK(cudaFree(pairs_array_single_step_gpu));


	
}






int main()
{
	

	printf("Starting...\n");

	//选择设备
	int dev = 0;
	cudaDeviceProp deviceProp;
	
	CHECK(cudaGetDeviceProperties(&deviceProp, dev));
	printf("Using Device %d: %s\n", dev, deviceProp.name);
	CHECK(cudaSetDevice(dev));


	//compute_on_gpu_all(Pairs_Total);
	//check_laser_array_on_gpu();
	every_step();
	return 0;
}



