#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "../include/Erorr_Check.hpp"
#include "../include/Sci_Constant.h"
#include "../include/Call_GPU.cuh"
#include "../include/PrintStruct.h"



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
	start = seconds();
	//计算
	Pairs_First_Step_Call_GPU(pairs_array_single_step_gpu, pairs);
	//保存
	SavePairsWhichOnGPU(pairs_array_single_step_gpu, pairs, first_file_name.c_str());
	//第一步完成
	elapse = seconds();
	printf("FirstStep compltete %lf\n", elapse - start);

	//---------------------------第二步----------------------------
	//计时
	start = seconds();
	//计算 后保存电离率
	Pairs_Second_Step_Whole_Call_GPU(pairs_array_single_step_gpu, pairs, Iter_Count);
	
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
	check_laser_array_on_gpu();
	return 0;
}



