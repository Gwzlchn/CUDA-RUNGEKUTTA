#ifndef __CALL_GPU_H
#define __CALL_GPU_H

#include "../include/Struct_Defines.h"
#include <vector_types.h>


dim3 get_pre_block(int dimx = 512);
dim3 get_compute_block(int dimx = 32);
dim3 get_grid(size_t size,const dim3& block);


void SaveArraysWhichOnGPU(double* gpu_array, size_t  size, const char* file_name);

void SaveLaserArraysWhichOnGPU(double* e1_array, double* e2_array, double* e_check_array, size_t size, const char* file_name);


void SavePairsWhichOnGPU(particle_pair* gpu_array, size_t  size,const char* file_name);

double* AllocArayOnGPU(size_t size);

particle_pair* AllocPairsOnGPU(size_t size);



//用于双核粒子的随机数化  初始化
void Pairs_Init_Call_GPU(particle_pair* pair_array_gpu, const size_t size);

void Pairs_First_Step_Call_GPU(particle_pair* pair_array_gpu, const size_t size);

void Prepare_Laser_QQ_array(double* qq_array_gpu);

void Prepare_Laser_E1_array(double* qq_array_gpu,double* e1_array_gpu);

void Prepare_Laser_E2_array(double* qq_array_gpu,double* e2_array_gpu);

void Prepare_Laser_E_Check_array(double* e1_array_gpu,double* e2_array_gpu,double* e_check_array);




void Pairs_Second_Step_Once_Call_GPU(particle_pair * pair_array_first_step_gpu, double* qq_array_gpu, const size_t size, const int index,
                                     unsigned long long& count_z_once, unsigned long long& count_zz_once);

void Pairs_Second_Step_Filter_Call_GPU(particle_pair* pair_array_sec_step_gpu, particle_pair* pair_array_filtered, size_t size,
                                       unsigned long long& count_z, unsigned long long& count_zz);


void Pairs_Second_Step_Whole_Call_GPU(particle_pair* pair_array_gpu, const size_t size, const int iter_times);

void Pairs_Second_Step_Once(particle_pair* pair_array_gpu, const size_t size);




#endif //__CALL_GPU_H