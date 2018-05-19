#ifndef __CALL_GPU_H
#define __CALL_GPU_H

#include "../include/Struct_Defines.h"
#include <vector_types.h>


dim3 get_pre_block(int dimx = 512);
dim3 get_compute_block(int dimx = 32);
dim3 get_grid(long size,const dim3& block);


void SaveArraysWhichOnGPU(double* array, long int size, char* file_name);

void SavePairsWhichOnGPU(particle_pair* array, long int size, char* file_name);

double* AllocArayOnGPU(long int size);

particle_pair* AllocPairsOnGPU(long int size);



//用于双核粒子的随机数化  初始化
void Pairs_Init_Call_GPU(particle_pair* pair_array_gpu, const long size);

void Pairs_First_Steo_Call_GPU(particle_pair* pair_array_gpu, const long size);


void Prepare_Laser_QQ_array(double* qq_array_gpu);

void Prepare_Laser_One_E1_array(double* e1_array_gpu, const double* qq_array_gpu, const int laser_index);

void Prepare_Laser_One_E2_array(double* e2_array_gpu, const double* qq_array_gpu, const int laser_index);

void Pairs_Second_Step_Call_GPU(particle_pair* pair_array_gpu, const long size);



#endif //__CALL_GPU_H