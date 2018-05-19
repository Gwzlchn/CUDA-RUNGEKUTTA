#ifndef __CALL_GPU_H
#define __CALL_GPU_H

#include "../include/Struct_Defines.h"


void SaveArraysWhichOnGPU(double* array, long int size, char* file_name);

void SavePairsWhichOnGPU(particle_pair* array, long int size, char* file_name);

double* AllocArayOnGPU(long int size);

particle_pair* AllocPairsOnGPU(long int size);







#endif //__CALL_GPU_H