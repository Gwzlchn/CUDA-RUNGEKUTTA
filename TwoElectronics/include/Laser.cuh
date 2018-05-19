#ifndef __LASER_CUH
#define __LASER_CUH

#include "device_launch_parameters.h"
#include "Struct_Defines.h"


__device__ double CalculationE1(const particle& first, const particle& second);
__device__ double CalculationE2(const particle& first, const particle& second);



#endif


