#ifndef __LASER_CUH
#define __LASER_CUH

#include "device_launch_parameters.h"
#include "Struct_Defines.h"


__device__ double CalculationE1(const nucleus& first, const nucleus& second);
__device__ double CalculationE2(const nucleus& first, const nucleus& second);



#endif


