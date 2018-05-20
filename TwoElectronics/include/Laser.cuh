#ifndef __LASER_CUH
#define __LASER_CUH

#include "device_launch_parameters.h"
#include "Struct_Defines.h"


__device__ double CalculationE1(const particle& first, const particle& second);
__device__ double CalculationE2(const particle& first, const particle& second);





__device__ double compute_qq_single(const size_t& now_step);

__device__ double compute_e_for_check(const size_t& now_step, const double& e1_single, const double& e2_single);

__device__ double compute_e1_single(const size_t& now_step,const double& qq_now_single, const double& EE0);

 __device__ double compute_e2_single(const size_t& now_step, const double& qq_now_single, const double& EE0);


 __host__ __device__ double compute_ee0_by_index(const int index);

#endif


