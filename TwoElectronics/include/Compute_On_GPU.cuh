#ifndef __CALL_CPU_CUH
#define __CALL_GPU_CUH

#include "./Struct_Defines.h"
#include "device_launch_parameters.h"



__global__ void pairs_init(particle_pair* pair_array, const long size,
	const double min_r, const double min_p);

__global__ void  pairs_first_step_on_gpu(particle_pair* first_setp_pair_array, const long size);


__global__ void pre_second_step_qq(double* QQ_array);


__global__ void pre_second_step_E_arr_check
(const double* E1_array, const double* E2_array, double* E_check_array);


__global__ void pre_second_step_e1_arr(const double* QQ_array, const double EE0, double* E1_array);


__global__ void pre_second_step_e2_arr(const double* QQ_array, const double EE0, double* E2_array);

__global__ void pairs_second_step_on_gpu
(particle_pair* second_arr, const long size, double* E1_array, double* E2_array);


__global__ void pairs_second_step_on_gpu_fliter
(const particle_pair* second_step_pair_array, particle_pair* second_step_pair_array_filter,
	const long size, unsigned long long* count_z, unsigned long long* count_zz);


#endif //CALL_CPU_CUH