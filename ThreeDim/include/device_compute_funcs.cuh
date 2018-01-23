#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH
#include <crt/host_defines.h>
#include "nucleus.hpp"
#include <thrust/device_vector.h>




__device__ double fx_step_one();
__device__ void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec);
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec);

#endif //DEVICE_COMPUTE_FUNCS_CUH

__device__ void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec)
{
	thrust::host_vector<nucleus> K1;

}
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec)