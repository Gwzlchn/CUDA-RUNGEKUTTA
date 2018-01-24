#ifndef DEVICE_COMPUTE_FUNCS_CUH
#define DEVICE_COMPUTE_FUNCS_CUH

#include "sci_const.cuh"
#include "nucleus.hpp"

#include <thrust/device_vector.h>




__device__ double fx_first(const nucleus& first, const nucleus& second,double t,double a,double a1 );
__device__ void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec);
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec);

#endif //DEVICE_COMPUTE_FUNCS_CUH

__device__  double fx_first(const nucleus& first, const nucleus& second, double t,double a,double a1)
{
	return (first.x - second.x) / sqrt(pow((pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2) + a1*a1), 3)) -
		(first.x - nuclear_spacing / 2.0*sin(PI*rotation)) / sqrt((pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2)));
}


void update_step_one(nucleus* step_one_fir, nucleus* step_one_sec)
{
	thrust::host_vector<nucleus> K1;
	//step_one_fir->x

}
__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec)


