#include "../include/Compute_On_GPU.cuh"
#include "../include/Sci_Constant.h"
#include "../include/Init_First_Second.cuh"
#include "../include/Runge_Kutta.cuh"
#include "../include/Laser.cuh"

#include <cstdlib>
#include <cuda_runtime.h>


__global__ void pairs_init(particle_pair* pair_array, const long size,
                           const double min_r, const double min_p)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size)
	{

		distribution(pair_array[idx].first, pair_array[idx].second, idx, min_r, min_p);
	}
	return;
}


__global__ void pairs_first_step_on_gpu(particle_pair* first_setp_pair_array, const long size)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx<size)
	{
		for (int i = 0; i < one_steps; i++)
			update_step_one(first_setp_pair_array[idx].first, first_setp_pair_array[idx].second);
	}


}

__global__ void pre_second_step_qq_arr(double * QQ_array)
{
	unsigned long idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		QQ_array[idx] = compute_qq_single(idx);
	}
}




__global__ void pre_second_step_E_arr_check
(const double* E1_array, const double* E2_array, double* E_check_array)
{
	unsigned long idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		E_check_array[idx] = compute_e_for_check(idx, E1_array[idx], E2_array[idx]);
	}
}



__global__ void pre_second_step_e1_arr(const double* QQ_array, const double EE0, double* E1_array)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		E1_array[idx] = compute_e1_single(idx, QQ_array[idx], EE0);
	}

}


__global__ void pre_second_step_e2_arr(const double* QQ_array, const double EE0, double* E2_array)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < 2 * two_steps)
	{
		E2_array[idx] = compute_e2_single(idx, QQ_array[idx], EE0);
	}

}




__global__ void pairs_second_step_on_gpu
(particle_pair* second_arr, const long size, double* E1_array, double* E2_array)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	double e1_laser_t1 = 0.0, e1_laser_t2 = 0.0, e1_laser_t3 = 0.0, e1_laser_t4 = 0.0;
	double e2_laser_t1 = 0.0, e2_laser_t2 = 0.0, e2_laser_t3 = 0.0, e2_laser_t4 = 0.0;
	int idx_of_ds = -1; // 相当于nn
	double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
	double now_t = 0.0; //当前时间，相当于t(1)
	if (idx<size)
	{
		for (int i = 0; i < two_steps; i++)
		{
			//第一个激光场强度
			t1 = now_t;
			if (t1 == 0)
			{
				e1_laser_t1 = 0.0;
				e2_laser_t1 = 0.0;
			}
			else
			{
				idx_of_ds = (2.0 * t1) / DX - 1;
				e1_laser_t1 = E1_array[idx_of_ds];
				e2_laser_t1 = E2_array[idx_of_ds];
			}
			//第二个激光场强度
			t2 = now_t + DX / 2.0;
			idx_of_ds = 2.0 * t2 / DX - 1;
			e1_laser_t2 = E1_array[idx_of_ds];
			e2_laser_t2 = E2_array[idx_of_ds];
			//第三个激光场强度
			t3 = now_t + DX / 2.0;
			idx_of_ds = 2 * t3 / DX - 1;
			e1_laser_t3 = E1_array[idx_of_ds];
			e2_laser_t3 = E2_array[idx_of_ds];
			//第四个激光场强度
			t4 = now_t + DX;
			idx_of_ds = 2.0 * t4 / DX - 1;
			e1_laser_t4 = E1_array[idx_of_ds];
			e2_laser_t4 = E2_array[idx_of_ds];
			double4 e1_laser = make_double4(e1_laser_t1, e1_laser_t2, e1_laser_t3, e1_laser_t4);
			double4 e2_laser = make_double4(e2_laser_t1, e2_laser_t2, e2_laser_t3, e2_laser_t4);
			update_step_two(second_arr[idx].first, second_arr[idx].second,
			                e1_laser, e2_laser);
			now_t = now_t + DX;
			/*if(idx_of_ds == -1 )
			update_step_two(second_arr[idx].first, second_arr[idx].second,
			0.0,DS[0],DS[0],DS[1]);
			else
			{
			update_step_two(second_arr[idx].first, second_arr[idx].second,
			DS[idx_of_ds], DS[idx_of_ds + 1], DS[idx_of_ds + 1], DS[idx_of_ds + 2]);
			}
			idx_of_ds += 2;*/

		}


	}
}



__global__ void pairs_second_step_on_gpu_fliter
(const particle_pair* second_step_pair_array, particle_pair* second_step_pair_array_filter,
 const long size, unsigned long long* count_z, unsigned long long* count_zz)
{
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size)
	{
		double ee1 = CalculationE1(second_step_pair_array[idx].first, second_step_pair_array[idx].second);
		double ee2 = CalculationE2(second_step_pair_array[idx].first, second_step_pair_array[idx].second);

		if (ee1*ee2 < 0)
		{
			atomicAdd(count_z, 1);
		}
		if ((ee1 > 0) && (ee2 > 0))
		{
			unsigned long long temp_idx = atomicAdd(count_zz, 1);
			/*nuclei temp;
			temp.first = second_arr[idx].first;
			temp.second = second_arr[idx].second;
			second__arr_filter[temp_idx - 1] = temp;*/
		}
	}

}
