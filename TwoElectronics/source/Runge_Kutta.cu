#include "../include/Runge_Kutta.cuh"
#include "../include/Init_First_Second.cuh"
#include "../include/Sci_Constant.h"






//第一个粒子 K1~K4 第一步循环
__device__ derivative fisrt_k_one_to_four_fisrt_step(const particle& first, const particle& second)
{
	//二阶导 三个数
	const double3 first_fx = fx_fy_fz_first_nucleus(first, second);
	//一阶导 三个数
	const double3 first_gx = gx_gy_gz_first_nucleus(first, second);

	derivative first_px_fx;
	first_px_fx.px = first_gx.x;
	first_px_fx.py = first_gx.y;
	first_px_fx.pz = first_gx.z;

	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y;
	first_px_fx.fz = first_fx.z;

	return first_px_fx;

}


__device__ derivative second_k_one_to_four_fisrt_step(const particle& first, const particle& second)
{
	//二阶导 三个数
	double3 second_fx = fx_fy_fz_second_nucleus(first, second);
	//一阶导 三个数
	double3 second_gx = gx_gy_gz_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second_gx.x;
	second_px_fx.py = second_gx.y;
	second_px_fx.pz = second_gx.z;

	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y;
	second_px_fx.fz = second_fx.z;
	return second_px_fx;

}



__device__ particle first_and_second_k_add_dx_raw(const derivative& k_one_to_four, const particle& raw_to_add)
{
	double now_dx = DX;

	particle k_add;
	k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	return k_add;
}

__device__ particle first_and_second_k_add_dx_div(const derivative& k_one_to_four, const particle& raw_to_add)
{
	double now_dx = DX / 2.0;

	particle k_add;
	k_add.x = raw_to_add.x + now_dx * k_one_to_four.px;
	k_add.y = raw_to_add.y + now_dx * k_one_to_four.py;
	k_add.z = raw_to_add.z + now_dx * k_one_to_four.pz;
	k_add.px = raw_to_add.px + now_dx * k_one_to_four.fx;
	k_add.py = raw_to_add.py + now_dx * k_one_to_four.fy;
	k_add.pz = raw_to_add.pz + now_dx * k_one_to_four.fz;

	return k_add;
}






__device__ void k_one_to_four_add(const derivative& K1, const derivative& K2, const derivative& K3, const derivative& K4,
	particle& raw_to_add)
{
	raw_to_add.x = raw_to_add.x + DX * (K1.px + 2.0*K2.px + 2.0*K3.px + K4.px) / 6.0;
	raw_to_add.y = raw_to_add.y + DX * (K1.py + 2.0*K2.py + 2.0*K3.py + K4.py) / 6.0;
	raw_to_add.z = raw_to_add.z + DX * (K1.pz + 2.0*K2.pz + 2.0*K3.pz + K4.pz) / 6.0;
	raw_to_add.px = raw_to_add.px + DX * (K1.fx + 2.0*K2.fx + 2.0*K3.fx + K4.fx) / 6.0;
	raw_to_add.py = raw_to_add.py + DX * (K1.fy + 2.0*K2.fy + 2.0*K3.fy + K4.fy) / 6.0;
	raw_to_add.pz = raw_to_add.pz + DX * (K1.fz + 2.0*K2.fz + 2.0*K3.fz + K4.fz) / 6.0;



	return;
}





__device__ void update_step_one(particle& step_one_first, particle& step_one_second)
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	const derivative second_k1 = second_k_one_to_four_fisrt_step(step_one_first, step_one_second);
	const particle first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	const particle second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	const derivative second_k2 = second_k_one_to_four_fisrt_step(first_k1_add, second_k1_add);
	const particle first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	const particle second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	const derivative second_k3 = second_k_one_to_four_fisrt_step(first_k2_add, second_k2_add);
	const particle first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	const particle second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);
	const derivative second_k4 = second_k_one_to_four_fisrt_step(first_k3_add, second_k3_add);

	k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);

	return;
}





//第一个粒子 K1~K4 第二步循环
__device__ derivative fisrt_k_one_to_four_second_step
(const particle& first, const particle& second, const double& e1_laser, const double& e2_laser)
{
	const double3 first_fx = fx_fy_fz_first_nucleus(first, second);
	const double3 first_gx = gx_gy_gz_first_nucleus(first, second);
	derivative first_px_fx;
	first_px_fx.px = first_gx.x;
	first_px_fx.py = first_gx.y;
	first_px_fx.pz = first_gx.z;
	first_px_fx.fx = first_fx.x;
	first_px_fx.fy = first_fx.y - e2_laser;
	first_px_fx.fz = first_fx.z - e1_laser;

	return first_px_fx;

}

//第二个粒子 K1~K4 第二步循环
__device__ derivative second_k_one_to_four_second_step
(const particle& first, const particle& second, const double& e1_laser, const double& e2_laser)
{
	const double3 second_fx = fx_fy_fz_second_nucleus(first, second);
	const double3 second_gx = gx_gy_gz_second_nucleus(first, second);
	derivative second_px_fx;
	second_px_fx.px = second_gx.x;
	second_px_fx.py = second_gx.y;
	second_px_fx.pz = second_gx.z;
	second_px_fx.fx = second_fx.x;
	second_px_fx.fy = second_fx.y - e2_laser;
	second_px_fx.fz = second_fx.z - e1_laser;
	return second_px_fx;

}







__device__ void update_step_two(particle& step_one_first, particle& step_one_second,
	const double4 e1_laser_now, const double4 e2_laser_now)
{
	//计算K1
	const derivative first_k1 = fisrt_k_one_to_four_second_step(step_one_first, step_one_second, e1_laser_now.x, e2_laser_now.x);
	const derivative second_k1 = second_k_one_to_four_second_step(step_one_first, step_one_second, e1_laser_now.x, e2_laser_now.x);
	const particle first_k1_add = first_and_second_k_add_dx_div(first_k1, step_one_first);
	const particle second_k1_add = first_and_second_k_add_dx_div(second_k1, step_one_second);

	//K2
	const derivative first_k2 = fisrt_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	const derivative second_k2 = second_k_one_to_four_second_step(first_k1_add, second_k1_add, e1_laser_now.y, e2_laser_now.y);
	const particle first_k2_add = first_and_second_k_add_dx_div(first_k2, step_one_first);
	const particle second_k2_add = first_and_second_k_add_dx_div(second_k2, step_one_second);

	//K3
	const derivative first_k3 = fisrt_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	const derivative second_k3 = second_k_one_to_four_second_step(first_k2_add, second_k2_add, e1_laser_now.z, e2_laser_now.z);
	const particle first_k3_add = first_and_second_k_add_dx_raw(first_k3, step_one_first);
	const particle second_k3_add = first_and_second_k_add_dx_raw(second_k3, step_one_second);

	//K4
	const derivative first_k4 = fisrt_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);
	const derivative second_k4 = second_k_one_to_four_second_step(first_k3_add, second_k3_add, e1_laser_now.w, e2_laser_now.w);

	k_one_to_four_add(first_k1, first_k2, first_k3, first_k4, step_one_first);
	k_one_to_four_add(second_k1, second_k2, second_k3, second_k4, step_one_second);
}

