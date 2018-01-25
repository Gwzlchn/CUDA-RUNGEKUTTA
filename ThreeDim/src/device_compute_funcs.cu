
#include <curand.h>
#include <curand_kernel.h>
#include <cmath>
#include <vector_types.h>
#include "../include/sci_const.h"
#include "../include/device_compute_funcs.h"

__device__ double E_kall(const nucleus& first, const nucleus& second)
{
	return E_total - (-1.0 / sqrt(pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
		pow((first.x - nuclear_spacing - 2.0*sin(PI*rotation)), 2) +
		first.y*first.y + elec_elec*elec_elec)) -
		(-1.0 / sqrt(pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec))
		- (1.0 / sqrt(nucleus_distance(first, second) + elec_nucl*elec_nucl)) -
		(-1.0 / sqrt(pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec)) -
			(-1.0 / sqrt(pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
				pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
				second.y*second.y + elec_elec*elec_elec));
}

__device__ void px_py_pz_distribution(nucleus& first, nucleus& second)
{
	double ekall = E_kall(first, second);
	curandStatePhilox4_32_10_t s;
	unsigned long long seed = 1;
	// seed a random number generator 
	curand_init(seed, 0, 0, &s);
	double2 random12 = curand_uniform2_double(&s);
	double2 random34 = curand_uniform2_double(&s);
	double random5 = curand_uniform_double(&s);

	double theta1 = random12.x*PI;
	double theta2 = random12.y*PI;
	double phi1 = random34.x * 2 * PI;
	double phi2 = random34.y * 2 * PI;

	first.px = sqrt(2.0*ekall*random5)*sin(theta1)*cos(phi1);
	first.py = sqrt(2.0*ekall*random5)*sin(theta1)*cos(phi1);
	first.pz = sqrt(2.0*ekall*random5)*cos(phi1);

	second.px = sqrt(2.0*ekall*(1 - random5))*sin(theta2)*cos(phi2);
	second.py = sqrt(2.0*ekall*(1 - random5))*sin(theta2)*cos(phi2);
	second.pz = sqrt(2.0*ekall*(1 - random5))*cos(phi2);



}


__device__  double nucleus_distance(const nucleus& first, const nucleus& second)
{
	return (pow((first.x - second.x), 2) + pow((first.y - second.y), 2) + pow((first.z - second.z), 2));
}




__device__  nucleus fx_first_nucleus(const nucleus& first, const nucleus& second)
{
	nucleus fx_first;

	fx_first.x = (first.x - second.x)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (first.x - nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- (first.x + nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.y = (first.y - second.y)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- first.y
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- first.y
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0*sin(PI*rotation)), 3) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.z = (first.z - second.z)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (first.z - nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((first.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3))
		- (first.z + nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((first.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((first.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			first.y*first.y + elec_elec*elec_elec), 3));

	fx_first.px = fx_first.py = fx_first.pz = 0;
	return fx_first;
}
__device__  nucleus fx_second_nucleus(const nucleus& first, const nucleus& second)
{
	nucleus fx_second;
	fx_second.x = (second.x - first.x)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (second.x - nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- (second.x + nuclear_spacing / 2.0 * sin(PI*rotation))
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3));

	fx_second.y = (second.y - first.y)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- second.y
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0*sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- second.y
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0*cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0*sin(PI*rotation)), 3) +
			second.y*second.y + elec_elec*elec_elec), 3));

	fx_second.z = (first.z - second.z)
		/ sqrt(pow((nucleus_distance(first, second) + elec_nucl*elec_nucl), 3))
		- (second.z - nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((second.z - nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x - nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3))
		- (second.z + nuclear_spacing / 2.0 * cos(PI*rotation))
		/ sqrt(pow((pow((second.z + nuclear_spacing / 2.0 * cos(PI*rotation)), 2) +
			pow((second.x + nuclear_spacing / 2.0 * sin(PI*rotation)), 2) +
			second.y*second.y + elec_elec*elec_elec), 3));
	fx_second.px = fx_second.py = fx_second.pz = 0;
	return fx_second;
}


void update_step_one(nucleus& step_one_fir, nucleus& step_one_sec)
{


}
//__device__ void update_step_two(nucleus* step_two_fir, nucleus* step_two_sec)

//生成双精度双正态分布随机数
__global__ void DoubleNormalRandomArrayD(nuclei* Array, const long Size)
{
	double A1, A2, A3, A4, Ekall;
	int i = threadIdx.x;
	double temp1 = 1;
	double temp2 = 1;
	curandState s;

	Ekall = -1;

	while (Ekall < 0)
	{
		A2 = A4 = 2;

		while (A2 > temp1 && A4 > temp2)
		{
			A1 = curand_uniform_double(&s);
			A2 = curand_uniform_double(&s);
			A3 = curand_uniform_double(&s);
			A4 = curand_uniform_double(&s);

			A1 = (A1 - 0.5) * 20;
			A3 = (A3 - 0.5) * 20;

			temp1 = exp((-pow((A1 - nuclear_spacing / 2.0), nuclear_spacing / 2.0)) /
				(nuclear_spacing / 2.0 * pow(stddev, nuclear_spacing / 2.0)))
				+ exp((-pow((A1 + nuclear_spacing / 2.0), nuclear_spacing / 2.0)) /
				(nuclear_spacing / 2.0 * pow(stddev, nuclear_spacing / 2.0)));
			temp2 = exp((-pow((A3 - nuclear_spacing / 2.0), nuclear_spacing / 2.0)) /
				(nuclear_spacing / 2.0 * pow(stddev, nuclear_spacing / 2.0)))
				+ exp((-pow((A3 + nuclear_spacing / 2.0), nuclear_spacing / 2.0)) /
				(nuclear_spacing / 2.0 * pow(stddev, nuclear_spacing / 2.0)));
		}

		Array[i].first.x = A1 * sin(rotation*PI);
		Array[i].first.y = 0;
		Array[i].first.z = A1 * cos(rotation*PI);
		Array[i].second.x = A3 * sin(rotation*PI);
		Array[i].second.y = 0;
		Array[i].second.z = A3 * cos(rotation*PI);
		Ekall = E_kall(Array[i].first, Array[i].second);
	}
	return;
}

//用于双核粒子的随机数化
extern "C" void NucleiRandomD(nuclei* Array, const long Size)
{
	int threadsPerBlock = 256;
	int threadsPerGrid = (2 * Size + threadsPerBlock - 1) / threadsPerBlock;
	DoubleNormalRandomArrayD << <threadsPerGrid, threadsPerBlock >> > (Array, Size);
}