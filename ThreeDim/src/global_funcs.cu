#include "../include/global_funcs.h"

void NormalRandomNuclei(nuclei* raw_nuclei,double* random_arr ,const long n)
{
	for(long i=0;i<(n/2);i++)
	{
		raw_nuclei[i].init_first.x = random_arr[i];
		raw_nuclei[i].init_second.x = random_arr[i + (n / 2)];
	}
}

void NormalRandom(nuclei* raw_nuclei, const long n)
{

}
void InitialNuclei(nuclei* randomed_nuclei, const long raw_count, long & left)
{

}
void FirstStep(nuclei* inited_nuclei, const long n)
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
	if(ix<n)
	{
		
	}
}