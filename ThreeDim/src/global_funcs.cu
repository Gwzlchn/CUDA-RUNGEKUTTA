#include "../include/global_funcs.h"
#include<curand.h>

void NormalRandomArray(double* arr, const long n2)
{
	curandGenerator_t gen;                                  //生成随机数变量
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);//步骤1：指定算法
	curandSetPseudoRandomGeneratorSeed(gen, 11ULL);         //步骤2：随机数初始化
	curandGenerateNormalDouble(gen, arr, n2, 0, 0.7);        //步骤3：生成随机数，存储到缓冲器中（第1个数字为均值，第二个为方差）
	curandDestroyGenerator(gen);                         	//释放指针
	return;
}
void NormalRandomNuclei(nuclei* raw_nuclei,double* random_arr ,const long n)
{
	for(long i=0;i<(n/2);i++)
	{
		raw_nuclei[i].init_first.x = random_arr[i];
		raw_nuclei[i].init_second.x = random_arr[i + (n / 2)];
	}
}