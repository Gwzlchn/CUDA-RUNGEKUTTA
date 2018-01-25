#include "./include/nucleus.hpp"
#include "./include/PrintStruct.hpp"
#include "./include/common.hpp"
#include "./include/Random.hpp"
#include <cuda_runtime.h>
#include <stdio.h>
extern "C" void NucleiRandomD(nuclei* Array, const long Size, double Angle);
int main()
{
	//计时
	double Start, Elaps;

	printf("Starting...\n");

	//选择设备
	int dev = 0;
	cudaDeviceProp deviceProp;
	CHECK(cudaGetDeviceProperties(&deviceProp, dev));
	printf("Using Device %d: %s\n", dev, deviceProp.name);
	CHECK(cudaSetDevice(dev));

	long pairs = 10000;
	long long nBytes = pairs * sizeof(nuclei);
	printf("Use %lld Bytes %lfMB\n", nBytes, nBytes / double(1024 * 1024));
	nuclei* test;
	nuclei* host;
	Start = seconds();
	cudaMalloc((void **)(&test), nBytes);
	host = (nuclei*)malloc(nBytes);
	
	NucleiRandomD(test, pairs, rotation);
	cudaMemcpy(host, test, nBytes, cudaMemcpyDeviceToHost);
	PrintStruct(host, pairs, 0);

	return 0;
	










}



