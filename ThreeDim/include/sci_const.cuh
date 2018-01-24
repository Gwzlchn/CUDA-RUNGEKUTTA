#pragma once
#ifndef SCI_CONST_CUH
#define SCI_CONST_CUH
#include <crt/host_defines.h>
//这个文件放科学常数


__global__ double rotation = 0.0;//转轴角度，对应之前的kk

__global__ double nuclear_spacing = 4.0;//核间距,对应之前的R

__global__ double PI = 3.1415926535897932384626433832795;//圆周率

__global__ double elec_nucl = 0.1;//电子与核之间参数，对应之前A1

__global__ double elec_elec = 1.25;//电子与电子之间参数。对应之前A

#endif
