
#ifndef SCI_CONST_H
#define SCI_CONST_H

#include <crt/host_defines.h>

__device__ double stddev = 0.7;//方差

__device__ double rotation = 0.0;//转轴角度，对应之前的kk

__device__ double nuclear_spacing = 4.0;//核间距,对应之前的R

__device__ double PI = 3.1415926535897932384626433832795;//圆周率

__device__ double elec_nucl = 0.1;//电子与核之间参数，对应之前A1

__device__ double elec_elec = 1.25;//电子与电子之间参数。对应之前A

__device__ double E_total = -1.0165 - 1.0 / 4.0;  //两个电子总能量 对应之前 E0



#endif
