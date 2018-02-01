
#ifndef SCI_CONST_H
#define SCI_CONST_H

#include <crt/host_defines.h>

//随机数有关
__device__ double stddev = 0.7;		//方差

__device__ double mean = 2.0;			//均值


//核粒子有关
__device__ double rotation = 0.0;	//转轴角度，对应之前的kk

__device__ double nuclear_spacing = 4.0;	//核间距,对应之前的R

__device__ double PI = 3.1415926535897932384626433832795;	//圆周率

__device__ double elec_nucl =1.25;	//电子与核之间参数，对应之前A

__device__ double pow_elec_nucl = elec_nucl * elec_nucl;//A*A

__device__ double elec_elec = 0.1;	//电子与电子之间参数。对应之前A1

__device__ double pow_elec_elec = elec_elec * elec_elec; //A1*A1

//!!注意这个4.0应该是R,但windows不支持。服务器可以
__device__ double E_total = -1.0165 - 1.0 / nuclear_spacing;	//两个电子总能量 对应之前 E0



//龙哥库塔用到变量
__device__ int  one_steps = 10000;

__device__ int	two_steps = 40000;

__device__ double DX = 0.0275438596;



#endif
