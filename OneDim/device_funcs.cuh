#ifndef DEVIVE_FUNCS_CUH
#define DEVIVE_FUNCS_CUH

//此处只有device函数。以及科学常量



__device__ const double PI=3.14159265358979323846; 
__device__ const double A=0.05;
__device__ const double E0=-0.5;
__device__ const double omega=0.057;
//步长DX 终点TIME
__device__ const double T0 = 2.0*PI/omega;
__device__ const int 	STEPS=40000;
__device__ const double DX=10.0*T0/STEPS;
__device__ const int 	STEPSFIRST=10000;
__device__ const int	STEPSSECOND=40000;


__device__ double fx(double x);
__device__ double Ekall(double x);
__device__ double Px(double x);

__device__ double EkallAtStepTwo(double t);
__device__ double fxAtStepTwo(double x,double time);

__device__ void updateXi(double& xi,double& pxi);
__device__ void updateXiAtStepTwo(double& xi,double& pxi,double time);
#endif

__device__ double fx(double x)
{

	return -x/sqrt((pow((pow(x,2.0)+pow(A,2.0)),3.0)));
}


__device__ double Ekall(double x)
{

	return E0+1.0/(sqrt(pow(x,2.0)+pow(A,2.0)));
}

__device__ double Px(double x)
{
	return sqrt(2*Ekall(x));
}

__device__ double EkallAtStepTwo(double t)

{
	return A*pow( (sin(omega/2.0/10.0*t)) , 2.0) * sin(omega * t);
}

__device__ double fxAtStepTwo(double x,double time)
{
	return fx(x)-EkallAtStepTwo(time);
}

__device__ void updateXi(double& xi,double& pxi)
{
	
	double K1  = pxi;
	double K11 = fx(xi);
	
	double K2  = pxi + K11/2.0*DX;
	double K22 = fx(xi + K1/2.0*DX);
	
	double K3  = pxi + K22/2.0*DX;
	double K33 = fx(xi + K2/2.0*DX);
	
	double K4  = pxi+K33*DX;
	double K44 = fx(xi+K3*DX);
	
	xi  = xi  + DX * (K1  + 2*K2  + 2*K3  + K4)/6.0;
	pxi = pxi + DX * (K11 + K22*2 + K33*2 + K44)/6.0;
	return;
}



__device__ void updateXiAtStepTwo(double& xi,double& pxi,double time)
{
	
	double K1  = pxi;
	double K11 = fxAtStepTwo(xi,time);
	
	double K2  = pxi + K11/2.0*DX;
	double K22 = fxAtStepTwo(xi + K1/2.0*DX,time+DX/2.0);
	
	double K3  = pxi + K22/2.0*DX;
	double K33 = fxAtStepTwo(xi + K2/2.0*DX,time+DX/2.0);
	
	double K4  = pxi+K33*DX;
	double K44 = fxAtStepTwo(xi + K3*DX,time+DX);
	
	xi  = xi  + DX * (K1  + 2*K2  + 2*K3  + K4)/6.0;
	pxi = pxi + DX * (K11 + K22*2 + K33*2 + K44)/6.0;
	return;
}
