#include"./test_define.cpp"
#include"./support/PrintStruct.h"
#include<cmath>

int main()
{
    double* QQ = new double[2*two_steps_in_host];
    double* t1_arr = new double[2*two_steps_in_host];

    for(int idx=0;idx!=2*two_steps_in_host;idx++){
        double t1 = 0.5 * DX * (idx+1);
        t1_arr[idx]=t1;
	    QQ[idx] = pow((sin(Omega1 / 2.0 / (2 * N1_const + N2_const)*t1)), 2);
    }

    PrintArray(t1_arr,2*two_steps_in_host,"t1_check",1);
    PrintArray(QQ,2*two_steps_in_host,"qq_cpu_check",1);
    return 0;
	    	
}