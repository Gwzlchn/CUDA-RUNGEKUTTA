#define E0  1.59 // 基态能量
#define PI 3.1415926535897932384626433832795
#define  E_sec_ion  1.065 // 第二电离能,Ip2

#define   Q_constant    1.225

#define   A_hardness    2.0 //硬度参数

#define   Omega1    0.0584

#define   Omega2    0.117

#define   N1_const    2.0

#define   N2_const    6.0

#define   T0_const    (2 * PI / Omega1)

//下面的10000为one_steps

#define   DX    ((2.0*N1_const + N2_const) * T0_const / (10000))
//#include<iostream>


int main(int argc, char const *argv[])
{
    int dx = DX;

    printf("%.10f\n",T0_const);
    printf("%.10f",dx);

    return 0;
}
