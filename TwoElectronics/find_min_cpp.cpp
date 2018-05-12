#include<iostream>
#include<cmath>
using namespace std;

int main()
{
    int nx=351, ny=301;
    double Q_constant = 1.225;
    double A_hardness = 2.0;
    double *P_Arr, *R_Arr;
	R_Arr = (double*)malloc(nx * sizeof(double));
	P_Arr = (double*)malloc(ny * sizeof(double));

	for (int i = 0; i < nx; i++)
		R_Arr[i] = 0.5 + 0.01 * i;
	for (int i = 0; i < ny; i++)
		P_Arr[i] = 0.0 + 0.01*i;


	double** mat = new double*[nx];
	for (int i = 0; i<nx; i++)
		mat[i] = new double[ny];

	double Vh, Vk, Ek;
	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			Vh = pow(Q_constant, 2) / (4.0*A_hardness*pow(R_Arr[i], 2)) *
				exp(A_hardness * (1.0 - pow((R_Arr[i] * P_Arr[j] / Q_constant), 4)));
			Vk = -2.0 / R_Arr[i];
			Ek = P_Arr[j] * P_Arr[j] / 2.0;
			mat[i][j] = Vh + Vk + Ek + 1.065;
		}
	}

	int min_x_index, min_y_index;
	double min = mat[0][0];
	for(int i = 0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			if( min > mat[i][j])
			{
				min = mat[i][j];
				min_x_index = i;
				min_y_index = j;
			}
		}
	}
    double min_r,min_p;
	min_r = R_Arr[min_x_index];
	min_p = P_Arr[min_y_index];

    printf("min_value %.10f, pos x: %d ,  pos y: %d\n",min,min_x_index,min_y_index);
    printf("min_r : %.10f \t min_p: %.10f \n",min_r,min_p);
    return 0;
}