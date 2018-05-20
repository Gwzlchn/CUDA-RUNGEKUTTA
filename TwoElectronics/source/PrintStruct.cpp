﻿#include <cstdio>
#include <cstring>
#include <cstdio>
#include <cstring>
#include <ctime>

#include "../include/PrintStruct.h"
#include "../include/Sci_Constant.h"
#define BUFLEN 255   

void PrintStruct(particle_pair* ToSaveNuclei, size_t n, const char* FileName)
{
	//format date and time. 
	struct tm *ptr;
	time_t lt;
	char str[9];
	lt = time(NULL);
	ptr = localtime(&lt);
	strftime(str, 9, "%m%d%H%M", ptr);

	FILE* file;
	switch (choose)
	{
		case 0:
		{
			char init_file_name[15] = "TestData/init_";
			///*strcat((char*)(FileName), tmpBuf);
			strcat(init_file_name, FileName);
			strcat(init_file_name, str);
			strcat(init_file_name, ".dat");
			file = fopen(init_file_name, "w");
			break; 
		}
		case 1:
		{
			char step_one_file_name[20] = "TestData/step_one_";
			strcat(step_one_file_name, FileName);
			strcat(step_one_file_name, str);
			strcat(step_one_file_name, ".dat");
			file = fopen(step_one_file_name, "w");
			break;  
		}
		case 2:
		{
			char step_two_file_name[23] = "TestData/step_two_ALL_";
			strcat(step_two_file_name, FileName);
			strcat(step_two_file_name, str);
			strcat(step_two_file_name, ".dat");
			file = fopen(step_two_file_name, "w");
			break;
		}
		case 3:
		{
			char step_two_file_name[26] = "TestData/step_two_fliter_";
			strcat(step_two_file_name, FileName);
			strcat(step_two_file_name, str);
			strcat(step_two_file_name, ".dat");
			file = fopen(step_two_file_name, "w");
			break;
		}

	}
	

	if (!file) perror("cannot open file");
	for (size_t i = 0; i < n; i++)
	{
		fprintf(file, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
			ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
			ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
			ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz
		);
		fprintf(file, "\n");
	}
	

	fclose(file);

	return;
}


void PrintArray(double* array, size_t n, const char* FileName)
{
	//format date and time. 
	struct tm *ptr;
	time_t lt;
	char str[9];
	lt = time(NULL);
	ptr = localtime(&lt);
	strftime(str, 9, "%m%d%H%M", ptr);

	FILE* file;
	switch (choose)
	{
		case 0:
		{
			char init_file_name[15] = "TestData/AW_";
			///*strcat((char*)(FileName), tmpBuf);
			strcat(init_file_name, FileName);
			strcat(init_file_name, str);
			strcat(init_file_name, ".dat");
			file = fopen(init_file_name, "w");
			break;
		}
		case 1:
		{
			char step_one_file_name[20] = "TestData/DS_";
			strcat(step_one_file_name, FileName);
			strcat(step_one_file_name, str);
			strcat(step_one_file_name, ".dat");
			file = fopen(step_one_file_name, "w");
			break;
		}
	}

	if (!file) perror("cannot open file");
	for (size_t i = 0; i < n; i++)
	{
		fprintf(file, "%-.10lf\n",array[i]);
	}


	fclose(file);


}


void Print_Count_Array(double* ee0_array,size_t * z_arr,size_t * zz_arr,int size,const char* file_name)
{
	FILE* file;
	file = fopen(file_name, "w");
	if (!file) perror("cannot open file");
	for(int i = 0;i<size;i++)
	{
		fprintf(file, "%d\t %.10lf\t %lld \t %lld \n", i, ee0_array[i], z_arr[i], zz_arr[i]);
	}
	fclose(file);
}

void PrintLaserArrays(double* e1_arr, double* e2_arr, double* e_check_arr, size_t size, const char* file_name)
{
	FILE* file;
	file = fopen(file_name, "w");
	if (!file) perror("cannot open file");


	for (int i = 0; i<size; i++)
	{
		double t1 = DX * i * 0.5;
		fprintf(file, "%.10lf \t %.10lf\t %.10lf \t %.10lf \n",(t1/T0_const) , e1_arr[i], e2_arr[i], e_check_arr[i]);
	}
	fclose(file);
}