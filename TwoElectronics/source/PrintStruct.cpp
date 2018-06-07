#include <cstdio>
#include <cstring>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <sys/types.h>  
#include <sys/stat.h>  

#include "../include/PrintStruct.h"
#include "../include/Sci_Constant.h"
#define BUFLEN 255   

void PrintStruct(particle_pair* ToSaveNuclei, size_t n, const char* FileName)
{	
	CreateDir(folder_name.c_str());
	//format date and time. 
	struct tm *ptr;
	time_t lt;
	char str[9];
	lt = time(NULL);
	ptr = localtime(&lt);
	strftime(str, 9, "%m%d%H%M", ptr);

	FILE* file;
	file = fopen(FileName,"w");
	

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
	CreateDir(folder_name.c_str());
	//format date and time. 
	struct tm *ptr;
	time_t lt;
	char str[9];
	lt = time(NULL);
	ptr = localtime(&lt);
	strftime(str, 9, "%m%d%H%M", ptr);

	FILE* file;
	file = fopen(FileName, "w");
	
	if (!file) perror("cannot open file");
	for (size_t i = 0; i < n; i++)
	{
		fprintf(file, "%-.10lf\n",array[i]);
	}


	fclose(file);


}


void Print_Count_Array(double* ee0_array, unsigned long long* z_arr, unsigned long long* zz_arr, int size, const char* file_name)
{
	CreateDir(folder_name.c_str());
	FILE* file;
	file = fopen(file_name, "w");
	if (!file) perror("cannot open file");
	for(int i = 0;i<size;i++)
	{
		fprintf(file, "%d\t %.10lf\t %lld \t %lld \n", i, ee0_array[i], z_arr[i], zz_arr[i]);
	}
	fclose(file);
}


void Print_Count_Array_Once(double ee0, unsigned long long z, unsigned long long zz, int size, const char* file_name)
{
	CreateDir(folder_name.c_str());
	FILE* file;
	file = fopen(file_name, "w");
	if (!file) perror("cannot open file");
	for (int i = 0; i<size; i++)
	{
		fprintf(file, "%d\t %.10lf\t %lld \t %lld \n", i, ee0, z, zz);
	}
	fclose(file);
}




void PrintLaserArrays(double* e1_arr, double* e2_arr, double* e_check_arr, size_t size, const char* file_name)
{
	CreateDir(folder_name.c_str());
	FILE* file;
	file = fopen(file_name, "w");
	if (!file) perror("cannot open file");


	for (size_t i = 0; i<size; i++)
	{
		double t1 = DX * i * 0.5;
		fprintf(file, "%.10lf \t %.10lf\t %.10lf \t %.10lf \n",(t1/T0_const) , e1_arr[i], e2_arr[i], e_check_arr[i]);
	}
	fclose(file);
}

void CreateDir(const char * dir_name)
{

	int isCreate = mkdir(dir_name, S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
	if (!isCreate)
		printf("create path:%s\n", dir_name);

}
