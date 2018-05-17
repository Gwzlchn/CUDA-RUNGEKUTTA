#include <cstdio>
#include <cstring>
#include "../include/PrintStruct.h"
#include <cstdio>
#include <cstring>
#include <ctime>
#define BUFLEN 255   

void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName, int choose)
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
	for (long i = 0; i < n; i++)
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


void PrintArray(double* array, long long n, const char* FileName, int choose)
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
	for (long i = 0; i < n; i++)
	{
		fprintf(file, "%-.10lf\n",array[i]);
	}


	fclose(file);


}



void PrintE1_E2(double *e1_e2,long long size,char* file_name)
{
	FILE *outFile;
	outFile = fopen(file_name, "w");
	if (!outFile)
	{
		perror("cannot open the file");
	}
	for(long i = 0;i<size;i=i+2)
	{
		fprintf(outFile,"%d\t%-.10lf\t%-.10lf\n",int(i/2),e1_e2[i],e1_e2[i+1])
	}
	fclose(outFile);
	return;

}

