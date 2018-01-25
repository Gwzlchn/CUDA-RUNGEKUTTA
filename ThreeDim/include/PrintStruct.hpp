#pragma once
#include "nucleus.hpp"
#include <cstdio>
#include <cstring>


//分三个文件保存,用最后一个整形变量选择文件名 
//仅初始化的数据："init_"+ 文件名   0
//第一步数据 "one_" + 文件名		1
//第二步数据 "two_" + 文件名		2

#ifndef PRINTSTRUCT_HPP
#define PRINTSTRUCT_HPP

void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName);

#endif //PRINTSTRUCT_HPP

void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName,int choose)
{
	switch(choose)
	{
	case 0:
		{char init_file_name[6] = "init_";
		strcat(init_file_name, FileName);
		FILE* init = fopen(init_file_name, "w");
		if (!init)
		{
			perror("cannot open file");
		}

		for (auto i = 0; i < n; i++)
		{

			fprintf(init, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
				ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
				ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
				ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
				ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz
			);

			fprintf(init, "\n");
		}
		fclose(init);
		break; }
	case 1:
		{char step_one_file_name[11] = "step_one_";
		strcat(step_one_file_name, FileName);
		FILE* one = fopen(step_one_file_name, "w");
		if (!one)
		{
			perror("cannot open file");
		}

		for (auto i = 0; i < n; i++)
		{

			fprintf(one, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
				ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
				ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
				ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
				ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz
			);

			fprintf(one, "\n");
		}
		fclose(one);
		break; }
	case 2:
		{char step_two_file_name[11] = "step_two_";
		strcat(step_two_file_name, FileName);
		FILE* two = fopen(step_two_file_name, "w");
		if (!two)
		{
			perror("cannot open file");
		}

		for (auto i = 0; i < n; i++)
		{

			fprintf(two, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
				ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
				ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
				ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
				ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz
			);

			fprintf(two, "\n");
		}
		fclose(two);
		break; }


		
	}
	//仅初始化的
	

	



	


	return;
}


