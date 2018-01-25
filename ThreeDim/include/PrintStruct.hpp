#pragma once
#include "nucleus.hpp"
#include <cstdio>
#include <cstring>


//分三个文件保存
//仅初始化的数据："init_"+ 文件名
//第一步数据 "one_" + 文件名
//第二步数据 "two_" + 文件名

#ifndef PRINTSTRUCT_HPP
#define PRINTSTRUCT_HPP

void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName);

#endif //PRINTSTRUCT_HPP

void PrintStruct(nuclei* ToSaveNuclei, long long n, const char* FileName)
{
	//仅初始化的
	char init_file_name[6] = "init_";
	strcat(init_file_name, FileName);
	FILE* init=fopen(init_file_name, "w");
	if (!init)
	{
		perror("cannot open file");
	}

	for (auto i = 0; i < n; i++)
	{

		fprintf(init, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].init_first.x, ToSaveNuclei[i].init_first.y, ToSaveNuclei[i].init_first.z,
			ToSaveNuclei[i].init_first.px, ToSaveNuclei[i].init_first.py, ToSaveNuclei[i].init_first.pz,
			ToSaveNuclei[i].init_second.x, ToSaveNuclei[i].init_second.y, ToSaveNuclei[i].init_second.z,
			ToSaveNuclei[i].init_second.px, ToSaveNuclei[i].init_second.py, ToSaveNuclei[i].init_second.pz
			);

		fprintf(init, "\n");
	}
	fclose(init);

	char step_one_file_name[11] = "step_one_";
	strcat(step_one_file_name, FileName);
	FILE* one = fopen(step_one_file_name, "w");
	if (!one)
	{
		perror("cannot open file");
	}

	for (auto i = 0; i < n; i++)
	{

		fprintf(one, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].step_one_first.x, ToSaveNuclei[i].step_one_first.y, ToSaveNuclei[i].step_one_first.z,
			ToSaveNuclei[i].step_one_first.px, ToSaveNuclei[i].step_one_first.py, ToSaveNuclei[i].step_one_first.pz,
			ToSaveNuclei[i].step_one_second.x, ToSaveNuclei[i].step_one_second.y, ToSaveNuclei[i].step_one_second.z,
			ToSaveNuclei[i].step_one_second.px, ToSaveNuclei[i].step_one_second.py, ToSaveNuclei[i].step_one_second.pz
		);

		fprintf(one, "\n");
	}
	fclose(one);



	char step_two_file_name[11] = "step_two_";
	strcat(step_two_file_name, FileName);
	FILE* two = fopen(step_two_file_name, "w");
	if (!two)
	{
		perror("cannot open file");
	}

	for (auto i = 0; i < n; i++)
	{

		fprintf(two, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].step_two_first.x, ToSaveNuclei[i].step_two_first.y, ToSaveNuclei[i].step_two_first.z,
			ToSaveNuclei[i].step_two_first.px, ToSaveNuclei[i].step_two_first.py, ToSaveNuclei[i].step_two_first.pz,
			ToSaveNuclei[i].step_two_second.x, ToSaveNuclei[i].step_two_second.y, ToSaveNuclei[i].step_two_second.z,
			ToSaveNuclei[i].step_two_second.px, ToSaveNuclei[i].step_two_second.py, ToSaveNuclei[i].step_two_second.pz
		);

		fprintf(two, "\n");
	}
	fclose(two);


	return;
}


