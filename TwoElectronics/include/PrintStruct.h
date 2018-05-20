//分三个文件保存,用最后一个整形变量选择文件名 
//仅初始化的数据："init_"+ 文件名   0
//第一步数据 "one_" + 文件名		1
//第二步数据 "two_" + 文件名		2

#ifndef PRINTSTRUCT_HPP
#define PRINTSTRUCT_HPP

#include "./Struct_Defines.h"


void PrintStruct(particle_pair* ToSaveNuclei, size_t n, const char* FileName);


//0 aw     1 ds
void PrintArray(double* array, size_t n, const char* FileName);

void Print_Count_Array(double* ee0_array, size_t * z_arr, size_t * zz_arr, int size, const char* file_name);


void PrintLaserArrays(double* e1_arr, double* e2_arr, double* e_check_arr, size_t size, const char* file_name);
#endif //PRINTSTRUCT_HPP


