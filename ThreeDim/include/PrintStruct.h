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


