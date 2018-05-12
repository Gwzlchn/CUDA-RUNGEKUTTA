#ifndef GLOBAL_FUNCS_H
#define GLOBAL_FUNCS_H

#include "nucleus.hpp"

void compute_on_gpu_one(const long pairs,const char* file_name);

void get_min_r_min_p(int nx, int ny, double& min_r, double& min_p);

void NucleiPreRandom(nuclei* Array, const long size);

void NucleiFisrtStep(nuclei* first_array, const long size);

void NucleiSecondStep(nuclei* second_array, nuclei* second_array_fliter, const long size, double* aw, double* ds, unsigned long long* count);

#endif //GLOBAL_FUNCS_H