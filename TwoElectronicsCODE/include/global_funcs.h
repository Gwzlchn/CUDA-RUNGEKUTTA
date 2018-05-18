#ifndef GLOBAL_FUNCS_H
#define GLOBAL_FUNCS_H

#include "nucleus.hpp"
#include <string>
void compute_on_gpu_one(const long pairs,const char* file_name);

void get_min_r_min_p(int nx, int ny, double& min_r, double& min_p);

void NucleiPreRandom(nuclei* Array, const long size);

void NucleiFisrtStep(nuclei* first_array, const long size);

void NucleiSecondStepPreQQ(double* QQ);
void NucleiSecondStepPreECheck(const double* QQ, const double EE0, double* E_check);


void NucleiSecondStepWholeLaserNoStream(nuclei* first_array, const long size, double* QQ);
void NucleiSecondStepWholeLaser(nuclei* first_array, const long size, double* QQ);

#endif //GLOBAL_FUNCS_H