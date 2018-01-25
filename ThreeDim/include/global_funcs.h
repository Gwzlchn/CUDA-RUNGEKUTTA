#ifndef GLOBAL_FUNS_H
#define GLOBAL_FUNS_H

#include <math.h>

#include "nucleus.hpp"
#include "device_launch_parameters.h"
#include "sci_const.cuh"
#include "Random.h"

void NormalRandom(nuclei* raw_nuclei, const long n);
void InitialNuclei(nuclei* randomed_nuclei, const long raw_count, long & left);

void FirstStep(nuclei* inited_nuclei, const long n);



#endif //GLOBAL_FUNS_H