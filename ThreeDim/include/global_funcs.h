#pragma once
#include "nucleus.hpp"
#ifndef GLOBAL_FUNS_H
#define GLOBAL_FUNS_H

void NormalRandom(nuclei* raw_nuclei, const long n);
void InitialNuclei(nuclei* randomed_nuclei, const long raw_count, long & left);

void FirstStep(nuclei* inited_nuclei, const long n);



#endif //GLOBAL_FUNS_H
