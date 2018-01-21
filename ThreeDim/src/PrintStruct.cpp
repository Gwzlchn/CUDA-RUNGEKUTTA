#include "../include/nucleus.hpp"
#include <cstdio>

void PrintStruct(nuclei* ToSaveNuclei,long long n,const char* FileName)
{
	FILE* fp;
	fp = fopen(FileName, "w");
	if (!fp)
	{
		perror("cannot open file");
	}

	for (auto i = 0; i < n; i++)
	{

		fprintf(fp, "%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t%-.10lf\t",
			ToSaveNuclei[i].first.x, ToSaveNuclei[i].first.y, ToSaveNuclei[i].first.z,
			ToSaveNuclei[i].first.px, ToSaveNuclei[i].first.py, ToSaveNuclei[i].first.pz,
			ToSaveNuclei[i].second.x, ToSaveNuclei[i].second.y, ToSaveNuclei[i].second.z,
			ToSaveNuclei[i].second.px, ToSaveNuclei[i].second.py, ToSaveNuclei[i].second.pz);
		fprintf(fp, "\n");
	}

	return;
}
