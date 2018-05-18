#ifndef __STRUCT_DEFINES_H
#define __STRUCT_DEFINES_H


struct derivative
{
	double px;
	double py;
	double pz;
	double fx;
	double fy;
	double fz;
};


struct nucleus {
	double x;
	double y;
	double z;

	double px;
	double py;
	double pz;
};

struct  nuclei {
	nucleus first;
	nucleus second;
};





#endif //__STRUCT_DEFINES_H
