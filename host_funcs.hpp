#ifndef HOSTFUNCTIONS_HPP
#define HOSTFUNCTIONS_HPP
#include <stdio.h>

 void StoreData(double *Matrix, const int NX,const int NY,const char name[]);



#endif // HOSTFUNCTIONS_HPP

 void StoreData(double *Matrix, const int NX,const int NY,const char name[])
{
   
	
	FILE* fp;
	fp = fopen(name, "w");
    if (!fp)
    {
        perror("cannot open file");
	}
	
    for (int i = 0; i < NX; i++)
    {
		
		for(int j=0;j<NY;j++){
			fprintf(fp,"%-.10lf\t\t",*(Matrix+j*NX+i));
			
		}
		fprintf(fp,"\n");
       
    }

	return;
}

