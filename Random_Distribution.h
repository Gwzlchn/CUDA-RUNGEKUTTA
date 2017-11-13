//
//  Random_Distribution.h
//  cuda
//
//  Created by 陆子旭 on 2017/11/13.
//  Copyright © 2017年 陆子旭. All rights reserved.
//

#ifndef Random_Distribution_h
#define Random_Distribution_h

#include <stdlib.h>
#include <math.h>
#define PI    3.14159265358979323846

//Equidistribution on [min,max]
double rand_m(double min, double max){
    return min+(max-min)*rand()/(RAND_MAX+1.0);
}

//Normal_Distribution(mean_value: miu; variance: sigma) at x
double normal(double x, double miu,double sigma){
    return 1.0/sqrt(2*PI)/sigma*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}

//Normal_Distribution by Rectangle area
double rand_normal_distribution(double miu,double sigma, double min ,double max){
    double x,y,dScope;
    do{
        x=rand_m(min,max);
        y=normal(x,miu,sigma);
        dScope=rand_m(0.0,normal(miu,miu,sigma));
    }while(dScope>y);
    return x;
}

#endif /* Random_Distribution_h */
