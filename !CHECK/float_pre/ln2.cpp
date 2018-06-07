#include<cmath>
#include<cstdio>
using namespace std;


double reciprocal(size_t i)
{
    return pow(-1,(i+1))*pow(2,i)/double(i);
}

int main(int argc, char const *argv[])
{
    size_t length = 100;
    double* arr = new double[length];
    
    
    for(size_t i = 1;i!=length + 1;i++){
        arr[i-1] = reciprocal(i);
        //printf("arr %d is %-.15lf",i,arr[i-1])
    }

    double answer = 0.0;

    for(size_t i = 0;i!=length;i++){
        answer +=arr[i];

    }

    printf("cpp answer is %-.15lf",answer);


    return 0;
}
