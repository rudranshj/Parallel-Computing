#include<iostream>
#include<stdio.h>
#include<omp.h>
#include<math.h>
using namespace std;

#define PI 3.14159265358

double func(double x){
  return (1.0 + sin(x));
}

double trapIntgPrl(int n, double a, double b){

	int myRank=omp_get_thread_num();
	int thread_count=omp_get_num_threads();

	double h,resCurr;
	int loc_n;
	double loc_a, loc_b;

	h=(b-a)/n; //same width of trapezoids whatever you do
	loc_n=n/thread_count;
	loc_a= a + myRank*loc_n*h;
	loc_b= loc_a + loc_n*h;
	resCurr=(func(loc_a)+func(loc_b))/2;
	for(int i=1; i<=loc_n-1; i++){
		resCurr+=func(loc_a+i*h);
	}
	resCurr*=h;

	return resCurr;
}


int main(int argc, char* argv[]){
	double a,b,finalRes;
	int n;
	int thread_count=1;

	if(argc==2){
		thread_count=strtol(argv[1], NULL, 10);
	}
	else{
		cout<<"No cmd line arg given. huh. quittind"<<endl;
		return 1;
	}

	n=100;
	a=0.0;
	b=PI;
	finalRes=0;


# pragma omp parallel num_threads(thread_count) reduction(+: finalRes)
    finalRes += trapIntgPrl(n, a, b); /* call to trapezoidal rule will be parallel */
	cout<<"The value of integral is:"<<finalRes<<endl;
	return 1;
}

// g++-10 -o intgParallelv2 -fopenmp numIntgParallel_v2.cpp




