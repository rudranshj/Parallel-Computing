#include <iostream>
#include <stdio.h>
#include <omp.h>
#include<math.h>
using namespace std;

// arctan(x)=Sum( (-1)^r * (x^(2r+1))/(2r+1) ); evaluate sum for x=1 you get value of pi/4

int main(int argc, char* argv[]){
	int thread_count=1;
	double res=0.0;


	if(argc==2){
		thread_count=strtol(argv[1], NULL, 10);
	}
	else{
		cout<<"No cmd line arg given. huh. quittind"<<endl;
		return -1;
	}

# pragma omp parallel for num_threads(thread_count) reduction(+: res)
	for (int i=0; i<1000; i++){
		res+=pow(-1,i)/(i*2 + 1);
	}
	res*=4;
	cout<<"estimation of pi:"<<res<<endl;


	int n=1000;
	double res2=0.0;
	int sign=1;
# pragma omp parallel for num_threads(thread_count) default(none) reduction(+: res2) private(sign) shared(n)//iterator i is by default private to loop
	for (int i=0; i<n; i++){ 
	// you could have written a simple if else statement also
		if(i%2) sign=1;
		else sign=-1;

		res2+=sign/(i*2 + 1);
	}
	res2*=4;
	cout<<"estimation of pi (using default directive):"<<res<<endl;

	return 1;
}
