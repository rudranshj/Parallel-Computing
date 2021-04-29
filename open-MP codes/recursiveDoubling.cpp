// Author: Rudransh Jaiswal
 Recursive doubling Parallel Program

#include<iostream>
#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<time.h>
#include<vector>
using namespace std;

// pass noOfThreads while executing this program
int main(int argc, char* argv[]){
	int thread_count=1;
	clock_t start, end;
	if(argc==2){
		thread_count=strtol(argv[1], NULL, 10);
	}
	else{
		cout<<"Expected thread_count as cmd line arg. Quitting!"<<endl;
		return 0;
	}
	
	start=clock();

	double a=0.0,b=3.0;
	int n=100,i,j;
	double h=(b-a)/n;
	n++;; // for including the point x=3.0 

	vector<double> fx(n);
	vector<double> dfdx(n);
	vector<double> dfdx_exact(n);
	vector<double> Aa(n);
	vector<double> Ab(n);
	vector<double> Ac(n);
	vector<double> Y(n);
	vector<double> alpha(n);
	vector<double> beta(n);

	#pragma omp parallel for default(shared) private(i)
	for(i=0; i<n; i++){
		fx[i]=sin(5*(a+i*h));
		Aa[i]=1.0;
		Ab[i]=4.0;
		Ac[i]=1.0;
	}

	#pragma omp parallel for default(shared) private(i)
	for(i=1; i<n-1; i++){
		Y[i]=3*(fx[i+1]-fx[i-1])/h;
	}

	Aa[0]=0.0;
	Aa[n-1]=2.0;
	Ab[0]=1.0;
	Ab[n-1]=1.0;
	Ac[0]=2.0;
	Ac[n-1]=0.0;
	Y[0] = (-2.5*fx[0] + 2*fx[1] + 0.5*fx[2])/h;
	Y[n-1] = (2.5*fx[n-1] - 2*fx[n-2] - 0.5*fx[n-3])/h;

	// Elimination Step
	int k_total = ceil(log2(n));
	int k;

	for(k=0; k<k_total; k++){
		int z=int(pow(2,k-1));

		#pragma omp parallel default(shared) private(i)
		{
		#pragma omp for
		for(i=z; i<n; i++){
			alpha[i] = -Aa[i]/Ab[i-z];
		}
		#pragma omp for
		for(i=0; i<n-z; i++){
			beta[i] = -Ac[i]/Ab[i+z];
		}
		#pragma omp for
		for(i=0; i<n; i++){
			Ab[i] += alpha[i]*Ac[i-z] + beta[i]*Aa[i+z];
			Y[i] += alpha[i]*Y[i-z] + beta[i]*Y[i+z];
		}

		#pragma omp single
		for(i=n-1; i>=2*z; i--){
			Aa[i] = alpha[i]*Aa[i-z];
		}	

		#pragma omp single
		for(i=0; i<(n-2*z); i++){
			Ac[i] = beta[i]*Ac[i+z];
		}

		#pragma omp for nowait
		for(i=0; i<2*z; i++){
			Aa[i]=0.0;
		}
		#pragma omp for nowait
		for(i=n-2*z; i<n; i++){
			Ac[i]=0.0;
		}	
		}
	}

	//Substitution step
	#pragma omp parallel for default(shared) private(i)
	for(i=0; i<n; i++){
		dfdx[i]=Y[i]/Ab[i];
		dfdx_exact[i]=5*cos(5*(a+i*h));
	}


	// for(i=0; i<n; i++){
	// 	// cout<<dfdx[i]<<",";
	// 	cout<<dfdx[i]<<" , "<<dfdx_exact[i]<<" :: "<<abs(dfdx_exact[i]-dfdx[i])<<endl;
	// }
	// cout<<endl<<endl;
	// for(i=0; i<n; i++){
	// 	cout<<dfdx[i]<<",";
	// }
	// cout<<endl<<endl;
	// for(i=0; i<n; i++){
	// 	cout<<dfdx_exact[i]<<",";
	// }


	end=clock();
	double toe = double(end - start) / double(CLOCKS_PER_SEC); 
	cout<<"Execution time with n="<<(n-1)<<", threads="<<thread_count<<" is: "<<toe<<"s"<<endl;
	return 1;
}