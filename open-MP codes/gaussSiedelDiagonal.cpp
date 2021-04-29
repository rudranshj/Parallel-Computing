// Author: Rudransh Jaiswal 

#include<iostream>
#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<time.h>
#include<vector>
using namespace std;


int main(int argc, char* argv[]){

	int thread_count=1;
	if(argc==2){
		thread_count=strtol(argv[1], NULL, 10);
	}
	else{
		cout<<"Expected thread_count as cmd line arg. Quitting!"<<endl;
		return -1;
	}

	double d=0.1; //d=dx=dy delta
	double xl=-1.0, xu=1.0, yl=-1.0, yu=1.0;
	int n = 1 + (xu-xl)/d;
	int i,j;

	vector<vector<double> > phi(n, vector<double>(n));
	vector<vector<double> > phi_exact(n, vector<double>(n));
	vector<vector<double> > q(n, vector<double>(n));
	double xc,yc,t1,t2;
	clock_t start, end;
	
	start=clock();

	#pragma omp parallel for num_threads(thread_count) collapse(2) default(shared) private(i,j,xc,yc)
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			xc = xl+d*i;
			yc = yl+d*j;
			phi[i][j]=0.0;
			phi_exact[i][j] = (xc*xc - 1)*(yc*yc - 1);
			q[i][j] = 2*(2 - xc*xc - yc*yc);
		}
	}

	int n_iter=0;
	int i_s,i_e;
	double mse=0.0;
	double mse_check;
	bool do_iter=true;

	#pragma omp parallel default(shared) private(i,j) reduction(+:mse)
	while (do_iter){

		#pragma omp single
		{
		n_iter+=1;
		mse=0.0;
		}

		for(int l=2; l<(2*n - 1); l++){
			i_s=1;
			i_e=l-1;
			if (l>=n){
				i_s=l-n+2;
				i_e=n-2;
			}
			#pragma omp for
			for(i=i_s; i<=i_e; i++){
				j=l-i;
				phi[i][j] = 0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1] + d*d*q[i][j]);
				mse+=(phi[i][j]-phi_exact[i][j])*(phi[i][j]-phi_exact[i][j]);
			}
		}

		#pragma omp barrier

		#pragma omp single
		{
		mse_check=sqrt(mse/(n*n));
		cout<<"mse: "<< mse_check<<endl;
		if(mse_check<0.01){
			do_iter=false;
		}
	}

		#pragma omp barrier


		// cout<<"mse is : "<<mse<<endl;
	
	}

	end=clock();
	double toe = double(end - start) / double(CLOCKS_PER_SEC);
	cout<<"Number of iterations "<<n_iter<<endl;
	cout<<"Time of execution: "<<toe<<endl;

	// for(i=0; i<21; i++){
	// 	cout<<phi[i][15]<<",";
	// }

	return 1;
}