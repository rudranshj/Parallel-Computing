// Numerical Derivative (MPI, n=100)
// Author: Rudransh Jaiswal ME17B063

#include <iostream>
#include <mpi.h>
#include <math.h>
using namespace std;

double fx(double x){
	return sin(5*x);
}
double fdx(double x){
	return 5*cos(5*x);
}

int main(int argc, char** argv){
	int myid,size;
	int n=100;
	double *dfdx_total, *dfdx, *fvals, a, b, h,start_t,end_t;
	a=0.0;
	b=3.0;
	h=(b-a)/n;
	MPI_Status status;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if(n%size > 0){
		if(myid==0){
			cout<<"ValueError: Ensure that n%n_procs = 0, here n="<<n<<", n_procs="<<size<<endl;
			cout<<"   Please set n and n_procs accordingly"<<endl;
			cout<<"   Quitting..."<<endl;
		}
		MPI_Finalize();
		return 0;
	}

	if(myid==0){
		start_t=MPI_Wtime();

		dfdx_total=(double*)malloc(n*sizeof(double));
		dfdx = (double*)malloc((n/size)*sizeof(double));
		fvals = (double*)malloc((n/size + 2)*sizeof(double));

		for(int i=0; i<n/size + 2; i++){ // from i=0 to n/size +2
			fvals[i] = fx(a + (i+ myid*(n/size))*h);
		}
		dfdx[0]=(fvals[1]-fvals[0])/h;
		dfdx[1]=(fvals[2]-fvals[1])/h;
		for(int i=2; i<n/size; i++){
			*(dfdx+i) = (-fvals[i+2] + 8*fvals[i+1] - 8*fvals[i-1] + fvals[i-2])/(12*h);
		}
	}

	else if(myid==size-1){
		dfdx = (double*)malloc((n/size)*sizeof(double));
		fvals = (double*)malloc((n/size + 2)*sizeof(double));

		for(int i=0; i<n/size + 2; i++){ //from i=-2 to i=n/size
			fvals[i] = fx(a + (i-2 + myid*(n/size))*h);
		}
		dfdx[n/size -2]=(fvals[n/size  ] - fvals[n/size -1])/h;
		dfdx[n/size -1]=(fvals[n/size + 1]-fvals[n/size ])/h;
		for(int i=0; i<n/size - 2; i++){
			*(dfdx+i) = (-fvals[i+4] + 8*fvals[i+3] - 8*fvals[i+1] + fvals[i])/(12*h);
		}
	}

	else{
		dfdx = (double*)malloc((n/size)*sizeof(double));
		fvals = (double*)malloc((n/size + 4)*sizeof(double));

		for(int i=0; i<n/size + 4; i++){ //from i=-2 to i=n/size +2
			fvals[i] = fx(a + (i-2 + myid*(n/size))*h);
		}
		for( int i= 0; i < n/size ; i++){
			*(dfdx+i) = (-fvals[i+4] + 8*fvals[i+3] - 8*fvals[i+1] + fvals[i])/(12*h);
		}
	}

	MPI_Gather(dfdx,n/size,MPI_DOUBLE,dfdx_total,n/size,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(myid==0){
		end_t=MPI_Wtime();
		cout<<"For n = "<<n<<", h = "<<h<<", n_procs = "<<size<<endl;
		cout<<"  Time of Execution: "<<end_t-start_t<<"s"<<endl;
	}

	// //print the values
	// if(myid==0){
	// 	cout<<"Numerical value of derivative:"<<endl;
	// 	for(int i=0; i<n; i++){
	// 		cout<<dfdx_total[i]<<", ";
	// 	}
	// 	cout<<endl;
	// }

	MPI_Finalize();
	return 0;
}