// Gauss-Seidel : Red-Black coloring approach
// Author: Rudransh Jaiswal

#include <iostream>
#include <mpi.h>
#include <math.h>
using namespace std;

double fq(double x, double y){
	return 2*(2 - x*x - y*y);
}
double phi_calc(double x, double y){
	return (pow(x,2)-1)*(pow(y,2)-1);
}

int main(int argc, char** argv){
	int myid,size,n,n_iters;
	double *Phi_local,*Phi,*qvals,d, *Phi_final,*Phi_exact,tol,tot_tol,eps,start_t,end_t;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	n=22; // n should be of the form size*k+2; size=no of processes
	d=2.0/(n-1);
	n_iters=0;
	eps=0.01;
	tot_tol=0.0;

	if(size==1 || (n-2)%size > 0){
		if(myid==0){
			cout<<"Value Error: Ensure that (n-2)%n_procs = 0 and n_procs>1, here n="<<n<<", n_procs="<<size<<endl;
			cout<<"   Please set n and n_procs accordingly"<<endl;
			if(size==1){cout<<"   for n_procs=1 run GS_red_black_serial.cpp"<<endl;	}
			cout<<"   Quitting..."<<endl;
		}
		MPI_Finalize();
		return 0;
	}

	if(myid==0){
		start_t=MPI_Wtime();
	}

	qvals=(double*)malloc((n*(n-2)/size)*n*sizeof(double)); //only local
	Phi_exact=(double*)malloc((n*(n-2)/size)*n*sizeof(double));
	for(int i=0; i<(n-2)/size; i++){
		int i1=i+1+myid*(n-2)/size;
		for(int j=0; j<n; j++){
			qvals[i*n+j] = fq(-1+i1*d,-1+j*d);
			Phi_exact[i*n + j] = phi_calc(-1+i1*d,-1+j*d);
		}
	}

	Phi=(double*)malloc(n*(2 + ((n-2)/size) )*sizeof(double));
	for(int i=0; i<(2 + (n-2)/size); i++){
		for(int j=0; j<n; j++){
			Phi[i*n+j]=0.0;
		}
	}
	Phi_local=(double*)malloc((n*(n-2)/size)*sizeof(double)); //only local, temp. stores the new values
	for(int i=0; i<n*(n-2)/size; i++){
		Phi_local[i]=0.0;
	}
	// MAIN ALGORITHM
	while(n_iters==0 or tot_tol>eps){
		tol=0.0;
		// first red (i+j)%2 > 0 and then black points
		if(n_iters>0){
			MPI_Recv(Phi+n,n*(n-2)/size,MPI_DOUBLE,myid,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self
			if(myid==0){				
				MPI_Recv(Phi+n*(1+(n-2)/size),n,MPI_DOUBLE,myid+1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			else if(myid==size-1){
				MPI_Recv(Phi,n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			else{
				MPI_Recv(Phi,n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self-1
				MPI_Recv(Phi+n*(1+(n-2)/size),n,MPI_DOUBLE,myid+1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self
			}
		}

		for(int i=0; i<(n-2)/size; i++){
			for(int j=1; j<n-1; j++){ //0 for j=0 or j=n-1
				if(((n-2)/size + myid +i +j)%2 > 0){
					Phi_local[i*n+j]= 0.25*(Phi[(i+2)*n+j] + Phi[i*n+j] + Phi[(i+1)*n+j+1] + Phi[(i+1)*n+j-1] + d*d*qvals[i*n+j]);
					tol += pow(Phi_local[i*n+j]-Phi_exact[i*n + j], 2);	
				}
			} 
		}

		MPI_Send(Phi_local,n*(n-2)/size,MPI_DOUBLE,myid,0,MPI_COMM_WORLD); //to same process
		if(myid>0){MPI_Send(Phi_local,n,MPI_DOUBLE,myid-1,2,MPI_COMM_WORLD);} //to previous process
		if(myid<size-1){MPI_Send(Phi_local+n*(-1 + (n-2)/size),n,MPI_DOUBLE,myid+1,1,MPI_COMM_WORLD);} //to next process

		MPI_Recv(Phi+n,n*(n-2)/size,MPI_DOUBLE,myid,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self
		if(myid==0){				
			MPI_Recv(Phi+n*(1+(n-2)/size),n,MPI_DOUBLE,myid+1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else if(myid==size-1){
			MPI_Recv(Phi,n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else{
			MPI_Recv(Phi,n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self-1
			MPI_Recv(Phi+n*(1+(n-2)/size),n,MPI_DOUBLE,myid+1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //from Phi_local of self
		}
		
		// black points now
		for(int i=0; i<(n-2)/size; i++){
			for(int j=1; j<n-1; j++){ //0 for j=0 or j=n-1
				if(((n-2)/size + myid +i +j)%2 == 0){
					Phi_local[i*n+j]= 0.25*(Phi[(i+2)*n+j] + Phi[i*n+j] + Phi[(i+1)*n+j+1] + Phi[(i+1)*n+j-1] + d*d*qvals[i*n+j]);
					tol += pow(Phi_local[i*n+j]-Phi_exact[i*n + j], 2);
				}
			} 
		}

		MPI_Send(Phi_local,n*(n-2)/size,MPI_DOUBLE,myid,0,MPI_COMM_WORLD); //to same process
		if(myid>0){MPI_Send(Phi_local,n,MPI_DOUBLE,myid-1,2,MPI_COMM_WORLD);} //to previous process
		if(myid<size-1){MPI_Send(Phi_local+n*(-1 + (n-2)/size),n,MPI_DOUBLE,myid+1,1,MPI_COMM_WORLD);} //to next process

		MPI_Allreduce(&tol, &tot_tol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		n_iters+=1;
		// if(myid==0){cout<<"n_iters, tot_tol: "<<n_iters<<", "<<tot_tol<<endl;}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0){
		Phi_final=(double*)malloc((n*n)*sizeof(double));
		for(int i=0; i<n; i++){
			Phi_final[i]=0.0;
			Phi_final[n*n-i-1]=0.0;
		}
	}
	MPI_Gather(Phi+n,n*(n-2)/size,MPI_DOUBLE,Phi_final+n,n*(n-2)/size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if(myid==0){
		end_t = MPI_Wtime();
		cout<<"For n = "<<n<<", delta = "<<d<<", n_procs = "<<size<<endl;
		cout<<"  n_iters: "<<n_iters<<endl;
		cout<<"  Time of Execution: "<<end_t-start_t<<"s"<<endl;

		// cout<<"  Time of execution: "<<toe<<endl;

		// cout<<"Value of Phi(0.52381,y): "<<endl;
		// for(int i=0; i<n; i++){
		// 	cout<<Phi_final[i*n +16]<<", ";
		// }
		// cout<<endl;	
	}

	MPI_Finalize();
	return 1;
}
