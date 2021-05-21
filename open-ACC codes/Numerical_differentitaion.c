// Author: Rudransh Jaiswal
// Pade Scheme LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// #include <openacc.h> //for compiling using gcc -fopenacc

#define N 1000
#define EPSILON 0.002
#define n_gangs 25

void printArr(double arr[N+1]){
	int i;
	for(i=0; i<=N; i++){
		printf("%lf,",arr[i]);
	}
	printf("\n");
}

void printMat(double a[][N+1]) {
	int ii,jj;
	for (ii = 0; ii <= N; ++ii) {
		printf("[");
		for (jj = 0; jj <= N; ++jj)
			printf("%.3f, ", a[ii][jj]);
		printf("],\n");
	}
}

double func(double x){
	return sin(5*x);
}

void init_fvals(double fx[N+1]){
	int i;
	double a=0.0;
	double h=3.0/N;
	// # pragma acc kernel
	#pragma acc parallel loop present(fx[0:N+1]) num_gangs(n_gangs)
	for(i=0; i<=N; i++){
		fx[i]=func(a+i*h);
	}
}

void init_Yarr(double Y[N+1], double fx[N+1]){
	int i;
	double h=(3.0)/N;
	Y[0]=(-2.5*fx[0] + 2*fx[1] + 0.5*fx[2])/h;
	Y[N]=(2.5*fx[N] - 2*fx[N-1] - 0.5*fx[N-2])/h;

	// # pragma acc kernel
	#pragma acc parallel loop present(Y[0:N+1]) num_gangs(n_gangs)
	for(i=1; i<N; i++){
		Y[i]=3*(fx[i+1]-fx[i-1])/h;
	}
}


// void init_A(double A[][N+1]){
// 	int i,j;
// 	#pragma acc parallel loop collapse(2) present(A[0:N+1][0:N+1])
// 	for(i=0; i<=N; i++){
// 		for(j=0; j<i-1; j++){
// 			A[i][j]=0.0;
// 		}
// 	}
// 	#pragma acc parallel loop collapse(2) present(A[0:N+1][0:N+1])
// 	for(i=0; i<=N; i++){
// 		for(j=i+2; j<=N; j++){
// 			A[i][j]=0.0;
// 		}
// 	}
// 	A[0][0]=1.0;
// 	A[0][1]=2.0;
// 	#pragma acc parallel loop present(A[0:N+1][0:N+1])
// 	for(int i=1; i<N; i++){
// 		A[i][i-1]=1.0;
// 		A[i][i]=4.0;
// 		A[i][i+1]=1.0;
// 	}
// 	A[N][N-1]=2.0;
// 	A[N][N]=1.0;
// }

// void lu_decomp(double A[][N+1]){
// 	int i,j,k;
// 	for(i=0; i<=N; i++){
// 		#pragma acc parallel loop present(A[0:N+1][0:N+1])
// 		for(j=0; j<i; j++){
// 			for(k=0; k<j; k++){
// 				A[i][j] -= A[i][k]*A[k][j];
// 			}
// 			A[i][j] /= (A[j][j] > EPSILON ? A[j][j] : 1);
// 		}
// 		#pragma acc parallel loop present(A[0:N+1][0:N+1]) collapse(2)
// 		for(j=0 ;j<=N; j++){
// 			for(k=0; k<i; k++){
// 				A[i][j]-=A[i][k]*A[k][j];
// 			}
// 		}
// 	}
// }

void init_Aarr(double Aa[N+1],double Ab[N+1],double Ac[N+1]){
	int i;
	Aa[0]=0.0;
	Ab[0]=1.0;
	Ac[0]=2.0;
	Aa[N]=2.0;
	Ab[N]=1.0;
	Ac[N]=0.0;
	#pragma acc parallel loop present(Aa[0:N+1],Ab[0:N+1],Ac[0:N+1]) num_gangs(n_gangs)
	for(i=1; i<N; i++){
		Aa[i]=1.0;
		Ab[i]=4.0;
		Ac[i]=1.0;
	}
}

void init_Aux(double l[N+1],double u[N+1],double Y_[N+1],double Aa[N+1],double Ab[N+1],double Ac[N+1], double Y[N+1]){
	int i;
	u[0]=Ab[0];

	for(i=0; i<N; i++){
		u[i+1]=Ab[i+1] - Ac[i]*Aa[i+1]/u[i];
	}

	#pragma acc parallel loop present(l[0:N+1]) num_gangs(n_gangs)
	for(i=0; i<N; i++){
		l[i+1]=Aa[i+1]/u[i];
	}

	Y_[0]=Y[0];
	for(i=1; i<=N; i++){
		Y_[i] = Y[i] - l[i]*Y_[i-1];
	}

}

// O(n) method and faster than dense parallel method
void solve_lud(double dfdx[N+1],double Y_[N+1],double Ac[N+1], double u[N+1]){
	int i;
	dfdx[N]=Y_[N]/u[N];
	for(i=N-1; i>=0; i--){
		dfdx[i]= (Y_[i] - Ac[i]*dfdx[i+1])/u[i];
	}
}

int main(int argc, char** argv){
	double Aa[N+1],Ab[N+1],Ac[N+1],Y[N+1],fx[N+1],dfdx[N+1];
	double l[N+1],u[N+1],Y_[N+1];
	clock_t start, end;
	start=clock();
	#pragma acc data create(Aa,Ab,Ac,Y,l,u,Y_,fx)
	{
		init_fvals(fx);
		init_Yarr(Y, fx);
		init_Aarr(Aa,Ab,Ac);
		init_Aux(l, u, Y_, Aa, Ab, Ac, Y);
		solve_lud(dfdx,Y_,Ac,u);
		printArr(dfdx);
	}
	end=clock();
	double toe = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Time of Execution: %lf s \n",toe);


	return 0;
}
