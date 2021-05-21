// Author: Rudransh Jaiswal
// Cholesky Decomposition parallel code

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE		float
#define N		    100
#define SMALLVALUE	0.001
#define n_gangs 25

void zero(TYPE mat[][N]){
	int i,j;
	#pragma acc parallel loop collapse(2) present(mat[0:N][0:N]) num_gangs(n_gangs)	
		for(i=0; i<N; i++){
			for(j=0; j<N; j++)
					mat[i][j]=0.0;
			}
}

void rand_init(TYPE mat[][N]){ //create a random matrix sequentially
	int i,j;
	for (i = 0; i<N; ++i)
		for(j=0; j<N; j++)
			mat[i][j] = ((rand())%13)/1.0 ;	
}

void pd_transf(TYPE mat_pd[][N],TYPE mat[][N]){ //matrix->A*A.T :positive definite
	int i,j,k;
	// #pragma acc parallel loop collapse(2) present(mat_pd[0:N][0:N],mat[0:N][0:N])
	#pragma acc kernel
		for (i = 0; i<N; ++i){
			for(j=0; j<=i; j++){
				for(k=0; k<N; k++)
					mat_pd[i][j] += mat[i][k]*mat[j][k];
			}
		}

	// #pragma acc parallel loop collapse(2) present(mat_pd[0:N][0:N])
	#pragma acc kernel
		for (i = 0; i<N; ++i){
			for(j=0; j<=i; j++){
				mat_pd[j][i]=mat_pd[i][j];
			}
		}
}

void printMat(TYPE a[][N]) {
	int ii,jj;
	for (ii = 0; ii < N; ++ii) {
		printf("[");
		for (jj = 0; jj < N; ++jj)
			printf("%.3f, ", a[ii][jj]);
		printf("],\n");
	}
}

void cholesky(TYPE mat[][N], TYPE lower[][N]) { //reordered loop for better parallelism
	int i,j,k;
	for(j=0; j<N; j++){

		#pragma acc parallel loop present(lower[0:N][0:N])
			for(k=0; k<j; k++)
				lower[j][j] += lower[j][k]*lower[j][k];
		lower[j][j]=sqrt(mat[j][j]-lower[j][j]);

		#pragma acc parallel loop collapse(2) present(lower[0:N][0:N]) num_gangs(n_gangs)
			for(i=j+1; i<N; i++){
				for(k=0; k<j; k++)
					lower[i][j] += lower[i][k]*lower[j][k];
				
				// lower[i][j] = (mat[i][j] - lower[i][j]);
				// if(lower[j][j]>SMALLVALUE)
				// 	lower[i][j]/=lower[j][j];
			}

		#pragma acc parallel loop present(lower[0:N][0:N]) num_gangs(n_gangs)
			for(i=j+1; i<N; i++){				
				lower[i][j] = (mat[i][j] - lower[i][j]);
				if(lower[j][j]>SMALLVALUE)
					lower[i][j]/=lower[j][j];
			}
	}
}
	
int main() {
	TYPE lower[N][N], mat_pd[N][N], mat_aux[N][N];
	clock_t start, end;
	rand_init(mat_aux); //initialise a random matrix serially using rand()
	start=clock();
	#pragma acc data create(mat_pd,lower) copyin(mat_aux)
	{
		zero(mat_pd);
		zero(lower);
		pd_transf(mat_pd,mat_aux); //mat_pd is pos. def.
		printf("Matrix:\n");
		printMat(mat_pd);
		cholesky(mat_pd,lower);
		printf("-------------------\n");
		printf("Lower Triangular Matrix:\n");
		printMat(lower);
	}
	end=clock();
	

	double toe = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Time of Execution: %lf s \n",toe);

	return 0;
}