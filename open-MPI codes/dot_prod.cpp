#include <iostream>
// #include <stdio.h>
#include <string.h>
#include <mpi.h>
using namespace std;

float dot_prod(float a[], float b[], int n){
	float ans=0.0;
	for(int i=0; i<n; i++){
		ans += a[i]*b[i];
	}
	return ans;
}


int main(int argc, char** argv){
	int myid, size;
	float *v1, *v2, total;
	float lsum,gsum;
	int n = 5;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	v1=(float*)malloc(n*sizeof(float));
	v2=(float*)malloc(n*sizeof(float));

	for(int i=0; i<n; i++){
		v1[i] = i+1+myid;
		v2[i] = 2*i + 2 + myid; 
	}
	lsum = dot_prod(v1,v2,n);
	MPI_Allreduce(&lsum,&gsum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	free(v1);
	free(v2);

	if(myid==0){
		cout<<"Inner product: "<<gsum<<endl;
	}
	
	MPI_Finalize();
	return 0;
}