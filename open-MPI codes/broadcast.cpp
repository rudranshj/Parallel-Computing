#include <iostream>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv){
	int i,size,myid,tag=100;
	int buf=0;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if(myid==0){
		buf=108;
	}
	MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myid==0){
		cout<<"\nBroadcasted values on processor are:"<<endl;
	}
	cout<<"\t myid: "<<myid<<", buf: "<<buf<<endl;
	MPI_Finalize();
	return 0;
}