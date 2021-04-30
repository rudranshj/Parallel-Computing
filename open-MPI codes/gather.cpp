// MPI_Gather send all the data to the root 

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv){
	int myid,size,tag=100;
	int send_buff, *recv_buff;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if(myid==0){
		recv_buff=(int*)malloc(size*sizeof(int));
	}

	send_buff = 100 + myid*myid;
	MPI_Gather(&send_buff,1,MPI_INT,recv_buff,1,MPI_INT,0,MPI_COMM_WORLD);

	if(myid==0){
		cout<<"Received values on root are:\n";
		for(int i=0; i<size; i++){
			cout<<"\t"<<recv_buff[i]<<endl;
		}
		free(recv_buff);
	}


	MPI_Finalize();
	return 0;
}