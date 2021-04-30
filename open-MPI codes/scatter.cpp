// MPI_Scatter scatters all the data from root to all the processes

#include <iostream>
#include <string.h>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv){
	int myid, size;
	int *send_buff, recv_buff;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if(myid==0){
		send_buff=(int*)malloc(size*sizeof(int));
		for(int i=0; i<size; i++){
			send_buff[i] = 100 + i*(5+i); 
		}
	}

	MPI_Scatter(send_buff,1,MPI_INT,&recv_buff,1,MPI_INT,0,MPI_COMM_WORLD); // all processes exec this

	if(myid==0){
		cout<<"Received values on processors are:"<<endl;
	}

	cout<<"\t myid: "<<myid<<", recv_buff: "<<recv_buff<<endl;

	if(myid==0){
		free(send_buff);
	}

	MPI_Finalize();
	return 0;
}