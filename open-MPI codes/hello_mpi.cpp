#include <iostream>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv){
  int i, myid,size,tag=100;
  char message_send[50], message_recv[50];
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(myid!=0){  // send the msg
    sprintf(message_send, "Hello from process id: %d \n",myid);
    MPI_Send(message_send,50, MPI_CHAR,0,tag,MPI_COMM_WORLD);
  }
  else{
    for(i=1; i<size; i++){
      MPI_Recv(message_recv, 50, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
      cout <<"\n"<<message_recv;
    }
    sprintf(message_send, "Hello from process id: %d \n",myid);
    cout <<"\n"<<message_send;
  }
  MPI_Finalize();
  return 0;
}


