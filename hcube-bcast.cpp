#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>

void Hcube_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  int rank, np;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &np);

  for (int direction = 1; direction < np; direction*=2) {
    int partner = rank ^ direction; // bit-wise XOR
    if (partner >= np) continue;

    if (rank < direction) {
      MPI_Send(buffer, count, datatype, partner, 0, comm);
    } else if (rank < 2*direction) {
      MPI_Recv(buffer, count, datatype, partner, 0, comm, MPI_STATUS_IGNORE);
    }
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  long N = 1000000;

  int rank;
  MPI_Comm_rank(comm, &rank);
  long* data = (long*) malloc(N*sizeof(long));
  if (rank == 0) { // initialize data
    for (long i = 0; i < N; i++) {
      data[i] = i;
    }
  }

  Hcube_Bcast(data, N, MPI_LONG, 0, comm);

  for (long i = 0; i < N; i++) { // Check
    if (data[i] != i) {
      fprintf(stderr, "Hcube_Bcast failed\n");
      exit(1);
    }
  }

  MPI_Finalize();
  return 0;
}

