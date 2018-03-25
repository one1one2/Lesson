#include <iostream>
#include <mpi.h>
#include <cstdlib>

int main(int argc, char *argv[]){
    int m, n, *buffer1, *buffer2;
    MPI_Init(&argc, &argv);
    MPI_Datatype buffer1_type, buffer2_type, medium_type;
    MPI_Status status;

    m = 4; n = 5;
    buffer1 = (int*) calloc(m*n, sizeof(int));
    buffer2 = (int*) calloc(m*n, sizeof(int));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            buffer1[i*n + j] = i*100 + j;
        }
    }
    
    MPI_Type_vector(m, n, -n, MPI_INT, &buffer1_type);
    MPI_Type_commit(&buffer1_type);
    MPI_Type_vector(n, 1, m, MPI_INT, &medium_type);
    MPI_Type_hvector(m, 1, sizeof(int), medium_type, &buffer2_type);
    MPI_Type_commit(&buffer2_type);
    MPI_Sendrecv(buffer1 + (m - 1)*n, 1, buffer1_type, 0, 99, 
                buffer2, 1, buffer2_type, 0, 99, MPI_COMM_SELF, &status);
    MPI_Type_free(&buffer1_type);
    MPI_Type_free(&buffer2_type);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout<<buffer2[i*n + j]<<" ";
        }
    }
    MPI_Finalize();
    return 0;
}
