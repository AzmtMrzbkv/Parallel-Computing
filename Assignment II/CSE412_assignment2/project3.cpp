#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    // start MPI
    int npes, myrank;
    char **grid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    int m, n, g;
    scanf("%d %d %d", &m, &n, &g);

    grid = new char *[m];
    for (int i = 0; i < m; i++)
    {
        grid[i] = new char[m];
        scanf("%s", grid[i]);
    }

    MPI_Finalize();
}