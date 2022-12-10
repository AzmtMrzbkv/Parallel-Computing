#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int npes, myrank;
    char **localGrid; // 1-D decomposition

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    // read the dimensions of the grids
    int m, n, g;
    scanf("%d %d %d", &m, &n, &g);

    // allocate localGrid and every process reads its local grid
    int k = m / npes + 1;
    localGrid = new char *[k];
    for (int i = 0; i < m; i++)
        if (i / npes == myrank)
        {
            localGrid[i] = new char[m];
            scanf("%s", localGrid[i]);
        }

    // loop game of life
    while (N--)
    {
        // TODO update localGrid and send info to other cells
    }

    // deallocate localGrid
    for (int i = 0; i < localGridSize; i++)
        delete[] grid[i];
    delete[] grid;

    MPI_Finalize();
}