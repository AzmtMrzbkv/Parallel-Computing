#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

char **alloc2D(int m, int n)
{
    char **arr = new char *[m];
    for (int i = 0; i < m; i++)
        arr[i] = new char[n];
    return arr;
}

void free2D(char **arr, int m)
{
    for (int i = 0; i < m; i++)
        delete[] arr[i];
    delete[] arr;
}

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

    if (myrank == 0)
    {
        scanf("%d %d %d", &m, &n, &g);
        // broadcast the dimensions to other processes
        MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
        // receive the dimensions from master process
        MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // allocate memory for local grid
    int k = m / npes + 1;
    if (myrank == npes - 1)
        // last process gets the remaining rows
        k = m % k;

    localGrid = alloc2D(k, m);

    if (myrank == 0)
    {
        // read and send initial grid
        char str[m];
        for (int i = 0; i < m; i++)
        {
            if (i / k == 0)
                scanf("%s", localGrid[i]);
            else
            {
                scanf("%s", str);
                MPI_Send(str, m, MPI_CHAR, i / k, i % k, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        // receive initial grid
        for (int i = 0; i < k; i++)
            MPI_Recv(localGrid[i], m, MPI_CHAR, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // loop game of life
    char **newlocalGrid = alloc2D(k, m);

    // while (n--)
    // {
    //     // TODO update localGrid and send info to other cells
    //     for (int i = 0; i < m; i++)
    //     {
    //         for (int j = 0; j < k; j++)
    //         {

    //             // count neighbors
    //             int neighbors = 0;
    //             if (i > 0 && localGrid[i - 1][j] == '#')
    //                 neighbors++;
    //             if (i < k - 1 && localGrid[i + 1][j] == '#')
    //                 neighbors++;
    //             if (j > 0 && localGrid[i][j - 1] == '#')
    //                 neighbors++;
    //             if (j < m - 1 && localGrid[i][j + 1] == '#')
    //                 neighbors++;

    //             if (i > 0 && j > 0 && localGrid[i - 1][j - 1] == '#')
    //                 neighbors++;
    //             if (i > 0 && j < m - 1 && localGrid[i - 1][j + 1] == '#')
    //                 neighbors++;
    //             if (i < k - 1 && j > 0 && localGrid[i + 1][j - 1] == '#')
    //                 neighbors++;
    //             if (i < k - 1 && j < m - 1 && localGrid[i + 1][j + 1] == '#')
    //                 neighbors++;

    //             // identify next state of cell
    //             if (neighbors < 2)
    //                 newlocalGrid[i][j] = '.';
    //             else if (neighbors == 3)
    //                 newlocalGrid[i][j] = '#';
    //             else if (neighbors > 3)
    //                 newlocalGrid[i][j] = '.';
    //             else
    //                 newlocalGrid[i][j] = localGrid[i][j];
    //         }

    //         // update localGrid
    //         for (int i = 0; i < k; i++)
    //             for (int j = 0; j < m; j++)
    //                 localGrid[i][j] = newlocalGrid[i][j];
    //     }
    // }

    // master process prints the final state
    if (myrank == 0)
    {
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < m; j++)
                printf("%c", localGrid[i][j]);
            printf("\n");
        }
    }

    free2D(newlocalGrid, k);
    free2D(localGrid, k);

    MPI_Finalize();
}