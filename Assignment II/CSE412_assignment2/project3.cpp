#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

char **alloc2D(int m, int n)
{
    // create additional 2 rows for communication
    char **arr = new char *[m + 2];
    for (int i = 0; i < m + 2; i++)
        arr[i] = new char[n];
    return arr;
}

void free2D(char **arr, int m)
{
    for (int i = 0; i < m + 2; i++)
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
    localGrid = alloc2D(k, m);

    if (myrank == 0)
    {
        // read and send initial grid
        char str[m];
        for (int i = 0; i < m; i++)
        {
            if (i / k == 0)
                scanf("%s", localGrid[i + 1]);
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
        if (myrank != npes - 1)
            for (int i = 0; i < k; i++)
                MPI_Recv(localGrid[i + 1], m, MPI_CHAR, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else
            for (int i = 0; i < m % k; i++)
                MPI_Recv(localGrid[i + 1], m, MPI_CHAR, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // loop game of life
    char **newlocalGrid = alloc2D(k, m);

    for (int t = 0; t < n; t++)
    {
        // send to upper neighbor
        if (myrank != 0)
            MPI_Send(localGrid[1], m, MPI_CHAR, myrank - 1, 0, MPI_COMM_WORLD);
        if (myrank != npes - 1)
            MPI_Recv(localGrid[k + 1], m, MPI_CHAR, myrank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // send to lower neighbor
        if (myrank != npes - 1)
            MPI_Send(localGrid[k], m, MPI_CHAR, myrank + 1, 0, MPI_COMM_WORLD);
        if (myrank != 0)
            MPI_Recv(localGrid[0], m, MPI_CHAR, myrank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 1; i <= k; i++)
        {
            if (myrank == npes - 1 && i > m % k)
                break;
            for (int j = 0; j < m; j++)
            {
                // count neighbors
                int neighbors = 0;
                // left right
                if (j > 0 && localGrid[i][j - 1] == '#')
                    neighbors++;
                if (j < m - 1 && localGrid[i][j + 1] == '#')
                    neighbors++;

                // upper
                if (myrank != 0 || i > 1)
                {
                    if (i > 0 && localGrid[i - 1][j] == '#')
                        neighbors++;
                    if (i > 0 && j > 0 && localGrid[i - 1][j - 1] == '#')
                        neighbors++;
                    if (i > 0 && j < m - 1 && localGrid[i - 1][j + 1] == '#')
                        neighbors++;
                }

                // lower
                if (myrank != npes - 1 || i < m % k)
                {
                    if (i < k - 1 && localGrid[i + 1][j] == '#')
                        neighbors++;
                    if (i < k - 1 && j > 0 && localGrid[i + 1][j - 1] == '#')
                        neighbors++;
                    if (i < k - 1 && j < m - 1 && localGrid[i + 1][j + 1] == '#')
                        neighbors++;
                }

                // identify next state of cell
                if (localGrid[i][j] == '#')
                {
                    if (neighbors < 2 || neighbors > 3)
                        newlocalGrid[i][j] = '.';
                    else
                        newlocalGrid[i][j] = '#';
                }
                else
                {
                    if (neighbors == 3)
                        newlocalGrid[i][j] = '#';
                    else
                        newlocalGrid[i][j] = '.';
                }
            }
        }

        // update localGrid
        for (int i = 1; i <= k; i++)
        {
            if (myrank == npes - 1 && i > m % k)
                break;
            for (int j = 0; j < m; j++)
                localGrid[i][j] = newlocalGrid[i][j];
        }

        // MPI_Barrier(MPI_COMM_WORLD);
    }

    if (myrank == 0)
    {
        // receive final grids parts from other processes
        char **finalGrid = alloc2D(m, m);
        for (int i = 0; i < k; i++)
            for (int j = 0; j < m; j++)
                finalGrid[i][j] = localGrid[i + 1][j];

        for (int i = k; i < m; i++)
            MPI_Recv(finalGrid[i], m, MPI_CHAR, i / k, i % k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // print final grid
        for (int i = 0; i < m; i++)
            printf("%s\n", finalGrid[i]);

        free2D(finalGrid, m);
    }
    else
    {
        // send final grid parts to master
        for (int i = 1; i <= k; i++)
        {
            if (myrank == npes - 1 && i > m % k)
                break;
            MPI_Send(localGrid[i], m, MPI_CHAR, 0, i - 1, MPI_COMM_WORLD);
        }
    }

    free2D(newlocalGrid, k);
    free2D(localGrid, k);

    MPI_Finalize();
}