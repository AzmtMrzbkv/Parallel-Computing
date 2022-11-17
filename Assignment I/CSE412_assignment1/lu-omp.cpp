#include <iostream>
#include <omp.h>

#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <chrono>
// #include <numa.h>

////////////////////////////////////////////

extern void srand48(long int);
extern double drand48(void);

double **matrix_getrn(int n);
void matrix_print(double **A, int n);
void matrix_LU(double **a, int *pi, double **L, double **U, int n);
double **matrix_mult(double **A, double **B, int n);
double **matrix_sub(double **A, double **B, int n);
double matrix_L21(double **A, int n);
void matrix_free(double **A, int n);

////////////////////////////////////////////

void usage(const char *name)
{
  std::cout << "usage: " << name
            << " matrix-size nworkers"
            << std::endl;
  exit(-1);
}

void matrix_free(double **A, int n)
{
  for (int i = 0; i < n; i++)
    delete[] A[i];
  delete[] A;
}

double **matrix_getrn(int n)
{
  double **res = new double *[n];
  for (int i = 0; i < n; i++)
    res[i] = new double[n];

  // fill the array with random numbers in interval [0.0, 1.0]
  srand48(1);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      res[i][j] = drand48();
    }
  }
  return res;
}

void matrix_print(double **A, int n)
{
  std::cout << std::setw(12) << std::setprecision(9) << std::fixed;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n - 1; j++)
      std::cout << A[i][j] << ' ';
    std::cout << A[i][n - 1] << '\n';
  }
}

// multiply to matrices
double **matrix_mult(double **A, double **B, int n)
{
  double **res = new double *[n];
  for (int i = 0; i < n; i++)
    res[i] = new double[n]{};

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        res[i][j] += A[i][k] * B[k][j];

  return res;
}

// subtraction of matrices
double **matrix_sub(double **A, double **B, int n)
{
  double **res = new double *[n];
  for (int i = 0; i < n; i++)
    res[i] = new double[n]{};

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      res[i][j] = A[i][j] - B[i][j];

  return res;
}

// compute L2,1 norm of a given matrix
double matrix_L21(double **A, int n)
{
  double L21 = 0, L21j;
  for (int i = 0; i < n; i++)
  {
    L21j = 0;
    for (int j = 0; j < n; j++)
      L21j += A[j][i] * A[j][i];
    L21 += std::sqrt(L21j);
  }

  return L21;
}

// pass a by value to keep original matrix and fill pi, l, u
void matrix_LU(double **a, int *pi, double **L, double **U, int n)
{
  int j;
  double jmax, jt;

  // initialization
  // #pragma omp parallel for ordered nowait shared(pi, L, n) default(none)
  for (int i = 0; i < n; i++)
  {
    pi[i] = i;
    L[i][i] = 1.0;
  }

  for (int k = 0; k < n; k++)
  {
    jmax = 0;
    j = k;

    for (int i = k; i < n; i++)
    {
      jt = std::fabs(a[i][k]);
      if (jt > jmax)
      {
        jmax = jt;
        j = i;
      }
    }

    if (jmax == 0)
    {
      // std::cout << "ERROR: matrix is singular!\n";
      // TODO: implement exit, mb need to deal with threading
      exit(-1);
    }

#pragma omp parallel sections shared(pi, a, U, k, j) default(none)
    {
#pragma omp section
      std::swap(pi[k], pi[j]);
#pragma omp section
      {
        std::swap(a[k], a[j]);
        U[k][k] = a[k][k];
      }
    }

#pragma omp parallel shared(a, pi, L, U, n, k, j) private(jt) default(none)
    {
#pragma omp for schedule(static)
      for (int z = 0; z < k; z++)
        std::swap(L[k][z], L[j][z]);

#pragma omp for schedule(static)
      for (int i = k + 1; i < n; i++)
        L[i][k] = a[i][k] / U[k][k];

#pragma omp for schedule(static)
      for (int i = k + 1; i < n; i++)
        U[k][i] = a[k][i];

#pragma omp for collapse(2) schedule(static)
      for (int i = k + 1; i < n; i++)
        for (int z = k + 1; z < n; z++)
          a[i][z] = a[i][z] - L[i][k] * U[k][z];
    }
  }
}

int main(int argc, char **argv)
{

  const char *name = argv[0];

  if (argc < 3)
    usage(name);

  int matrix_size = atoi(argv[1]);

  int nworkers = atoi(argv[2]);

  std::cout << name << ": "
            << matrix_size << " " << nworkers
            << std::endl;

  omp_set_num_threads(nworkers);

  // initialize vectors in serial
  double **A = matrix_getrn(matrix_size);
  int *pi = new int[matrix_size];
  double **L = new double *[matrix_size];
  double **U = new double *[matrix_size];
  for (int i = 0; i < matrix_size; i++)
  {
    L[i] = new double[matrix_size]{};
    U[i] = new double[matrix_size]{};
  }

  // create copy of A for L2,1 calculation
  double **AA = new double *[matrix_size];
  for (int i = 0; i < matrix_size; i++)
  {
    AA[i] = new double[matrix_size];
    for (int j = 0; j < matrix_size; j++)
      AA[i][j] = A[i][j];
  }

  // #pragma omp parallel
  // {
  //   // uncomment the line below to remove the data race
  //   #pragma omp critical
  //   std::cout << "hello world from thread "
  //             << omp_get_thread_num() << std::endl;
  // }

  auto start_time = std::chrono::high_resolution_clock::now();
  matrix_LU(AA, pi, L, U, matrix_size);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << std::endl;

  double **P = new double *[matrix_size];
  for (int i = 0; i < matrix_size; i++)
    P[i] = new double[matrix_size]{};

  for (int i = 0; i < matrix_size; i++)
    P[i][pi[i]] = 1;

  auto PA = matrix_mult(P, A, matrix_size), LU = matrix_mult(L, U, matrix_size);
  auto R = matrix_sub(PA, LU, matrix_size);
  std::cout << "L2,1: " << matrix_L21(R, matrix_size) << '\n';

  // deallocate memory
  delete[] pi;
  matrix_free(A, matrix_size);
  matrix_free(AA, matrix_size);
  matrix_free(L, matrix_size);
  matrix_free(U, matrix_size);
  matrix_free(P, matrix_size);
  matrix_free(PA, matrix_size);
  matrix_free(LU, matrix_size);
  matrix_free(R, matrix_size);

  return 0;
}
