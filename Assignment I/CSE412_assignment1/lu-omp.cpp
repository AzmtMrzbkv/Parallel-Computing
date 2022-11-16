#include <iostream>
#include <omp.h>

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <time.h>

////////////////////////////////////////////

// extern void srand48(long int);
extern double drand48(void);

std::vector<std::vector<double>> matrix_getrn(int n);
void matrix_print(std::vector<std::vector<double>> &A);
void matrix_LU(std::vector<std::vector<double>> a, std::vector<int> &pi, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U);
std::vector<std::vector<double>> matrix_mult(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
std::vector<std::vector<double>> matrix_sub(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
double matrix_L21(std::vector<std::vector<double>> &A);

////////////////////////////////////////////

void usage(const char *name)
{
  std::cout << "usage: " << name
            << " matrix-size nworkers"
            << std::endl;
  exit(-1);
}

std::vector<std::vector<double>> matrix_getrn(int n)
{
  std::vector<std::vector<double>> res(n, std::vector<double>(n));

  // fill the array with random numbers in interval [0.0, 1.0]
  // srand48(0);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      res[i][j] = drand48();
    }
  }
  return res;
}

void matrix_print(std::vector<std::vector<double>> &A)
{
  int n = A.size(), m = A[0].size();

  std::cout << std::setw(12) << std::setprecision(9) << std::fixed;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m - 1; j++)
      std::cout << A[i][j] << ' ';
    std::cout << A[i][m - 1] << '\n';
  }
}

// multiply to matrices
std::vector<std::vector<double>> matrix_mult(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
  int p = A[0].size(), p1 = B.size();
  if (p != p1)
  {
    std::cerr << "ERROR: cannot multiply these matrices\n";
    exit(-1);
  }

  // find result matrix size n x m
  int n = A.size(), m = B[0].size();

  std::vector<std::vector<double>> res(n, std::vector<double>(m, 0));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      for (int k = 0; k < p; k++)
        res[i][j] += A[i][k] * B[k][j];

  return res;
}

// subtraction of matrices
std::vector<std::vector<double>> matrix_sub(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
  int n1 = A.size(), m1 = A[0].size();
  int n2 = B.size(), m2 = B[0].size();

  if (n1 != n2 || m1 != m2)
  {
    std::cerr << "ERROR: cannot subtract these matrices\n";
    exit(-1);
  }

  std::vector<std::vector<double>> res(n1, std::vector<double>(m1));
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < m1; j++)
      res[i][j] = A[i][j] - B[i][j];

  return res;
}

// compute L2,1 norm of a given matrix
double matrix_L21(std::vector<std::vector<double>> &A)
{
  double L21 = 0, L21j;

  int n = A.size(), m = A[0].size();
  for (int i = 0; i < m; i++)
  {
    L21j = 0;
    for (int j = 0; j < n; j++)
      L21j += A[j][i] * A[j][i];
    L21 += std::sqrt(L21j);
  }

  return L21;
}

// pass a by value to keep original matrix and fill pi, l, u
void matrix_LU(std::vector<std::vector<double>> a, std::vector<int> &pi, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U)
{
  int n = a.size(), j;
  double jmax;

// initialization
#pragma omp parallel shared(pi, L, n) default(none)
  {
#pragma omp for nowait schedule(static)
    for (int i = 0; i < n; i++)
      pi[i] = i;

#pragma omp for nowait schedule(static)
    for (int i = 0; i < n; i++)
      L[i][i] = 1.0;
  }

  for (int k = 0; k < n; k++)
  {
    jmax = 0;
    j = k;

    for (int i = k; i < n; i++)
      if (jmax < std::fabs(a[i][k]))
      {
        jmax = std::fabs(a[i][k]);
        j = i;
      }

    if (jmax == 0)
    {
      // std::cout << "ERROR: matrix is singular!\n";
      // TODO: implement exit, mb need to deal with threading
      exit(-1);
    }

    std::swap(pi[k], pi[j]);
    std::swap(a[k], a[j]);
    U[k][k] = a[k][k];

#pragma omp parallel shared(a, pi, L, U, n, k, j) default(none)
    {
#pragma omp for nowait schedule(static)
      for (int z = 0; z < k; z++)
        std::swap(L[k][z], L[j][z]);

#pragma omp for nowait schedule(static)
      for (int i = k + 1; i < n; i++)
        L[i][k] = a[i][k] / U[k][k];

#pragma omp for schedule(static)
      for (int i = k + 1; i < n; i++)
        U[k][i] = a[k][i];

#pragma omp for collapse(2) schedule(static)
      for (int i = k + 1; i < n; i++)
        for (int z = k + 1; z < n; z++)
          a[i][z] -= L[i][k] * U[k][z];
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
  std::vector<std::vector<double>> A = matrix_getrn(matrix_size);
  std::vector<int> pi(matrix_size);
  std::vector<std::vector<double>> L(matrix_size, std::vector<double>(matrix_size, 0));
  std::vector<std::vector<double>> U(matrix_size, std::vector<double>(matrix_size, 0));

  // #pragma omp parallel
  // {
  //   // uncomment the line below to remove the data race
  //   #pragma omp critical
  //   std::cout << "hello world from thread "
  //             << omp_get_thread_num() << std::endl;
  // }
  auto start_time = clock();
  matrix_LU(A, pi, L, U);
  std::cout << (clock() - start_time) << std::endl;

  std::vector<std::vector<double>> P(matrix_size, std::vector<double>(matrix_size, 0));
  for (int i = 0; i < matrix_size; i++)
    P[i][pi[i]] = 1;

  auto PA = matrix_mult(P, A), LU = matrix_mult(L, U);
  auto R = matrix_sub(PA, LU);
  std::cout << "L2,1: " << matrix_L21(R) << '\n';

  return 0;
}
