#include <iostream>
#include <omp.h>

#include <vector>
#include <stdlib.h>
#include <cmath>

// extern void srand48(long int);
extern double drand48(void);

void usage(const char *name)
{
  std::cout << "usage: " << name
            << " matrix-size nworkers"
            << std::endl;
  exit(-1);
}

////////////////////////////////////////////

std::vector<std::vector<double>> matrix_getrn(int n)
{
  std::vector<std::vector<double>> res(n, std::vector<double>(n));

  // fill the array with random numbers
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

void matrix_LU(std::vector<std::vector<double>> &a, std::vector<int> &pi, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U)
{
  int n = a.size();
  // initialization
  pi = std::vector<int>(n);
  // TODO what is the value outside of this borders ?
  L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
  U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));

  for (int i = 0; i < n; i++)
  {
    pi[i] = i;
    L[i][i] = 1.0;
  }

  // decomposition
  double k, maxN;
  for (int i = 0; i < n; i++)
  {
    maxN = 0, k = 0;
    for (int j = i; j < n; j++)
    {
      if (maxN < std::abs(a[j][i]))
      {
        maxN = std::abs(a[j][i]);
        k = j;
      }
    }

    if (maxN == 0)
      std::cerr << "ERROR\n";

    std::swap(pi[i], pi[k]);
    std::swap(a[i], a[k]);

    std::swap(L[i], L[k]);
    std::swap(L[i][0], L[k][0]);

    U[i][i] = a[i][i];

    for (int j = i + 1; j < n; j++)
    {
      L[j][i] = a[j][i] / U[i][i];
      U[i][j] = a[i][j];
    }
    for (int j = i + 1; j < n; j++)
    {
      for (int z = i + 1; z < n; z++)
        a[j][z] -= L[j][i] * U[i][z];
    }
  }
}

// multiply to matrices
std::vector<std::vector<double>> matrix_mult(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
  int p = A[0].size(), p1 = B.size();
  if (p != p1)
    std::cerr << "ERROR: can not multiply these matrices\n";

  // find result matrix size n x m
  int n = A.size(), m = B[0].size();

  std::vector<std::vector<double>> res(n, std::vector<double>(m, 0));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      for (int k = 0; k < p; k++)
        res[i][j] += A[i][k] * B[k][j];
    }
  }

  return res;
}

// subtraction of matrices
std::vector<std::vector<double>> matrix_sub(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
  int n1 = A.size(), m1 = A[0].size();
  int n2 = B.size(), m2 = B[0].size();

  if (n1 != n2 || m1 != m2)
    std::cerr << "ERROR: cannot subtract these matrices\n";

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
////////////////////////////////////////////

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

#pragma omp parallel
  {
// uncomment the line below to remove the data race
// #pragma omp critical
// std::cout << "hello world from thread "
//           << omp_get_thread_num() << std::endl;
#pragma omp single
    {
      std::vector<std::vector<double>> A = matrix_getrn(matrix_size);

      std::vector<int> pi;
      std::vector<std::vector<double>> L, U;
      matrix_LU(A, pi, L, U); // after this pi, L, U are pointing to the result of LU decomposition

      std::vector<std::vector<double>> P(matrix_size, std::vector<double>(matrix_size, 0));
      for (int i = 0; i < matrix_size; i++)
        P[i][pi[i]] = 1;

      auto PA = matrix_mult(P, A), LU = matrix_mult(L, U);
      auto R = matrix_sub(PA, LU);
      std::cout << "L2,1: " << matrix_L21(R) << '\n';
    }
  }
  return 0;
}
