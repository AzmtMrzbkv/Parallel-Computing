#include <iostream>
#include <omp.h>

#include <vector>

void usage(const char *name)
{
  std::cout << "usage: " << name
            << " matrix-size nworkers"
            << std::endl;
  exit(-1);
}

////////////////////////////////////////////

std::vector<std::vector<double>> get_matrix(int n)
{
  return {};
}

void LU_decompose(std::vector<std::vector<double>> &a, std::vector<double> &pi, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U)
{
  int n = a.size();
  // initialization
  // TODO what is the value outside of this borders ?
  for (int i = 0; i < n; i++)
  {
    pi[i] = i + 1;

    for (int j = 0; j < n; j++)
    {
      if (i == j)
        L[i][j] = 1.0;

      // below diagonal
      else if (i > j)
        U[i][j] = 0.0;

      // above diagonal (j > i)
      else
        L[i][j] = 0.0;
    }
  }

  // decomposition
  double maxN, k; // TODO what is the initial value of this functions ?
  for (int i = 0; i < n; i++)
  {
    maxN = 0;
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
    // TODO continue
  }
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
    std::cout << "hello world from thread "
              << omp_get_thread_num() << std::endl;
  }
  return 0;
}
