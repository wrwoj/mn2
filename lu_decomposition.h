#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H
#pragma omp parallel


#include <vector>

double LUDecomposition(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x);

#endif // LU_DECOMPOSITION_H
