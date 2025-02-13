#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include <vector>
#include <cmath>

double GSMethod(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x, double max_norm, std::vector<double>& residue_norms, int N);

#endif // GAUSS_SEIDEL_H
