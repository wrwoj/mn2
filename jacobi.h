#ifndef JACOBI_H
#define JACOBI_H

#include <vector>
#include <cmath>

double JacobiMethod(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x, double max_norm, std::vector<double>& residue_norms, int N);

#endif // JACOBI_H
