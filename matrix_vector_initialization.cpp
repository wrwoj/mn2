#include "matrix_vector_initialization.h"
#include <cmath>

std::vector<std::vector<double>> initializeMatrixA(int size, double a1, double a2, double a3) {
    std::vector<std::vector<double>> A(size, std::vector<double>(size, 0));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                A[i][j] = a1;
            }
            else if (std::abs(i - j) == 1) {
                A[i][j] = a2;
            }
            else if (std::abs(i - j) == 2) {
                A[i][j] = a3;
            }
        }
    }
    return A;
}

std::vector<double> initializeVectorB(int size) {
    std::vector<double> b(size, 0);
    for (int i = 0; i < size; ++i) {
        b[i] = std::sin((i) * 3.0);
    }
    return b;
}
