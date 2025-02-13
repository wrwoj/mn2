#pragma omp parallel
#include "lu_decomposition.h"
#include <chrono>
#include <iostream>

double LUDecomposition(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x) {
    auto start = std::chrono::steady_clock::now();

    int N = b.size();
    int iterations = 0;
    double residual_norm = 0;


    // Initialize matrices L and U
    std::vector<std::vector<double>> L(N, std::vector<double>(N, 0));
    std::vector<std::vector<double>> U(N, std::vector<double>(N, 0));

    // Perform LU decomposition
    for (int i = 0; i < N; ++i) {
        L[i][i] = 1; // Diagonal elements of L are 1
        for (int j = 0; j < N; ++j) {
            if (i <= j) {
                double sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            }
            if (i > j) {
                double sum = 0;
                for (int k = 0; k < j; ++k) {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = (A[i][j] - sum) / U[j][j];
            }
        }
    }

    // Solve LY = B for Y using forward substitution
    std::vector<double> Y(N, 0);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * Y[j];
        }
        Y[i] = (b[i] - sum) / L[i][i];
    }

    // Solve UX = Y for X using backward substitution
    for (int i = N - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < N; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (Y[i] - sum) / U[i][i];
    }

    // Calculate the residual vector
    std::vector<double> residual(N);
    for (int i = 0; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j] * x[j];
        }
        residual[i] = sum - b[i];
    }

    // Calculate the norm of the residual vector
    residual_norm = 0;
    for (double res : residual) {
        residual_norm += res * res;
    }
    residual_norm = std::sqrt(residual_norm);

    std::cout << "LU decomposition method:\n";
    std::cout << "\tnorm of residual vector: " << residual_norm << std::endl << std::endl;
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "\tduration [ms]: " << duration << std::endl;
    return duration;
}
