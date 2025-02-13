#include "gauss_seidel.h"
#include <chrono>
#include <iostream>

double GSMethod(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x, double max_norm, std::vector<double>& residue_norms, int N) {
    auto start = std::chrono::steady_clock::now();

    int iterations = 0;
    double residual_norm = 0;
    std::vector<double> x_prev(N, 1);

    while (true) {
        for (int i = 0; i < N; ++i) {
            double S = 0;
            for (int j = 0; j < i; ++j) {
                S += A[i][j] * x[j];
            }
            for (int j = i + 1; j < N; ++j) {
                S += A[i][j] * x_prev[j];
            }
            x[i] = (b[i] - S) / A[i][i];
        }

        iterations++;
        std::vector<double> residual(N);
        for (int i = 0; i < N; ++i) {
            double sum = 0;
            for (int j = 0; j < N; ++j) {
                sum += A[i][j] * x[j];
            }
            residual[i] = sum - b[i];
        }
        residual_norm = 0;
        for (double res : residual) {
            residual_norm += res * res;
        }
        residual_norm = std::sqrt(residual_norm);
        residue_norms.push_back(residual_norm);

        if (std::isnan(residual_norm) || residual_norm <= max_norm || iterations >= 10000) {
            break;
        }

        x_prev = x;
    }

    std::cout << "\nGauss-Seidel method:\n";
    std::cout << "\tMax norm: " << max_norm << std::endl;
    std::cout << "\tnumber of iterations: " << iterations << std::endl;
    if (std::isnan(residual_norm)) {
        std::cout << "\tWarning: Residual norm became NaN." << std::endl;
    }
    else {
        std::cout << "\tnorm of residual vector: " << residual_norm << std::endl;
    }
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "\tduration [ms]: " << duration << std::endl;
    return duration;
}
