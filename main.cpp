    //cpp libraries
    #include <iostream>
    #include <vector>
    #include <cmath>
    #include <chrono>

    //my own code
    #include "jacobi.h"
    #include "gauss_seidel.h"
    #include "lu_decomposition.h"
    #include "matrix_vector_initialization.h"
    #include "plotting.h"


    //matplotlibcpp taken from https://github.com/lava/matplotlib-cpp
    #include "matplotlibcpp.h"


    namespace plt = matplotlibcpp;

    #define MAX_ITERATIONS 100000
    #define INITIAL_GUESS 1
    #define a2 -1
    #define a3 -1
    #define MAX_NORM 1e-9
    #define MORE_DEMANDING_NORM 1e-14


    void taskE();

    int main() {
        auto start = std::chrono::high_resolution_clock::now();


        std::vector<double> jacobi_residue_norms;
        std::vector<double> gs_residue_norms;

        const int N = 981;
        double a1 = 10;


        std::vector<std::vector<double>> A = initializeMatrixA(N, a1, a2, a3);
        std::vector<double> b = initializeVectorB(N);

        std::vector<double> x(N, INITIAL_GUESS);
        std::cout << "\================= Task 1" << " =================\n";

        JacobiMethod(A, b, x, MAX_NORM, jacobi_residue_norms, N);
        GSMethod(A, b, x, MAX_NORM, gs_residue_norms, N);
        LUDecomposition(A, b, x);

        plotResidueNorm(jacobi_residue_norms, gs_residue_norms, "Jacobi Method");
        gs_residue_norms.clear();
        jacobi_residue_norms.clear();
        a1 = 3;
        A = initializeMatrixA(N, a1, a2, a3);

        std::cout << "\================= Task 2" << " =================\n";

        JacobiMethod(A, b, x, MAX_NORM, jacobi_residue_norms, N);
        GSMethod(A, b, x, MAX_NORM, gs_residue_norms, N);
        LUDecomposition(A, b, x);

        plotResidueNorm(jacobi_residue_norms, gs_residue_norms, "Jacobi Method");
    
        taskE();


        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "\================= =====" << " =================\n";

        std::cout << "Elapsed time in main function: " << elapsed_seconds.count() << " seconds\n";

        return 0;
    }



    void taskE() {
        double a1 = 10;

        std::vector<int> N_values;
        for (int i = 100; i <= 3100; i += 200) {
            N_values.push_back(i);
        }

        std::vector<double> jacobi_times;
        std::vector<double> gs_times;

        std::vector<double> jacobi_times_2;
        std::vector<double> gs_times_2;

        std::vector<double> lu_times;




        for (int z : N_values) {
            std::cout << "\================= N == "<< z << " =================\n";

            std::vector<std::vector<double>> A = initializeMatrixA(z, a1, a2, a3);
            std::vector<double> b = initializeVectorB(z);
            std::vector<double> x(z, INITIAL_GUESS);

            std::vector<double> jacobi_residue_norms;
            std::vector<double> gs_residue_norms; 

            jacobi_times.push_back(JacobiMethod(A, b, x, MAX_NORM, jacobi_residue_norms, z));
            gs_times.push_back(GSMethod(A, b, x, MAX_NORM, gs_residue_norms, z));

            jacobi_times_2.push_back(JacobiMethod(A, b, x, MORE_DEMANDING_NORM, jacobi_residue_norms, z));
            gs_times_2.push_back(GSMethod(A, b, x, MORE_DEMANDING_NORM, gs_residue_norms, z));

            lu_times.push_back(LUDecomposition(A, b, x));
        }

        std::vector<double> jacobi_times_log;
        for (const auto& norm : jacobi_times) {
            jacobi_times_log.push_back(std::log(norm));
        }

        std::vector<double> gs_times_log;
        for (const auto& norm : gs_times) {
            gs_times_log.push_back(std::log(norm));
        }

        std::vector<double> lu_times_log;
        for (const auto& norm : lu_times) {
            lu_times_log.push_back(std::log(norm));
        }

        std::vector<double> jacobi_times_2_log;
        for (const auto& norm : jacobi_times_2) {
            jacobi_times_2_log.push_back(std::log(norm));
        }

        std::vector<double> gs_times_2_log;
        for (const auto& norm : gs_times_2) {
            gs_times_2_log.push_back(std::log(norm));
        }


        plt::named_plot("Jacobi", N_values, jacobi_times_log);
        plt::named_plot("Gauss", N_values, gs_times_log);

        plt::named_plot("Jacobi for res max=1e-14", N_values, jacobi_times_2_log);
        plt::named_plot("Gauss for res max=1e-14", N_values, gs_times_2_log);

        plt::named_plot("LU", N_values, lu_times_log);
        plt::title(" - Time vs. Matrix Size");
        plt::xlabel("Matrix Size");
        plt::ylabel("Logarithm of Time[ms])");
        plt::legend();
        plt::show();


        plt::named_plot("Jacobi", N_values, jacobi_times);
        plt::named_plot("Gauss", N_values, gs_times);


        plt::named_plot("LU", N_values, lu_times);
        plt::title(" - Time vs. Matrix Size");
        plt::xlabel("Matrix Size");
        plt::ylabel("Time[ms]");
        plt::legend();
        plt::show();
    }

    