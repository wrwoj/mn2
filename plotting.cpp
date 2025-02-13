#include "matplotlibcpp.h"
#include "plotting.h"

namespace plt = matplotlibcpp;

void plotResidueNorm(const std::vector<double>& residue_norms1, const std::vector<double>& residue_norms2, const std::string& method_name) {
    // Transforming y-values into log scale
    // (it is necessary 'cause yscale is not included in matplotlibcpp, at least as far as i know)
    std::vector<double> log_residue_norms1;
    for (const auto& norm : residue_norms1) {
        log_residue_norms1.push_back(std::log(norm));
    }
    std::vector<double> log_residue_norms2;
    for (const auto& norm : residue_norms2) {
        log_residue_norms2.push_back(std::log(norm));
    }

    // Plotting residue norms
    plt::named_plot("Jacobi", log_residue_norms1);
    plt::named_plot("Gauss", log_residue_norms2);

    plt::title("Residue Norm vs. Iteration");    
    plt::xlabel("Iteration");
    plt::ylabel("Logarithm   of Residue Norm");
    plt::legend();
    plt::show();
}
