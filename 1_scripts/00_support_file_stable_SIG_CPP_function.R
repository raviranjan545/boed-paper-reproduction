# Numerically stable implementation of SIG utility function

library(Rcpp)

# Create the stable C++ function using basic Rcpp (no Armadillo dependency)
cppFunction('
NumericVector SIGnlmcpp_stable(NumericMatrix Y, NumericMatrix MU1, NumericMatrix MU2, 
                               NumericVector SIG1, NumericVector SIG2) {
    
    int B = Y.nrow();
    int n = Y.ncol();
    
    NumericVector out(B);
    NumericVector log_likelihoods(B);
    
    for (int i = 0; i < B; i++) {
        // Compute log-likelihood under scenario 1 (numerator)
        double ss_resid_1 = 0.0;
        for (int k = 0; k < n; k++) {
            double diff = Y(i, k) - MU1(i, k);
            ss_resid_1 += diff * diff;
        }
        double log_like_1 = ss_resid_1 * (-0.5 / SIG1[i]) + (-0.5 * n * std::log(SIG1[i]));
        
        // Compute log-likelihoods under all scenario 2 models
        for (int j = 0; j < B; j++) {
            double ss_resid_2 = 0.0;
            for (int k = 0; k < n; k++) {
                double diff = Y(i, k) - MU2(j, k);
                ss_resid_2 += diff * diff;
            }
            log_likelihoods[j] = ss_resid_2 * (-0.5 / SIG2[j]) + (-0.5 * n * std::log(SIG2[j]));
        }
        
        // Use log-sum-exp to compute log of average likelihood
        double max_log_like = log_likelihoods[0];
        for (int j = 1; j < B; j++) {
            if (log_likelihoods[j] > max_log_like) {
                max_log_like = log_likelihoods[j];
            }
        }
        
        double sum_exp = 0.0;
        for (int j = 0; j < B; j++) {
            sum_exp += std::exp(log_likelihoods[j] - max_log_like);
        }
        
        double log_avg_likelihood = max_log_like + std::log(sum_exp) - std::log(B);
        
        // Compute SIG utility: log L1 - log(avg L2)
        out[i] = log_like_1 - log_avg_likelihood;
    }
    
    return out;
}')

# Usage in your utilitynlmTemp function:
# Replace the .Call line in the SIG criterion with:
# as.vector(SIGnlmcpp_stable(y, mu1, mu2, 
#                            sam[1:B, colnames(sam) == "sig2"], 
#                            sam[-(1:B), colnames(sam) == "sig2"]))