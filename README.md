R Code to reproduce the numerical results (including the application) in the manuscript "Frequentist-calibrated Bayesian group sequential design with dynamic borrowing".

**Repository structure and usage**


1. *numerical_assessment.R*: reproduces the fixed-analysis-prior numerical assessment of the manuscript (L = 3 stages, cumulative sample sizes nℓ​ = (100, 200, 300)) — simply source and run the file to reproduce an example where the true effect theta = 0.3, analysis prior N(mu_a = 0.3, σ2/n_a) with n_a = 0 (non-informative prior).

2. *application.R*: reproduces the SSTARLET phase III application (Section 5): a 4-stage design (nℓ​ = c(100, 200, 300, 400)) with dynamic borrowing from historical control- and treatment-arm data. Historical data, design parameters, and both borrowing settings (τ2 = 0 and tau = "upper_bound", i.e. τ2=σ2/nℓ​) are set at the top of the script — simply source and run the file to reproduce all reported results and the ERR_curves.pdf figure (Figure 1).

**Package Requirements**
gsDesign (O'Brien-Fleming group sequential bounds)

**Disclaimer**
The scripts contain the main functions used to implement the simulation studies described in the manuscript. They can be modified to accommodate alternative design parameters (e.g., number of interim analyses, sample sizes, historical data, or borrowing settings).
