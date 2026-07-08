
# Reproduce results of Section 4 (main manuscript) and of Web Appendix C and D
# of the Supplementary Material

# -----------------------------------------------------------------------------
# Function: dh
# -----------------------------------------------------------------------------
# Computes the dissimilarity between the current study data and a historical
# dataset using a Hellinger-based distance between two Gaussian
# distributions with equal variance but different sample sizes.
#
# Arguments:
#   x_c  : observed mean in the current sample at the current stage
#   x_h  : observed mean in the historical dataset
#   s_c  : standard deviation of the current sample ( = sqrt(sigma2))
#   s_h  : standard deviation of the historical dataset ( = sqrt(sigma2))
#   n_c  : sample size of the current study at the current stage
#   n_h  : sample size of the historical dataset
#
# Returns:
#   A scalar in [0, 1], where 0 indicates perfect similarity and 1 indicates
#   maximum dissimilarity.
# -----------------------------------------------------------------------------
dh = function(x_c, x_h, s_c, s_h, n_c, n_h){
  s_c = s_c/sqrt(n_c)
  s_h = s_h/sqrt(n_h)
  d = 1 - sqrt( (2*s_c*s_h)/ (s_c^2 + s_h^2))*
    exp(-1/4*(x_c-x_h)^2/(s_c^2 + s_h^2))
  
  return(sqrt(d))
}


# -----------------------------------------------------------------------------
# Function: PO
# -----------------------------------------------------------------------------
# Computes the Posterior Odds (PO) in favor of H_1: theta > theta_0 given the
# current observed test statistic, the current sample size, and the parameters
# of the analysis prior.
#
# Arguments:
#   t    : standardized test statistic at the current stage,
#          t = sqrt(n / sigma2) * (x_bar - theta_0)
#   n    : cumulative sample size at the current stage
#   mu_a : mean of the analysis prior
#   n_a  : effective sample size (ESS) of the analysis prior;
#          n_a = 0 corresponds to a non-informative prior
#
# Returns:
#   A positive scalar: the posterior odds in favor of H_1.
# -----------------------------------------------------------------------------
PO = function(t, n, mu_a, n_a){
  
  t_a = sqrt(n_a)*(mu_a - theta_0)/sqrt(sigma2)
  
  t_post = (sqrt(n_a)*t_a + sqrt(n)*t) / sqrt(n_a + n)
  
  omega_post = pnorm(t_post)/(1-pnorm(t_post))
  
  omega_post
}


# -----------------------------------------------------------------------------
# Function: PO_inv
# -----------------------------------------------------------------------------
# Computes the inverse of PO: given a target threshold on the Posterior Odds
# scale (kb), returns the corresponding frequentist threshold on the
# standardized test statistic scale (k_f).
#
# Arguments:
#   kb   : target threshold on the Posterior Odds scale
#          (e.g., the Bayesian-calibrated bound k_b)
#   n    : cumulative sample size at the current stage
#   mu_a : mean of the analysis prior
#   n_a  : effective sample size (ESS) of the analysis prior
#
# Returns:
#   k_f: the standardized test statistic value such that PO(k_f, n, mu_a, n_a) = kb.
# -----------------------------------------------------------------------------
PO_inv = function(kb, n, mu_a, n_a){
  
  t_a = sqrt(n_a)*(mu_a - theta_0)/sqrt(sigma2)
  
  kf = sqrt((n_a + n)/n)*qnorm(kb/(kb+1)) - sqrt(n_a/n)*t_a
  
  kf
}


# -----------------------------------------------------------------------------
# Function: update_analysis_prior
# -----------------------------------------------------------------------------
# Updates the analysis prior parameters (mu_a, n_a) at stage l >= 2 by
# borrowing information from historical datasets. Borrowing is weighted by
# the similarity between the current data and each historical dataset (via dh),
# and may be penalized by a heterogeneity parameter tau.
#
# If no historical data are provided, the prior parameters are returned
# unchanged.
#
# Arguments:
#   l          : current stage index (integer >= 1)
#   mu_a       : prior mean at the previous stage
#   n_a        : prior ESS at the previous stage
#   current    : vector c(x_l, n_l): observed mean and cumulative sample size
#                at stage l-1 (used to compute similarity with historical data)
#   historical : list of vectors c(x_h, n_h), one per historical dataset;
#                if NULL, no borrowing is performed
#   tau        : heterogeneity penalization parameter:
#                  - numeric value >= 0: fixed penalization
#                    (tau = 0 means no extra penalization beyond similarity)
#                  - "others" (non-numeric): tau is estimated adaptively
#                    as sigma2 / n_l (variance of the current MLE)
#
# Returns:
#   A named list with:
#     mu_a   : updated prior mean (weighted average of historical means)
#     n_a    : updated prior ESS
#     w      : vector of weights assigned to each historical dataset
#     s      : vector of similarity scores (1 - dh) for each historical dataset
#     ess    : vector of effective sample size contributions per historical dataset
#     tau_sq : the tau value used in the current update
# -----------------------------------------------------------------------------
update_analysis_prior = function(l, mu_a, n_a, current = NULL, historical = NULL, tau = 0){
  if(is.null(historical)){
    list(mu_a = as.numeric(mu_a), n_a = as.numeric(n_a), w = NULL, s = NULL, ess = NULL, tau_sq = NA)
  } else{ 
    x_l = current[1]
    n_l = current[2]
    
    x_h  <- sapply(historical, `[`, 1)   
    n_h <- sapply(historical, `[`, 2) 
    H = length(historical)
    
    if(is.numeric(tau)){
      simil = apply(rbind(x_h,n_h), 2, function(x){
        
        s = (1- round(dh(x_l, x[1], sqrt(sigma2), sqrt(sigma2), n_l,x[2]),3))
        w = s*x[2]/(sigma2 + s*x[2]*tau)
        ess = s*x[2]
        
        rbind(w,s,ess)
      })
      w = simil[1,]
      s = simil[2,]
      ess = simil[3,]
      
      mu_a = sum(w*x_h)/sum(w)
      n_a = sigma2/(1/sum(w) + tau)
    } else{
      # Adaptive tau: estimated as the variance of the current MLE
      tau_sq = sigma2/n_l
      
      simil = apply(rbind(x_h,n_h), 2, function(x){
        
        s = (1- round(dh(x_l, x[1], sqrt(sigma2), sqrt(sigma2), n_l,x[2]),3))
        w = s*x[2]/(sigma2 + s*x[2]*tau_sq)
        ess = s*x[2]
        
        rbind(w,s,ess)
      })
      w = simil[1,]
      s = simil[2,]
      ess = simil[3,]
      
      mu_a = sum(w*x_h)/sum(w)
      n_a = sigma2/(1/sum(w) + tau_sq)
    }
    tau_sq = tau
    list(mu_a = as.numeric(mu_a), n_a = as.numeric(n_a), w = w, s = s, ess = ess, tau_sq = as.numeric(tau_sq))
  }
}


# -----------------------------------------------------------------------------
# Function: PO_sequential
# -----------------------------------------------------------------------------
# Computes the Posterior Odds sequentially across all planned interim stages.
# Can operate in two modes:
#
#   (1) Trial mode (compute_threshold = FALSE):
#       Given a vector of observed MLEs (mle_l), computes the PO at each stage.
#       Also returns Bayesian decision thresholds (k_b), frequentist thresholds
#       (k_f) on the MLE scale (c_b), and the updated prior parameters.
#
#   (2) Threshold mode (compute_threshold = TRUE):
#       Given a vector of frequentist bounds (k_l), computes the corresponding
#       PO thresholds without observed data, returning the Bayesian-calibrated
#       bounds k_b.
#
# In both modes, the analysis prior is optionally updated at each stage via
# dynamic historical borrowing (see update_analysis_prior).
#
# Arguments:
#   mle_l             : vector of cumulative MLEs at each stage; required when
#                       compute_threshold = FALSE
#   k_l               : vector of frequentist decision bounds (e.g., O'Brien-
#                       Fleming bounds) at each stage; required when
#                       compute_threshold = TRUE, and optionally used alongside
#                       mle_l to compute k_b and k_f
#   n_l               : vector of cumulative sample sizes at each planned stage
#   mu_a              : initial prior mean (at stage 1)
#   n_a               : initial prior ESS (at stage 1); n_a = 0 means
#                       non-informative prior
#   compute_threshold : logical; if FALSE (default), computes PO from observed
#                       MLEs; if TRUE, computes PO thresholds from k_l bounds
#   historical        : list of historical datasets, each c(x_h, n_h); passed
#                       to update_analysis_prior; NULL means no borrowing
#   tau               : heterogeneity penalization; passed to
#                       update_analysis_prior (see that function for details)
#
# Returns (compute_threshold = FALSE):
#   A named list with:
#     po    : vector of PO values at each stage
#     mu_a  : vector of prior means used at each stage (after borrowing updates)
#     n_a   : vector of prior ESS values used at each stage
#     w     : matrix (H x L) of borrowing weights per historical source and stage
#     s     : matrix (H x L) of similarity scores per historical source and stage
#     ess   : matrix (H x L) of ESS contributions per historical source and stage
#     tau_sq: matrix (1 x L) of tau values used at each stage
#     k_b   : vector of Bayesian PO thresholds (derived from k_l via non-
#             informative prior), one per stage; NULL if k_l not supplied
#     k_f   : vector of Bayesian PO thresholds (derived from k_l via analysis
#             prior), one per stage; NULL if k_l not supplied
#     c_b   : vector of frequentist MLE thresholds corresponding to k_b,
#             translated via PO_inv; NULL if k_l not supplied
#
# Returns (compute_threshold = TRUE):
#   A named vector k_b: the PO value at each stage when evaluated at the
#   frequentist bound k_l, using the (possibly updated) analysis prior.
# -----------------------------------------------------------------------------
PO_sequential <- function(mle_l = NULL, k_l = NULL, n_l, mu_a, n_a, compute_threshold = F, historical = NULL, tau = 0){
  
  weight = matrix(nrow = length(historical), ncol = length(n_l))
  s = matrix(nrow = length(historical), ncol = length(n_l))
  ess = matrix(nrow = length(historical), ncol = length(n_l))
  tau_sq = matrix(nrow = 1, ncol = length(n_l))
  
  po = c()
  l = 1
  mu_a_l = c(mu_a)
  n_a_l = c(n_a)
  if(is.null(k_l)){
    k_b = NULL
    k_f = NULL
    c_b = NULL
  }else{
    k_b = c(PO(k_l[1], n_l[1], mu_a_l[1], n_a = 0))
    k_f = c(PO(k_l[1], n_l[1], mu_a_l[1], n_a_l[1]))
    c_b = c(PO_inv(k_b[1], n_l[1], mu_a_l[1], n_a_l[1]))
  }
  
  
  if(compute_threshold == F){ 
    t_l = c()
    
    # Stage l = 1
    t_l[1] = sqrt(n_l[1]/sigma2)*(mle_l[1]-theta_0)
    po[1] <- PO(t_l[1], n_l[1], mu_a_l[1], n_a_l[1])
    
    
    # Stage l >= 2
    for (l in 2:length(n_l)){
      
      current = c(mle_l[l-1], n_l[l-1])
      # Update the analysis prior using borrowing from historical data
      analysis_prior = update_analysis_prior(l, 
                                             mu_a = mu_a_l[l-1],
                                             n_a = n_a_l[l-1],
                                             current = current, historical = historical, tau = tau)
      mu_a_l[l] = analysis_prior$mu_a
      n_a_l[l] = analysis_prior$n_a
      weight[,l] = analysis_prior$w
      s[,l] = analysis_prior$s
      ess[,l] = analysis_prior$ess
      tau_sq[,l] = analysis_prior$tau_sq
      
      t_l[l] = sqrt(n_l[l]/sigma2)*(mle_l[l]-theta_0)
      
      po[l] <- PO(t_l[l], n_l[l], mu_a_l[l], n_a_l[l])
      
      if(!is.null(k_l)){
        k_b = c(k_b, PO(k_l[l], n_l[l], mu_a_l[l], n_a = 0))
        k_f = c(k_f, PO(k_l[l], n_l[l], mu_a_l[l], n_a_l[l]))
        c_b = c(c_b, PO_inv(k_b[l], n_l[l], mu_a_l[l], n_a_l[l]))
      }
    }
    return(list(po = po, mu_a = mu_a_l, n_a = n_a_l, w = weight, s = s, ess = ess, tau_sq = tau_sq, k_b = k_b, k_f = k_f, c_b = c_b))
  }
  
  if(compute_threshold == T){ # requires k_l
    
    # Stage l = 1
    po[1] <- PO(k_l[1], n_l[1], mu_a_l[1], n_a_l[1])
    
    # Stage l >= 2
    for (l in 2:length(n_l)){
      
      analysis_prior = update_analysis_prior(l, 
                                             mu_a = mu_a_l[l-1],
                                             n_a = n_a_l[l-1],
                                             current = NULL, historical = historical, tau = tau)
      mu_a_l[l] = analysis_prior$mu_a
      n_a_l[l] = analysis_prior$n_a
      
      po[l] <- PO(k_l[l], n_l[l], mu_a_l[l], n_a_l[l])
    }
    return(rbind(k_b = po))
  }
  
}



# =============================================================================
# Numerical Assessment
# =============================================================================

require(gsDesign)

# -----------------------------------------------------------------------------
# Global parameters
# sigma2   : known population variance (assumed equal across all studies)
# theta_0  : null hypothesis value for the treatment effect
# n_l      : vector of cumulative sample sizes at each planned interim stage
# -----------------------------------------------------------------------------
sigma2 = 4
theta_0 = 0
n_l = c(100, 200, 300)

# -----------------------------------------------------------------------------
# O'Brien-Fleming group sequential bounds (one-sided, alpha = 0.05, 80% power)
# These serve as the frequentist reference thresholds k_l at each stage.
# -----------------------------------------------------------------------------
obf <- gsDesign(k = 3,          # number of interim stages
                test.type = 1,  # one-sided test
                alpha = 0.05, 
                beta = 0.2,
                sfu = sfLDOF)$upper$bound 

# -----------------------------------------------------------------------------
# Compute Bayesian PO thresholds k_b calibrated to the OBF bounds
# (using a non-informative prior, n_a = 0)
# -----------------------------------------------------------------------------
k_b = PO_sequential(k_l = obf, n_l = n_l, mu_a = 0, n_a = 0, compute_threshold = T)


# =============================================================================
# Simulation Study 1: Fixed Analysis Prior (No Borrowing)
# =============================================================================
# Parameters:
#   mu_a  : prior mean (set to the assumed true effect)
#   n_a   : prior ESS (moderate informativeness)
#   theta : true treatment effect used to generate data
# =============================================================================

mu_a = 0.3
n_a = 30
theta = 0.3

res = list()
for(i in 1:10000){
  
  set.seed(i)
  
  # Generate a full sample of size max(n_l) from N(theta, sigma2)
  x = rnorm(max(n_l), mean = theta, sd = sqrt(sigma2))
  
  # Compute cumulative MLEs at each stage
  mle_l = c()
  for(l in 1:3){
    mle_l[l] = mean(x[1:n_l[l]])
  }
  
  # Run sequential PO computation with fixed analysis prior
  r = PO_sequential(mle_l = mle_l, k_l = obf, n_l = n_l, mu_a = mu_a, n_a = n_a)
  
  res$po_values[[i]] = round(r$po,3)
  res$k_f_values[[i]] = round(r$k_f,3)
  res$k_b_values[[i]] = round(r$k_b,3)
  res$c_b_values[[i]] = round(r$c_b,3)
  res$freq_power[[i]] = r$po > r$k_f   # Frequentist-calibrated Bayesian decision
  res$bayes_power[[i]] = r$po > r$k_b  # Pure Bayesian decision (non-informative threshold)
}

# Collect simulation results into matrices (rows = replications, cols = stages)
po_mat      <- do.call(rbind, res$po_values)
k_f_mat     <- do.call(rbind, res$k_f_values)
k_b_mat     <- do.call(rbind, res$k_b_values)
c_b_mat     <- do.call(rbind, res$c_b_values)
freq_mat    <- do.call(rbind, res$freq_power)
bayes_mat   <- do.call(rbind, res$bayes_power)

# Summary table: median PO, median k_f, k_b, c_b, and empirical power at each stage
tab <- data.frame(
  row.names   = c("Stage 1", "Stage 2", "Stage 3"),
  po          = round(apply(po_mat,    2, median), 3),
  k_f         = round(apply(k_f_mat,   2, median), 3),
  k_b         = round(apply(k_b_mat,   2, mean),   3),
  c_b         = round(apply(c_b_mat,   2, mean),   3),
  freq_power  = round(apply(freq_mat,  2, mean),   3),
  bayes_power = round(apply(bayes_mat, 2, mean),   3)
)

t(tab)


# =============================================================================
# Simulation Study 2: Dynamic Borrowing from Historical Data
# =============================================================================
# Four historical datasets are defined, with varying means and sample sizes,
# representing both consistent and inconsistent historical evidence.
#
# historical : list of c(x_h, n_h) vectors
#   - Datasets 1-2: consistent with theta = 0.3 (same as current true effect)
#   - Datasets 3-4: inconsistent (mean = 0, i.e., no effect)
#
# tau   : heterogeneity penalization (tau = 0 means borrowing depends solely
#         on similarity, without additional variance inflation)
# theta : true treatment effect for data generation
# =============================================================================

theta_0 = 0
sigma2 = 4
tau = "upper_bound" # tau = "upper_bound" to set tau^2 = sigma^2/n_l (adaptive, 
                    # estimated from current data)
                    # tau = 0 to set tau^2 = 0 
                    # tau = numeric > 0  to set a fixed heterogeneity penalty

historical = list(c(0.3, 100),   # historical dataset 1: consistent, n=100
                  c(0.3, 300),   # historical dataset 2: consistent, n=300
                  c(0,   100),   # historical dataset 3: inconsistent, n=100
                  c(0,   300))   # historical dataset 4: inconsistent, n=300
theta = 0.3

res = list()
for(i in 1:10000){
  
  set.seed(i)
  
  x = rnorm(max(n_l), mean = theta, sd = sqrt(sigma2))
  mle_l = c()
  for(l in 1:length(n_l)){
    mle_l[l] = mean(x[1:n_l[l]])
  }
  
  # Run sequential PO with dynamic borrowing
  po = PO_sequential(mle_l = mle_l, k_l = obf, n_l = n_l, mu_a = 0, n_a = 0, historical = historical, tau = tau)
  
  po$freq_decision  = po$po > po$k_f
  po$bayes_decision = po$po > po$k_b
  
  res[[i]] = po
}

# Aggregate simulation results across replications
median_po   = apply(sapply(res, function(x) x$po),    1, median)
median_k_f  = apply(sapply(res, function(x) x$k_f),   1, median)
k_b         = rowMeans(sapply(res, function(x) x$k_b))
mean_mu_a   = rowMeans(sapply(res, function(x) x$mu_a))
mean_n_a    = rowMeans(sapply(res, function(x) x$n_a))
mean_tau    = rowMeans(sapply(res, function(x) x$tau_sq))

# Average borrowing weights per historical source (rows) and stage (cols)
mean_w = apply(simplify2array(lapply(res, function(x) x$w)), c(1,2), mean, na.rm = TRUE)
rownames(mean_w) = paste0("w_", 1:nrow(mean_w))

# Average similarity scores per historical source and stage
mean_s = apply(simplify2array(lapply(res, function(x) x$s)), c(1,2), mean, na.rm = TRUE)
rownames(mean_s) = paste0("s_", 1:nrow(mean_s))

# Average ESS contributions per historical source and stage
mean_ess = apply(simplify2array(lapply(res, function(x) x$ess)), c(1,2), mean, na.rm = TRUE)
rownames(mean_ess) = paste0("ess_", 1:nrow(mean_ess))

mean_freq_power  = rowMeans(sapply(res, function(x) x$freq_decision))
mean_bayes_power = rowMeans(sapply(res, function(x) x$bayes_decision))

# Compile results into a single summary matrix
res_borrowing = rbind(median_po,
                      median_k_f,
                      k_b, 
                      mean_mu_a,
                      mean_n_a,
                      mean_tau,
                      mean_w,
                      mean_s,
                      mean_ess,
                      mean_freq_power,
                      mean_bayes_power)
colnames(res_borrowing) = c("Stage 1", "Stage 2", "Stage 3")

# Print results 
round(res_borrowing, 3)

