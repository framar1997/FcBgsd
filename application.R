# Reproduce results of Section 5 (application) of the main manuscript
# and of the corresponding Web Appendix

# -----------------------------------------------------------------------------
# Function: dh
# -----------------------------------------------------------------------------
# Computes the dissimilarity between the current study data and a historical
# dataset using a Hellinger-based distance between two Gaussian
# distributions with equal variance but different sample sizes.
#
# Arguments:
#   x_c  : observed mean in the current sample at the current stage
#           (on the logit scale)
#   x_h  : observed mean in the historical dataset (on the logit scale)
#   s_c  : standard deviation of the current sample ( = sqrt(sigma2))
#   s_h  : standard deviation of the historical dataset ( = sqrt(sigma2))
#   n_c  : sample size of the current study at the current stage
#   n_h  : sample size of the historical dataset
#
# Returns:
#   A scalar in [0, 1], where 0 indicates perfect similarity and 1 indicates
#   maximum dissimilarity. Internally, s_c and s_h are divided by sqrt(n) to
#   obtain standard errors, and the square root of the Hellinger distance is
#   returned.
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
# of the analysis prior for theta.
#
# Arguments:
#   t         : standardized test statistic at the current stage,
#               t = sqrt(n / nu2_theta) * (mle_theta - theta_0)
#   n         : cumulative (total) sample size at the current stage
#   mu_a      : mean of the analysis prior for theta
#   n_a       : effective sample size (ESS) of the analysis prior;
#               n_a = 0 corresponds to a non-informative prior
#   nu2_theta : variance of the MLE of theta = logit(eta0) - logit(eta1),
#               equal to 4 * sigma2 (delta method approximation)
#
# Returns:
#   A positive scalar: the posterior odds in favor of H_1.
# -----------------------------------------------------------------------------
PO = function(t, n, mu_a, n_a, nu2_theta){
  
  t_a = sqrt(n_a)*(mu_a - theta_0)/sqrt(nu2_theta)
  
  t_post = (sqrt(n_a)*t_a + sqrt(n)*t) / sqrt(n_a + n)
  
  omega_post = exp(pnorm(t_post, log.p = TRUE) - pnorm(t_post, lower.tail = FALSE, log.p = TRUE))
  
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
#   kb        : target threshold on the Posterior Odds scale
#   n         : cumulative sample size at the current stage
#   mu_a      : mean of the analysis prior for theta
#   n_a       : effective sample size (ESS) of the analysis prior
#   nu2_theta : variance of the MLE of theta (see PO for definition)
#
# Returns:
#   k_f: the standardized test statistic value such that PO(k_f, n, mu_a, n_a) = kb.
# -----------------------------------------------------------------------------
PO_inv = function(kb, n, mu_a, n_a, nu2_theta){
  
  t_a = sqrt(n_a)*(mu_a - theta_0)/sqrt(nu2_theta)
  
  kf = sqrt((n_a + n)/n)*qnorm(kb/(kb+1)) - sqrt(n_a/n)*t_a
  
  kf
}


# -----------------------------------------------------------------------------
# Function: update_analysis_prior
# -----------------------------------------------------------------------------
# Updates the analysis prior parameters (mu_a, n_a) for a single arm (either
# eta0 or eta1) at stage l >= 2, by borrowing information from historical
# datasets for that arm. Borrowing is weighted by the similarity between the
# current arm data and each historical dataset (via dh), and may be penalized
# by a heterogeneity parameter tau.
#
# If no historical data are provided, the prior parameters are returned
# unchanged.
#
# Arguments:
#   l          : current stage index (integer >= 1)
#   mu_a       : prior mean at the previous stage (on the logit scale)
#   n_a        : prior ESS at the previous stage
#   current    : vector c(x_l, n_l): observed logit-scale mean and arm-level
#                sample size at stage l-1
#   historical : list of vectors c(x_h, n_h), one per historical dataset for
#                this arm (on the logit scale); if NULL, no borrowing is performed
#   tau        : heterogeneity penalization parameter:
#                  - numeric value >= 0: fixed penalization
#                    (tau = 0 means no extra penalization beyond similarity)
#                  - "upper_bound" (non-numeric): tau is set adaptively as
#                    sigma2 / n_l (variance of the current arm MLE)
#
# Returns:
#   A named list with:
#     mu_a : updated prior mean (weighted average of historical logit means)
#     n_a  : updated prior ESS
#     w    : vector of weights assigned to each historical dataset
#     s    : vector of similarity scores for each historical dataset
#     ess  : vector of ESS contributions per historical dataset
# -----------------------------------------------------------------------------
update_analysis_prior = function(l, mu_a, n_a, current = NULL, historical = NULL, tau = 0){
  if(is.null(historical)){
    list(mu_a = as.numeric(mu_a), n_a = as.numeric(n_a), w = NULL)
  } else{ 
    x_l = current[1]
    n_l = current[2]
    x_l_plogis = plogis(current[1])
    
    x_h  <- sapply(historical, `[`, 1)
    n_h <- sapply(historical, `[`, 2) 
    
    x_h_plogis = sapply(x_h, function(x) plogis(x))
    
    if(is.numeric(tau)){
      simil = apply(rbind(x_h, n_h, x_h_plogis), 2, function(x){
        s = (1 - round(dh(x_l, as.numeric(x[1]), sqrt(sigma2), sqrt(sigma2), n_l, as.numeric(x[2])), 3))
        w = s*x[2]/(sigma2 + s*x[2]*tau)
        ess = s*x[2]
        rbind(w,s,ess)
      })
      w = simil[1,]
      s = simil[2,]
      ess = simil[3,]
      
      mu_a = if(sum(w) == 0) mu_a else sum(w*x_h)/sum(w)
      n_a = sigma2/(1/sum(w) + tau)
    } else{
      # Adaptive tau: estimated as the variance of the current arm MLE
      tau_sq = sigma2/n_l
      simil = apply(rbind(x_h, n_h, x_h_plogis), 2, function(x){
        s = (1 - round(dh(x_l, as.numeric(x[1]), sqrt(sigma2), sqrt(sigma2), n_l, as.numeric(x[2])), 3))
        w = s*x[2]/(sigma2 + s*x[2]*tau_sq)
        ess = s*x[2]
        rbind(w,s,ess)
      })
      w = simil[1,]
      s = simil[2,]
      ess = simil[3,]
      
      mu_a = if(sum(w) == 0) mu_a else sum(w*x_h)/sum(w)
      n_a = sigma2/(1/sum(w) + tau_sq)
    }
    
    list(mu_a = as.numeric(mu_a), n_a = as.numeric(n_a), w = w, s = s, ess = ess)
  }
}


# -----------------------------------------------------------------------------
# Function: combine_priors_theta
# -----------------------------------------------------------------------------
# Combines the separate analysis priors for eta0 and eta1 into a single
# Normal prior for the treatment effect theta = logit(eta0) - logit(eta1),
# using the delta method (independence approximation).
#
# Arguments:
#   prior_eta0 : named list with mu_a and n_a for the control arm prior
#   prior_eta1 : named list with mu_a and n_a for the treatment arm prior
#
# Returns:
#   A named list with:
#     mu_a : prior mean for theta = mu_a0 - mu_a1
#     n_a  : prior ESS for theta, computed as 4 / (1/n_a0 + 1/n_a1)
#            (accounts for the variance of the difference on the logit scale)
# -----------------------------------------------------------------------------
combine_priors_theta = function(prior_eta0, prior_eta1) {
  
  mu_a_theta  = prior_eta0$mu_a - prior_eta1$mu_a
  
  n_a_theta = 4*((1/prior_eta0$n_a + 1/prior_eta1$n_a)^(-1))
  
  list(mu_a = mu_a_theta, n_a = n_a_theta)
}


# -----------------------------------------------------------------------------
# Function: PO_sequential
# -----------------------------------------------------------------------------
# Computes the Posterior Odds sequentially across all planned interim stages
# for the two-arm log-odds setting. Can operate in two modes:
#
#   (1) Trial mode (compute_threshold = FALSE):
#       Given vectors of observed MLEs for theta, eta0, and eta1, computes
#       the PO at each stage. At each stage l >= 2, the arm-level priors are
#       updated via dynamic borrowing, then combined into a prior for theta.
#
#   (2) Threshold mode (compute_threshold = TRUE):
#       Given a vector of frequentist bounds (k_l), computes the corresponding
#       Bayesian PO thresholds k_b without observed data.
#
# Arguments:
#   mle_l           : vector of cumulative MLEs for theta at each stage;
#                     required when compute_threshold = FALSE
#   mle_eta0_l      : vector of cumulative MLEs for logit(eta0) at each stage
#   mle_eta1_l      : vector of cumulative MLEs for logit(eta1) at each stage
#   k_l             : vector of frequentist decision bounds (e.g., OBF bounds)
#   n_l             : vector of cumulative total sample sizes at each stage
#                     (each arm contributes n_l[l]/2 observations)
#   mu_a, n_a       : initial prior mean and ESS for theta at stage 1
#   mu_a0, n_a0     : initial prior mean and ESS for logit(eta0) at stage 1
#   mu_a1, n_a1     : initial prior mean and ESS for logit(eta1) at stage 1
#   compute_threshold : logical; if FALSE (default), computes PO from observed
#                       MLEs; if TRUE, computes PO thresholds from k_l bounds
#   historical_eta0 : list of historical datasets for the control arm,
#                     each c(logit_mean, n); passed to update_analysis_prior
#   historical_eta1 : list of historical datasets for the treatment arm,
#                     each c(logit_mean, n); passed to update_analysis_prior
#   tau             : heterogeneity penalization; passed to update_analysis_prior:
#                       - tau = 0            : no extra penalization (tau^2 = 0)
#                       - tau = "upper_bound": adaptive, tau^2 = sigma^2 / n_l
#                       - tau = numeric > 0  : fixed heterogeneity penalty
#   nu2_theta       : variance of the MLE of theta ( = 4 * sigma2)
#
# Returns (compute_threshold = FALSE):
#   A named list with:
#     po                    : vector of PO values at each stage
#     mu_a, n_a             : vectors of prior parameters for theta at each stage
#     mu_a0, n_a0           : vectors of prior parameters for eta0 at each stage
#     mu_a1, n_a1           : vectors of prior parameters for eta1 at each stage
#     w_eta0, s_eta0, ess_eta0 : matrices of borrowing weights, similarity scores,
#                                and ESS contributions for the control arm
#     w_eta1, s_eta1, ess_eta1 : same for the treatment arm
#     k_b                   : vector of Bayesian PO thresholds (non-informative prior)
#     k_f                   : vector of Bayesian PO thresholds (analysis prior)
#     c_b                   : vector of frequentist MLE thresholds via PO_inv
#
# Returns (compute_threshold = TRUE):
#   A named vector k_b: the PO value at each stage when evaluated at k_l.
# -----------------------------------------------------------------------------
PO_sequential <- function(mle_l = NULL, k_l = NULL, n_l, 
                          mu_a, n_a,
                          mu_a0 = 0, n_a0 = 0,
                          mu_a1 = 0, n_a1 = 0,
                          compute_threshold = F, 
                          historical_eta0 = NULL,
                          historical_eta1 = NULL,
                          mle_eta0_l = NULL,
                          mle_eta1_l = NULL,
                          tau = 0,
                          nu2_theta) {
  
  weight_eta0 = matrix(nrow = max(1, length(historical_eta0)), ncol = length(n_l))
  weight_eta1 = matrix(nrow = max(1, length(historical_eta1)), ncol = length(n_l))
  s_eta0   = matrix(nrow = max(1, length(historical_eta0)), ncol = length(n_l))
  s_eta1   = matrix(nrow = max(1, length(historical_eta1)), ncol = length(n_l))
  ess_eta0 = matrix(nrow = max(1, length(historical_eta0)), ncol = length(n_l))
  ess_eta1 = matrix(nrow = max(1, length(historical_eta1)), ncol = length(n_l))
  
  po = c()
  mu_a_l  = c(mu_a);   n_a_l  = c(n_a)
  mu_a0_l = c(mu_a0);  n_a0_l = c(n_a0)
  mu_a1_l = c(mu_a1);  n_a1_l = c(n_a1)
  
  if (compute_threshold == F) {
    
    t_l = c()
    t_l[1] = sqrt(n_l[1] / nu2_theta) * (mle_l[1] - theta_0)
    po[1]  = PO(t_l[1], n_l[1], mu_a_l[1], n_a_l[1], nu2_theta = nu2_theta)
    
    k_b = c(PO(k_l[1], n_l[1], mu_a_l[1], n_a = 0, nu2_theta = nu2_theta))
    k_f = c(PO(k_l[1], n_l[1], mu_a_l[1], n_a_l[1], nu2_theta = nu2_theta))
    c_b = c(PO_inv(k_b[1], n_l[1], mu_a_l[1], n_a_l[1], nu2_theta = nu2_theta))
    
    for (l in 2:length(n_l)) {
      
      # Update prior for eta0 using historical control arm data
      prior_eta0 = update_analysis_prior(l,
                                         mu_a = mu_a0_l[l-1], n_a = n_a0_l[l-1],
                                         current = c(mle_eta0_l[l-1], n_l[l-1]/2),
                                         historical = historical_eta0, tau = tau)
      mu_a0_l[l] = prior_eta0$mu_a
      n_a0_l[l]  = prior_eta0$n_a
      weight_eta0[, l] = prior_eta0$w
      s_eta0[, l]      = prior_eta0$s
      ess_eta0[, l]    = prior_eta0$ess
      
      # Update prior for eta1 using historical treatment arm data
      prior_eta1 = update_analysis_prior(l,
                                         mu_a = mu_a1_l[l-1], n_a = n_a1_l[l-1],
                                         current = c(mle_eta1_l[l-1], n_l[l-1]/2),
                                         historical = historical_eta1, tau = tau)
      mu_a1_l[l] = prior_eta1$mu_a
      n_a1_l[l]  = prior_eta1$n_a
      weight_eta1[, l] = prior_eta1$w
      s_eta1[, l]      = prior_eta1$s
      ess_eta1[, l]    = prior_eta1$ess
      
      # Combine arm-level priors into a single prior for theta
      prior_theta = combine_priors_theta(prior_eta0, prior_eta1)
      mu_a_l[l]  = prior_theta$mu_a
      n_a_l[l]   = prior_theta$n_a
      
      t_l[l] = sqrt(n_l[l] / nu2_theta) * (mle_l[l] - theta_0)
      po[l]  = PO(t_l[l], n_l[l], mu_a_l[l], n_a_l[l], nu2_theta = nu2_theta)
      
      k_b = c(k_b, PO(k_l[l],     n_l[l], mu_a_l[l], n_a = 0, nu2_theta = nu2_theta))
      k_f = c(k_f, PO(k_l[l],     n_l[l], mu_a_l[l], n_a_l[l], nu2_theta = nu2_theta))
      c_b = c(c_b, PO_inv(k_b[l], n_l[l], mu_a_l[l], n_a_l[l], nu2_theta = nu2_theta))
    }
    
    return(list(
      po      = po,
      mu_a    = mu_a_l,  n_a    = n_a_l,
      mu_a0   = mu_a0_l, n_a0   = n_a0_l,
      mu_a1   = mu_a1_l, n_a1   = n_a1_l,
      w_eta0  = weight_eta0, s_eta0 = s_eta0, ess_eta0 = ess_eta0,
      w_eta1  = weight_eta1, s_eta1 = s_eta1, ess_eta1 = ess_eta1,
      k_b     = k_b, k_f = k_f, c_b = c_b
    ))
  }
  
  if (compute_threshold == T) {
    
    po[1] = PO(k_l[1], n_l[1], mu_a_l[1], n_a_l[1], nu2_theta = nu2_theta)
    mu_a0_l <- mu_a1_l <-  rep(0, length(n_l))
    n_a0_l <- n_a1_l <- rep(0, length(n_l))
    for(l in 2:length(n_l)){
      
      prior_eta0 = update_analysis_prior(l,
                                         mu_a = mu_a0_l[l-1], n_a = n_a0_l[l-1],
                                         current = NULL, historical = historical_eta0,
                                         tau = tau)
      prior_eta1 = update_analysis_prior(l,
                                         mu_a = mu_a1_l[l-1], n_a = n_a1_l[l-1],
                                         current = NULL, historical = historical_eta1,
                                         tau = tau)
      
      prior_theta = combine_priors_theta(prior_eta0, prior_eta1)
      mu_a_l[l]  = prior_theta$mu_a
      n_a_l[l]   = prior_theta$n_a
      
      po[l] = PO(k_l[l], n_l[l], mu_a_l[l], n_a_l[l], nu2_theta = nu2_theta)
    }
    
    return(rbind(k_b = po))
  }
}


# =============================================================================
# Application Analysis
# =============================================================================

library(gsDesign)

# -----------------------------------------------------------------------------
# Global parameters
# sigma2    : known variance on the logit scale (approximation)
# nu2_theta : variance of the MLE of theta = logit(eta0) - logit(eta1),
#             equal to 4 * sigma2 by the delta method (two equal-sized arms)
# theta_0   : null hypothesis value (superiority trial: theta_0 = 0)
# n_l       : vector of cumulative total sample sizes at each planned stage
# -----------------------------------------------------------------------------
sigma2    = 4
nu2_theta = 4 * sigma2
theta_0   = 0
n_l       = c(100, 200, 300, 400)

# -----------------------------------------------------------------------------
# O'Brien-Fleming group sequential bounds (one-sided, alpha = 0.05, 80% power)
# These serve as the frequentist reference thresholds k_l at each stage.
# -----------------------------------------------------------------------------
fr <- gsDesign(k = 4,          # number of interim stages
               test.type = 1,  # one-sided test
               alpha = 0.05, 
               beta = 0.2,
               sfu = sfLDOF)   # Lan-DeMets O'Brien-Fleming

obf = fr$upper$bound

# -----------------------------------------------------------------------------
# Historical data: incidence rates of the binary endpoint in the control arm
# (eta0) and one historical observation for the treatment arm (eta1).
# All rates are converted to the logit scale for use in the prior.
#
# eta_0       : observed event rates in 3 historical control arm studies
#               (hypothesized phase II SSTARLET trial, adult trial, and
#               pediatric trial; the pediatric estimate uses a continuity
#               correction since 0 events were observed)
# n_0         : sample sizes of the 3 historical control arm studies
# eta_1       : observed event rate in 1 historical treatment arm study
#               (hypothesized phase II SSTARLET trial)
# n_1         : sample size of the historical treatment arm study
# -----------------------------------------------------------------------------
eta_0 = c(2/2000, 4/3443, (0+0.5)/(358+1))
n_0   = c(2000, 3443, 358)
eta_1 = 1/2000
n_1   = 2000

# Convert to logit scale
logit_eta0 = log(eta_0 / (1 - eta_0))
logit_eta1 = log(eta_1 / (1 - eta_1))

# Historical datasets for each arm (logit mean, sample size)
historical_eta0 = list(
  c(logit_eta0[1], n_0[1]),
  c(logit_eta0[2], n_0[2]),
  c(logit_eta0[3], n_0[3])
)

historical_eta1 = list(
  c(logit_eta1, n_1)
)

# -----------------------------------------------------------------------------
# Initial (non-informative) priors for eta0, eta1, and theta at stage 1
# -----------------------------------------------------------------------------
mu_a0_init = 0;  n_a0_init = 0
mu_a1_init = 0;  n_a1_init = 0
mu_a_init  = mu_a0_init - mu_a1_init   # = 0
n_a_init   = 0

# -----------------------------------------------------------------------------
# Hypothesized true values for the simulation
# eta_0_true : assumed true event rate in the control arm
# eta_1_true : assumed true event rate in the treatment arm (40% reduction)
# -----------------------------------------------------------------------------
eta_0_true = 0.001
eta_1_true = 0.6 * eta_0_true


# =============================================================================
# Simulation Study 1: tau^2 = 0 
# =============================================================================

tau = 0  # tau = 0            : tau^2 = 0 (maximum borrowing given similarity)
# tau = "upper_bound": adaptive, tau^2 = sigma^2 / n_l
# tau = numeric > 0  : fixed heterogeneity penalty

res = list()

for(i in 1:10000){
  
  set.seed(i)
  
  # Generate arm-level data on the logit scale
  x_eta0 = rnorm(max(n_l/2), mean = qlogis(eta_0_true), sd = sqrt(sigma2))
  x_eta1 = rnorm(max(n_l/2), mean = qlogis(eta_1_true), sd = sqrt(sigma2))
  
  mle_l      = c()   # MLE of theta = logit(eta0) - logit(eta1) at each stage
  mle_eta0_l = c()
  mle_eta1_l = c()
  
  for(l in 1:length(n_l)){
    n_per_arm      = n_l[l] / 2
    mle_eta0_l[l]  = mean(x_eta0[1:n_per_arm])
    mle_eta1_l[l]  = mean(x_eta1[1:n_per_arm])
    mle_l[l]       = mle_eta0_l[l] - mle_eta1_l[l]
  }
  
  # Run sequential PO with dynamic borrowing (tau = 0)
  po = PO_sequential(
    mle_l           = mle_l,
    mle_eta0_l      = mle_eta0_l,
    mle_eta1_l      = mle_eta1_l,
    k_l             = obf,
    n_l             = n_l,
    mu_a            = mu_a_init,
    n_a             = n_a_init,
    mu_a0           = mu_a0_init,
    n_a0            = n_a0_init,
    mu_a1           = mu_a1_init,
    n_a1            = n_a1_init,
    historical_eta0 = historical_eta0,
    historical_eta1 = historical_eta1,
    tau             = tau,
    nu2_theta       = nu2_theta
  )
  
  po$freq_decision  = po$po > po$k_f
  po$bayes_decision = po$po > po$k_b
  
  res[[i]] = po
}

# Aggregate simulation results across replications
median_po    = apply(sapply(res, function(x) x$po), 1, median)
mean_po      = rowMeans(sapply(res, function(x) x$po))
mean_k_f     = rowMeans(sapply(res, function(x) x$k_f))
median_k_f   = apply(sapply(res, function(x) x$k_f), 1, median)
k_b          = rowMeans(sapply(res, function(x) x$k_b))
mean_mu_a    = rowMeans(sapply(res, function(x) x$mu_a))
mean_n_a     = rowMeans(sapply(res, function(x) x$n_a))
mean_mu_a0   = rowMeans(sapply(res, function(x) x$mu_a0))
mean_n_a0    = rowMeans(sapply(res, function(x) x$n_a0))
mean_mu_a1   = rowMeans(sapply(res, function(x) x$mu_a1))
mean_n_a1    = rowMeans(sapply(res, function(x) x$n_a1))

# Average borrowing weights, similarity scores, and ESS for the control arm
mean_w_eta0   = apply(simplify2array(lapply(res, function(x) x$w_eta0)),   c(1,2), mean, na.rm = TRUE)
mean_s_eta0   = apply(simplify2array(lapply(res, function(x) x$s_eta0)),   c(1,2), mean, na.rm = TRUE)
mean_ess_eta0 = apply(simplify2array(lapply(res, function(x) x$ess_eta0)), c(1,2), mean, na.rm = TRUE)
rownames(mean_w_eta0)   = paste0("w0_",   1:nrow(mean_w_eta0))
rownames(mean_s_eta0)   = paste0("s0_",   1:nrow(mean_s_eta0))
rownames(mean_ess_eta0) = paste0("ess0_", 1:nrow(mean_ess_eta0))

# Average borrowing weights, similarity scores, and ESS for the treatment arm
mean_w_eta1   = apply(simplify2array(lapply(res, function(x) x$w_eta1)),   c(1,2), mean, na.rm = TRUE)
mean_s_eta1   = apply(simplify2array(lapply(res, function(x) x$s_eta1)),   c(1,2), mean, na.rm = TRUE)
mean_ess_eta1 = apply(simplify2array(lapply(res, function(x) x$ess_eta1)), c(1,2), mean, na.rm = TRUE)
rownames(mean_w_eta1)   = paste0("w1_",   1:nrow(mean_w_eta1))
rownames(mean_s_eta1)   = paste0("s1_",   1:nrow(mean_s_eta1))
rownames(mean_ess_eta1) = paste0("ess1_", 1:nrow(mean_ess_eta1))

mean_freq_power  = rowMeans(sapply(res, function(x) x$freq_decision))
mean_bayes_power = rowMeans(sapply(res, function(x) x$bayes_decision))

# Compile results into a single summary matrix
case_tau0 = rbind(
  median_po,
  median_k_f,
  k_b,
  mean_mu_a,   mean_n_a,
  mean_mu_a0,  mean_n_a0,
  mean_mu_a1,  mean_n_a1,
  mean_w_eta0, mean_s_eta0, mean_ess_eta0,
  mean_w_eta1, mean_s_eta1, mean_ess_eta1,
  mean_freq_power,
  mean_bayes_power
)
options(scipen = 999)

# Print results (excluding Stage 1, as borrowing starts from Stage 2)
case_tau0 = round(case_tau0[, -1], 3)
colnames(case_tau0) = paste("Stage", 1:ncol(case_tau0))

# =============================================================================
# Simulation Study 2: tau^2 = sigma^2 / n_l 
# =============================================================================

tau = "upper_bound"  # tau = "upper_bound": adaptive, tau^2 = sigma^2 / n_l
# tau = 0            : tau^2 = 0
# tau = numeric > 0  : fixed heterogeneity penalty

res = list()

for(i in 1:10000){
  
  set.seed(i)
  
  # Generate arm-level data on the logit scale
  x_eta0 = rnorm(max(n_l/2), mean = qlogis(eta_0_true), sd = sqrt(sigma2))
  x_eta1 = rnorm(max(n_l/2), mean = qlogis(eta_1_true), sd = sqrt(sigma2))
  
  mle_l      = c()
  mle_eta0_l = c()
  mle_eta1_l = c()
  
  for(l in 1:length(n_l)){
    n_per_arm      = n_l[l] / 2
    mle_eta0_l[l]  = mean(x_eta0[1:n_per_arm])
    mle_eta1_l[l]  = mean(x_eta1[1:n_per_arm])
    mle_l[l]       = mle_eta0_l[l] - mle_eta1_l[l]
  }
  
  # Run sequential PO with dynamic borrowing (tau = "upper_bound")
  po = PO_sequential(
    mle_l           = mle_l,
    mle_eta0_l      = mle_eta0_l,
    mle_eta1_l      = mle_eta1_l,
    k_l             = obf,
    n_l             = n_l,
    mu_a            = mu_a_init,
    n_a             = n_a_init,
    mu_a0           = mu_a0_init,
    n_a0            = n_a0_init,
    mu_a1           = mu_a1_init,
    n_a1            = n_a1_init,
    historical_eta0 = historical_eta0,
    historical_eta1 = historical_eta1,
    tau             = tau,
    nu2_theta       = nu2_theta
  )
  
  po$freq_decision  = po$po > po$k_f
  po$bayes_decision = po$po > po$k_b
  
  res[[i]] = po
}

# Aggregate simulation results across replications
median_po    = apply(sapply(res, function(x) x$po), 1, median)
mean_po      = rowMeans(sapply(res, function(x) x$po))
mean_k_f     = rowMeans(sapply(res, function(x) x$k_f))
median_k_f   = apply(sapply(res, function(x) x$k_f), 1, median)
k_b          = rowMeans(sapply(res, function(x) x$k_b))
mean_mu_a    = rowMeans(sapply(res, function(x) x$mu_a))
mean_n_a     = rowMeans(sapply(res, function(x) x$n_a))
mean_mu_a0   = rowMeans(sapply(res, function(x) x$mu_a0))
mean_n_a0    = rowMeans(sapply(res, function(x) x$n_a0))
mean_mu_a1   = rowMeans(sapply(res, function(x) x$mu_a1))
mean_n_a1    = rowMeans(sapply(res, function(x) x$n_a1))

mean_w_eta0   = apply(simplify2array(lapply(res, function(x) x$w_eta0)),   c(1,2), mean, na.rm = TRUE)
mean_s_eta0   = apply(simplify2array(lapply(res, function(x) x$s_eta0)),   c(1,2), mean, na.rm = TRUE)
mean_ess_eta0 = apply(simplify2array(lapply(res, function(x) x$ess_eta0)), c(1,2), mean, na.rm = TRUE)
rownames(mean_w_eta0)   = paste0("w0_",   1:nrow(mean_w_eta0))
rownames(mean_s_eta0)   = paste0("s0_",   1:nrow(mean_s_eta0))
rownames(mean_ess_eta0) = paste0("ess0_", 1:nrow(mean_ess_eta0))

mean_w_eta1   = apply(simplify2array(lapply(res, function(x) x$w_eta1)),   c(1,2), mean, na.rm = TRUE)
mean_s_eta1   = apply(simplify2array(lapply(res, function(x) x$s_eta1)),   c(1,2), mean, na.rm = TRUE)
mean_ess_eta1 = apply(simplify2array(lapply(res, function(x) x$ess_eta1)), c(1,2), mean, na.rm = TRUE)
rownames(mean_w_eta1)   = paste0("w1_",   1:nrow(mean_w_eta1))
rownames(mean_s_eta1)   = paste0("s1_",   1:nrow(mean_s_eta1))
rownames(mean_ess_eta1) = paste0("ess1_", 1:nrow(mean_ess_eta1))

mean_freq_power  = rowMeans(sapply(res, function(x) x$freq_decision))
mean_bayes_power = rowMeans(sapply(res, function(x) x$bayes_decision))

case_taumax = rbind(
  mean_po,
  mean_k_f,
  k_b,
  mean_mu_a,   mean_n_a,
  mean_mu_a0,  mean_n_a0,
  mean_mu_a1,  mean_n_a1,
  mean_w_eta0, mean_s_eta0, mean_ess_eta0,
  mean_w_eta1, mean_s_eta1, mean_ess_eta1,
  mean_freq_power,
  mean_bayes_power
)
options(scipen = 999)

case_taumax = round(case_taumax[, -1], 3)
colnames(case_taumax) = c(paste("Stage", 1:ncol(case_taumax)))

# =============================================================================
# Final results summary
# =============================================================================

case_tau0
case_taumax


# =============================================================================
# ERR Power Curves
# =============================================================================
# For each value of epsilon in a grid, data are generated with eta_1 shifted
# by epsilon relative to eta_0, and the empirical ERR_f and ERR_b are computed
# at each stage. Results are saved to disk and then plotted.
#
# epsilon_grid : grid of effect sizes (differences on the probability scale)
# n_sim        : number of Monte Carlo replications per epsilon value
# =============================================================================

n_sim        = 10000
epsilon_grid = seq(0, 0.0008, 0.00005)


# -----------------------------------------------------------------------------
# Power curves: tau^2 = 0
# -----------------------------------------------------------------------------

tau = 0

freq_power  = matrix(NA, nrow = length(epsilon_grid), ncol = length(n_l))
bayes_power = matrix(NA, nrow = length(epsilon_grid), ncol = length(n_l))

for(e in seq_along(epsilon_grid)){
  
  epsilon = epsilon_grid[e]
  
  freq_stop  = numeric(length(n_l))
  bayes_stop = numeric(length(n_l))
  
  theta_0 = 0
  for(i in 1:n_sim){
    
    set.seed(i)
    
    x_eta0 = rnorm(max(n_l/2), mean = qlogis(eta_0_true),           sd = sqrt(sigma2))
    x_eta1 = rnorm(max(n_l/2), mean = qlogis(eta_0_true - epsilon), sd = sqrt(sigma2))
    
    mle_l      = c()
    mle_eta0_l = c()
    mle_eta1_l = c()
    
    for(l in 1:length(n_l)){
      n_per_arm      = n_l[l] / 2
      mle_eta0_l[l]  = mean(x_eta0[1:n_per_arm])
      mle_eta1_l[l]  = mean(x_eta1[1:n_per_arm])
      mle_l[l]       = mle_eta0_l[l] - mle_eta1_l[l]
    }
    
    po = PO_sequential(
      mle_l           = mle_l,
      mle_eta0_l      = mle_eta0_l,
      mle_eta1_l      = mle_eta1_l,
      k_l             = obf,
      n_l             = n_l,
      mu_a            = mu_a_init,
      n_a             = n_a_init,
      mu_a0           = mu_a0_init,
      n_a0            = n_a0_init,
      mu_a1           = mu_a1_init,
      n_a1            = n_a1_init,
      historical_eta0 = historical_eta0,
      historical_eta1 = historical_eta1,
      tau             = tau,
      nu2_theta       = nu2_theta
    )
    
    freq_stop  = freq_stop  + as.numeric(po$po > po$k_f)
    bayes_stop = bayes_stop + as.numeric(po$po > po$k_b)
  }
  
  freq_power[e, ]  = freq_stop  / n_sim
  bayes_power[e, ] = bayes_stop / n_sim
}

saveRDS(list(epsilon_grid = epsilon_grid,
             freq_power   = freq_power,
             bayes_power  = bayes_power,
             tau          = tau,
             n_sim        = n_sim),
        file = "power_results_tau0.rds")


# -----------------------------------------------------------------------------
# Power curves: tau^2 = sigma^2 / n_l (adaptive)
# -----------------------------------------------------------------------------

tau = "upper_bound"

freq_power  = matrix(NA, nrow = length(epsilon_grid), ncol = length(n_l))
bayes_power = matrix(NA, nrow = length(epsilon_grid), ncol = length(n_l))

for(e in seq_along(epsilon_grid)){
  
  epsilon = epsilon_grid[e]
  
  freq_stop  = numeric(length(n_l))
  bayes_stop = numeric(length(n_l))
  theta_0 = 0
  for(i in 1:n_sim){
    
    x_eta0 = rnorm(max(n_l/2), mean = qlogis(eta_0_true),           sd = sqrt(sigma2))
    x_eta1 = rnorm(max(n_l/2), mean = qlogis(eta_0_true - epsilon), sd = sqrt(sigma2))
    
    mle_l      = c()
    mle_eta0_l = c()
    mle_eta1_l = c()
    
    for(l in 1:length(n_l)){
      n_per_arm      = n_l[l] / 2
      mle_eta0_l[l]  = mean(x_eta0[1:n_per_arm])
      mle_eta1_l[l]  = mean(x_eta1[1:n_per_arm])
      mle_l[l]       = mle_eta0_l[l] - mle_eta1_l[l]
    }
    
    po = PO_sequential(
      mle_l           = mle_l,
      mle_eta0_l      = mle_eta0_l,
      mle_eta1_l      = mle_eta1_l,
      k_l             = obf,
      n_l             = n_l,
      mu_a            = mu_a_init,
      n_a             = n_a_init,
      mu_a0           = mu_a0_init,
      n_a0            = n_a0_init,
      mu_a1           = mu_a1_init,
      n_a1            = n_a1_init,
      historical_eta0 = historical_eta0,
      historical_eta1 = historical_eta1,
      tau             = tau,
      nu2_theta       = nu2_theta
    )
    
    freq_stop  = freq_stop  + as.numeric(po$po > po$k_f)
    bayes_stop = bayes_stop + as.numeric(po$po > po$k_b)
  }
  
  freq_power[e, ]  = freq_stop  / n_sim
  bayes_power[e, ] = bayes_stop / n_sim
}

saveRDS(list(epsilon_grid = epsilon_grid,
             freq_power   = freq_power,
             bayes_power  = bayes_power,
             tau          = tau,
             n_sim        = n_sim),
        file = "power_results_tauUpperBound.rds")


# -----------------------------------------------------------------------------
# Plot ERR curves (both tau settings, stages 2-4)
# -----------------------------------------------------------------------------

pdf("ERR_curves.pdf", width = 10, height = 7)

col_freq  = "#E41A1C"
col_bayes = "#377EB8"
stages_plot = c(2, 3, 4)

# --- tau^2 = 0 ---
res = readRDS("power_results_tau0.rds")
epsilon_grid = res$epsilon_grid
freq_power   = res$freq_power
bayes_power  = res$bayes_power

par(mfrow = c(2, 3),
    mar  = c(3.2, 3.2, 2.5, 0.8),
    mgp  = c(1.8, 0.5, 0),
    oma  = c(0, 0, 2, 0),
    cex.axis = 1.35,
    cex.lab  = 1.6)

for(s in stages_plot){
  
  plot(epsilon_grid, freq_power[, s],
       type = "l", col = col_freq, lwd = 2,
       ylim = c(0, 1),
       xlab = expression(epsilon),
       ylab = if(s == 2) expression("ERR for" ~ tau^2 == 0) else "",
       main = paste0("Stage ", s),
       las  = 1, bty = "l",
       cex.lab = 1.7)
  
  lines(epsilon_grid, bayes_power[, s], col = col_bayes, lwd = 2)
  
  abline(v = eta_0_true - eta_1_true, lty = 3, col = "grey50")
  abline(h = 0.05, lty = 3, col = "grey50")
}

mtext("ERR curves", outer = TRUE, cex = 1.3, font = 2)

# --- tau^2 = sigma^2 / n_l ---
res = readRDS("power_results_tauUpperBound.rds")
epsilon_grid = res$epsilon_grid
freq_power   = res$freq_power
bayes_power  = res$bayes_power

for(s in stages_plot){
  
  plot(epsilon_grid, freq_power[, s],
       type = "l", col = col_freq, lwd = 2,
       ylim = c(0, 1),
       xlab = expression(epsilon),
       ylab = if(s == 2) expression(paste("ERR for ", tau^2, " = ", sigma^2/n["ℓ"])) else "",
       main = paste0("Stage ", s),
       las  = 1, bty = "l",
       cex.lab = 1.7)
  
  lines(epsilon_grid, bayes_power[, s], col = col_bayes, lwd = 2)
  
  abline(v = eta_0_true - eta_1_true, lty = 3, col = "grey50")
  abline(h = 0.05, lty = 3, col = "grey50")
  
  if(s == 4){
    legend("bottomright",
           legend = c("Frequentist", "Bayesian"),
           col    = c(col_freq, col_bayes),
           lwd    = 2, bty = "n", cex = 1.7)
  }
}

dev.off()