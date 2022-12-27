library(MASS)
library(mvnfast)
library(OwenQ)

# Some silly micro-optimization instead of using one more general log(pnorm(b) - pnorm(a)), 
# should also be more accurate though.
# log_dtruncnorm_left <- function(x, a) {
#   dnorm(x, log = T) - pnorm(a, lower.tail = F, log.p = T)
# }
# 
# log_dtruncnorm_right <- function(x, b) {
#   dnorm(x, log = T) - pnorm(b, log.p = T)
# }

# The next four functions are basically the same ones as before.

risk_reduction_lowest = function(r2,K,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  # risk/K, so if we want overall error < 1e-5, we need risk/K < 1e-5, so the integral error < K * 1e-5
  # Let's do 4 digits?
  risk = integrate(integrand_lowest,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # If abs risk: K-risk
  return(reduction)
}

risk_reduction_exclude = function(r2,K,q,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_t = function(t,u)
  {
    y = dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r2),lower.tail=F)
    return(y)
  }
  # Replacement for integrand_t, should be numerically better.
  inner_integral <- function(x, a, b) {
    if (is.infinite(x)) {
      OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
        # OwenT(x, b) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
        OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
    }
    else {
      OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
        OwenT(x, (a+b*x)/x) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
    }
  }
  
  integrand_u = function(us)
  {
    y = numeric(length(us))
    beta_vec = zq*sqrt(2)-us
    denom = pnorm(beta_vec)
    dnorm_u <- dnorm(us)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    
    a <- (zk - r/sqrt(2) * us) / sqrt(1-r2)
    b <- -(r/sqrt(2)) / sqrt(1-r2)
    
    for (i in seq_along(us))
    {
      u = us[i]
      beta <- beta_vec[i]
      # internal_int = integrate(integrand_t,-Inf,beta,u, abs.tol = 1e-10)$value
      
      internal_int = ifelse(r2 != 1, pnorm(beta) - inner_integral(beta, a[i], b) + inner_integral(-Inf, a[i], b),
                            ifelse(beta < sqrt(2) / r * zk - u, 0, pnorm(beta) - pnorm(sqrt(2) / r * zk - u)))
      numer = dnorm_u[i]*(1-(1-denom[i])^n) * internal_int
      term1 = numer/denom[i]
      
      # Probably doesn't matter, but also numerically more stable
      # term1 <- exp(dnorm(u, log = T) + VGAM::log1mexp(-n*pnorm(beta, lower.tail = F, log.p = T)) - 
      #                pnorm(beta, log.p = T)) * internal_int
      
      # internal_int = integrate(integrand_t,beta,Inf,u, abs.tol = 1e-10)$value
      internal_int <- ifelse(r2 != 1, 1 - pnorm(beta) - inner_integral(Inf, a[i], b) + inner_integral(beta, a[i], b),
                             ifelse(beta < sqrt(2) / r * zk - u, 1-pnorm(sqrt(2) / r * zk - u), 1-pnorm(beta)))
      term2 = dnorm_u[i]*(1-denom[i])^(n-1) * internal_int
      y[i] = term1 + term2
      # print(sprintf("%f", y[i]))
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # K-risk
  return(reduction)
}

risk_reduction_lowest_conditional = function(r2,K,n,qf,qm,relative=T,parental_avg_given=F)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  if (parental_avg_given)
  {
    # It is assumed that what is given is directly the parental average, so that the paternal and maternal quantiles are the same (both representing the quantile of the parental average)
    c = zqf * r/sqrt(2)
  } else {
    c = (zqf+zqm)/2 * r
  }
  baseline = pnorm((zk-c)/sqrt(1-r2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  # risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  risk = integrate(integrand_lowest_cond,-Inf,Inf, abs.tol = baseline * 1e-4)$value
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

risk_reduction_exclude_conditional = function(r2,K,q,n,qf,qm,relative=T)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  gamma = zq*sqrt(2) - c/(r/sqrt(2))
  
  integrand_t = function(t)
  {
    y = dnorm(t)*pnorm((zk-t*r/sqrt(2)-c)/sqrt(1-r2),lower.tail=F)
    return(y)
  }
  
  denom = pnorm(gamma)
  err <- (1e-4 * baseline * denom) / (2 * (1-(1-denom)^n))
  internal_int = integrate(integrand_t,-Inf,gamma, abs.tol = err)$value
  numer = (1-(1-denom)^n) * internal_int
  term1 = numer/denom
  
  err <- (1e-4 * baseline) / (2 * (1-(1-denom)^n))
  internal_int = integrate(integrand_t,gamma,Inf)$value
  term2 = (1-denom)^(n-1) * internal_int
  
  risk = term1 + term2
  
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

# The next two are using monte carlo simulation + control variate

baseline_risk <- function(r2, h2, K, df, dm) {
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  
  posterior = function(gm,gf)
  {
    y = 1
    y = y * dnorm(gm/h)/h
    y = y * dnorm(gf/h)/h
    arg = (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
  
  integrand_gm = function(gms,gf)
  {
    y = numeric(length(gms))
    for (i in seq_along(gms))
    {
      gm = gms[i]
      arg = (zk - (gm+gf)/2) / sqrt(1-h2/2)
      y[i] = pnorm(arg, lower.tail=F)
      post = posterior(gm,gf)
      y[i] = y[i] * post
    }
    return(y)
  }
  
  integrand_gf = function(gfs)
  {
    y = numeric(length(gfs))
    for (i in seq_along(gfs))
    {
      gf = gfs[i]
      y[i] = integrate(integrand_gm,-Inf,Inf,gf)$value
    }
    return(y)
  }
  
  risk_baseline = integrate(integrand_gf,-Inf,Inf)$value
  return(risk_baseline)
}

risk_reduction_lowest_family_history = function(r2, h2, K, n, df, dm, n_samples = 10000)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail = F)
  
  integrand_lowest_given_parents = function(t,gm,gf)
  {
    arg = (zk - (gm+gf)/2 - t*r/sqrt(2)) / sqrt(1-h2/2-r2/2)
    y = log(n) + dnorm(t, log = T) + (n-1) * pnorm(t, lower.tail = F, log.p = T) + pnorm(arg, lower.tail = F, log.p = T)
    return(y)
  }
  
  posterior = function(gm,gf)
  {
    y = dnorm(gm, sd = h, log = T)
    y = y + dnorm(gf, sd = h, log = T)
    arg = (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y <- y + pnorm(arg, lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    arg = (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y <- y + pnorm(arg, lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    return(y)
  }
  
  integrand_gm = function(t, gms,gfs)
  {
    y <- integrand_lowest_given_parents(t, gms, gfs)
    post <- posterior(gms,gfs)
    y <- y + post
    return(y)
  }  
  
  opt <- optim(c(1, 1, 1), function(theta) -integrand_gm(theta[1], theta[2], theta[3]), hessian = T)
  
  # Sample points from t distribution
  data <- mvtnorm::rmvt(n_samples, sigma = solve(opt$hessian), df = 3, delta = opt$par)
  y <- exp(integrand_gm(data[, 1], data[, 2], data[, 3]))

  risk_selection = mean(y / mvtnorm::dmvt(data, solve(opt$hessian), df = 3, delta = opt$par, log = F))
  sd <- sqrt(mean((y / mvtnorm::dmvt(data, solve(opt$hessian), df = 3, delta = opt$par, log = F) - risk_selection)^2) / n_samples)
  risk_baseline <- baseline_risk(r2, h2, K, df, dm)
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction, sd))
}

risk_reduction_exclude_family_history <- function(r2, h2, K, q, n, df, dm, n_samples = 10000)
{
  # We need h > r, so if not, subtract epsilon from r. Unless h is zero, and then we add epsilon to h
  if (h2 == 0) {h2 <- h2 + 0.0001}
  else if (h2 == r2) {r2 <- r2 - 0.0001}
  
  r <- sqrt(r2)
  h <- sqrt(h2)
  zk <- qnorm(K, lower.tail = F)
  zq <- qnorm(q, lower.tail = F)
  
  posterior <- function(gm,gf)
  {
    y <- dnorm(gm, sd = h, log = T)
    y <- y + dnorm(gf, sd = h, log = T)
    arg <- (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y <- y + pnorm(arg,lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    arg <- (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y <- y + pnorm(arg,lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    return(y)
  }
  
  integrand_t <- function(t,gm,gf)
  {
    arg <- (zk-t*r/sqrt(2)-(gm+gf)/2)/sqrt(1-h2/2-r2/2)
    y <- dnorm(t, log = T) + pnorm(arg, lower.tail = F, log.p = T)
    return(y)
  }
  
  integrand_c <- function(cs, gm, gf, t)
  {
    gamma = zq*sqrt(2) - cs/(r/sqrt(2))
    denom <- pnorm(gamma)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    f <- integrand_t(t, gm, gf)
    
    numer <- log(1-(1-denom)^n)
    # print(sum(t<=gamma))
    y <- ifelse(t <= gamma, numer-log(denom), (n-1) * log(1-denom)) + f
    y <- y + dnorm(cs, mean=r2/h2 * (gm+gf)/2, sd=r/h * sqrt((h2-r2)/2), log = T)
    
    return(y)
  }
  
  integrand_gm = function(gms,gfs, c, t)
  {
    y <- integrand_c(c, gms, gfs, t)
    post <- posterior(gms, gfs)
    y <- y + post
    return(y)
  }
  # opt <- optim(c(1, 1, 1, zq*sqrt(2) - 1/(r/sqrt(2)) - 1), function(theta) -integrand_gm(theta[1], theta[2], theta[3], theta[4]), hessian = T)
  opt <- optim(c(1, 1, 1, 1), function(theta) -integrand_gm(theta[1], theta[2], theta[3], theta[4]), hessian = T)
  
  # print(opt)
  # print(solve(opt$hessian))
  
  # Sample points from t distribution
  data <- mvtnorm::rmvt(n_samples, sigma = solve(opt$hessian), df = 3, delta = opt$par)
  y <- exp(integrand_gm(data[, 1], data[, 2], data[, 3], data[, 4]))
  
  risk_selection = mean(y / mvtnorm::dmvt(data, solve(opt$hessian), df = 3, delta = opt$par, log = F))
  sd <- sqrt(mean((y / mvtnorm::dmvt(data, solve(opt$hessian), df = 3, delta = opt$par, log = F) - risk_selection)^2) / n_samples)
  risk_baseline <- baseline_risk(r2, h2, K, df, dm)
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction, sd))
}

simulate_lowest_risk_two_traits = function(r2A,r2B,rho,KA,KB,ns,nfam=10000)
{
  # results = array(0,c(length(ns),4))
  results <- matrix(nrow = length(ns), ncol = 4)
  
  disease_countA = numeric(length(ns))
  disease_count_randomA = numeric(length(ns))
  disease_countB = numeric(length(ns))
  disease_count_randomB = numeric(length(ns))
  
  baselineA = numeric(length(ns))
  baselineB = numeric(length(ns))
  
  t_sickA = qnorm(1-KA)
  t_sickB = qnorm(1-KB)
  
  rA = sqrt(r2A)
  rB = sqrt(r2B)
  
  mu = c(0,0)
  sigma = 1/2 * matrix(c(r2A,rho*rA*rB,rho*rA*rB,r2B),nrow=2)
  
  # chol_sigma <- chol(sigma)
  cs <- matrix(0, nrow = nfam, ncol = 2)
  
  for (i in seq_along(ns))
  {
    # cat('\r',i)
    
    n = ns[i]
    rmvn(nfam, mu, sigma, A = cs)
    
    xs <- rmvn(nfam*n, mu, sigma)
    xsA = matrix(xs[,1],nrow=nfam)
    xsB = matrix(xs[,2],nrow=nfam)
    
    envsA = rnorm(nfam,0,sqrt(1-r2A))
    envsB = rnorm(nfam,0,sqrt(1-r2B))
    
    scoresA <- xsA + cs[,1]#csA
    scoresB <- xsB + cs[,2]#csB
    
    selected_ind <- max.col(-scoresA)
    
    liabA <- scoresA[cbind(1:nfam, selected_ind)] + envsA
    liabB <- scoresB[cbind(1:nfam, selected_ind)] + envsB
    disease_countA <- sum(liabA > t_sickA)
    disease_countB <- sum(liabB > t_sickB)
    
    liab_randomA <- scoresA[, 1] + envsA
    disease_count_randomA <- sum(liab_randomA > t_sickA)
    liab_randomB <- scoresB[, 1] + envsB
    disease_count_randomB <- sum(liab_randomB > t_sickB)
    
    rrrA = (disease_count_randomA - disease_countA)/disease_count_randomA
    arrA = (disease_count_randomA - disease_countA)/nfam
    rrrB = (disease_count_randomB - disease_countB)/disease_count_randomB
    arrB = (disease_count_randomB - disease_countB)/nfam
    
    results[i,1] = rrrA
    results[i,2] = arrA
    results[i,3] = rrrB
    results[i,4] = arrB
  }  
  return(results)
}


