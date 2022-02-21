library(MASS)
library(mvnfast)

# Some silly micro-optimization instead of using one more general log(pnorm(b) - pnorm(a)), 
# should also be more accurate though.
log_dtruncnorm_left <- function(x, a) {
  dnorm(x, log = T) - pnorm(a, lower.tail = F, log.p = T)
}

log_dtruncnorm_right <- function(x, b) {
  dnorm(x, log = T) - pnorm(b, log.p = T)
}

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
  risk = integrate(integrand_lowest,-Inf,Inf)$value
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
  integrand_u = function(us)
  {
    y = numeric(length(us))
    beta_vec = zq*sqrt(2)-us
    denom = pnorm(beta_vec)
    dnorm_u <- dnorm(us)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    
    for (i in seq_along(us))
    {
      u = us[i]
      beta <- beta_vec[i]
      internal_int = integrate(integrand_t,-Inf,beta,u)$value
      numer = dnorm_u[i]*(1-(1-denom[i])^n) * internal_int
      term1 = numer/denom[i]
      
      internal_int = integrate(integrand_t,beta,Inf,u)$value
      term2 = dnorm_u[i]*(1-denom[i])^(n-1) * internal_int
      y[i] = term1 + term2
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf)$value
  reduction = (K-risk)/K
  # K-risk
  return(reduction)
}

risk_reduction_lowest_conditional = function(r2,K,n,qf,qm,relative=T,parental_avg_given=F)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
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
  risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
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
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  gamma = zq*sqrt(2) - c/(r/sqrt(2))
  
  integrand_t = function(t)
  {
    y = dnorm(t)*pnorm((zk-t*r/sqrt(2)-c)/sqrt(1-r2),lower.tail=F)
    return(y)
  }
  
  internal_int = integrate(integrand_t,-Inf,gamma)$value
  denom = pnorm(gamma)
  numer = (1-(1-denom)^n) * internal_int
  term1 = numer/denom
  
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
risk_reduction_lowest_family_history2 = function(r2,h2,K,n,df,dm, n_samples = 10000)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  
  integrand_lowest_given_parents = function(t,gm,gf)
  {
    # The sample is with respect to t ~ Normal(0, 1), so I've removed the dnorm term.
    arg = (zk - (gm+gf)/2 - t*r/sqrt(2)) / sqrt(1-h^2/2-r^2/2)
    y = n * pnorm(t,lower.tail=F)^(n-1) * pnorm(arg, lower.tail=F)
    return(y)
  }
  
  posterior = function(gm,gf)
  {
    # Again the the dnorm terms are removed + I've used the change of variable
    # gf' = gf/h
    # gm' = gm/h
    # so it will work with N(0, 1). Though it did force me to use h*gf' in every place.
    y <- 1
    arg = (zk-h*gm)/sqrt(1-h2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-h*gf)/sqrt(1-h2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
  
  integrand_gm = function(gms,gf,baseline, x2)
  {
    if (baseline) {
      arg = (zk - (h*gms+h*gf)/2) / sqrt(1-h^2/2)
      y <- pnorm(arg, lower.tail = F)
    }
    else {
      y <- (integrand_lowest_given_parents(x2, h*gms, h*gf))
    }
    post = posterior(gms,gf)
    y = y * post
    return(y)
  }
  
  integrand_gf = function(gfs,baseline)
  {
    # Probably needs some refactoring, basically we estimate the integral with control variate
    # by regression. .lm.fit is quite a bit faster than then normal lm.
    # Perhaps I should use some multivariate basis function instead of those?
    x <- rnorm(n_samples)
    x2 <- rnorm(n_samples)
    y <- integrand_gm(x, gfs, baseline, x2)
    mat <- cbind(1, x, tanh(x), gfs, tanh(gfs), x2, tanh(x2))
    fit <- .lm.fit(mat, y)
    std <- sqrt(sum((fit$resid)^2) / (length(x) - fit$rank))
    return(list(y = fit$coef[1], std = sqrt(diag(solve(crossprod(mat)))[1]) * std))
  }
  risk_selection = integrand_gf(rnorm(n_samples), F)
  risk_baseline = integrand_gf(rnorm(n_samples), T)

  relative_reduction = (risk_baseline$y-risk_selection$y)/risk_baseline$y
  abs_reduction = risk_baseline$y-risk_selection$y
  return(c(risk_baseline$y,risk_selection$y,relative_reduction,abs_reduction,
           risk_baseline$std, risk_selection$std))
}

risk_reduction_exclude_family_history2 = function(r2,h2,K,q,n,df,dm, n_samples = 10000)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  # Most of the same comments applies here.
  
  posterior = function(gm,gf)
  {
    y <- 1
    arg = (zk-h*gm)/sqrt(1-h2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-h*gf)/sqrt(1-h2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
  
  integrand_t = function(t,gm,gf)
  {
    arg = (zk-t*r/sqrt(2)-(gm+gf)/2)/sqrt(1-h2/2-r2/2)
    y = pnorm(arg,lower.tail=F, log.p = T)
    return(y)
  }
  
  integrand_c = function(cs,gm,gf)
  {
    # Again did change of variable so it will be N(0, 1).
    c <- cs
    modifier <- r/h * sqrt((h2-r2)/2)
    c = modifier * c + r2/h2 * (gm+gf)/2
    
    name <- (r/h * sqrt((h2-r2)/2)) * c + r2/h2 * (gm+gf)/2
    
    gamma = zq*sqrt(2) - name/(r/sqrt(2))
    denom = pnorm(gamma)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    
    samples <- qnorm(runif(length(cs), 0, denom))
    samples2 <- samples
    f1 <- integrand_t(samples, gm, gf) + dnorm(samples, log = T) - log_dtruncnorm_right(samples, b=gamma)
    internal_int <- exp(f1)
    
    numer = (1-(1-denom)^n) * internal_int
    term1 = numer/denom
    
    samples <- qnorm(runif(length(cs), denom, 1))
    f2 <- integrand_t(samples, gm, gf) + dnorm(samples, log = T) - log_dtruncnorm_left(samples, a=gamma)
    internal_int <- exp(f2)
    term2 = (1-denom)^(n-1) * internal_int
    y <- term1 + term2
    return(list(y=y, x1 = samples2 + dnorm(gamma) / denom, x2 = samples - dnorm(gamma) / (1-denom)))
  }
  
  integrand_gm = function(gms,gf,baseline)
  {
    if(baseline) {
      arg <- (zk - (h*gms+h*gf)/2) / sqrt(1-h2/2)
      y <- list(y=pnorm(arg, lower.tail = F), x1=NULL, x2=NULL)
    }
    else {
      y <- integrand_c(rnorm(length(gms)), h*gms, h*gf)
    }
    post = posterior(gms,gf)
    y <- list(y=y$y * post, x1 = y$x1, x2 = y$x2)
    return(y)
  }
  
  integrand_gf = function(gfs,baseline)
  {
    x3 <- rnorm(length(gfs))
    y <- c(integrand_gm(x3, gfs, baseline), x3 = list(x3))
    return(y)
  }
  
  # From few tests, it seems that other than x, the rest don't really reduce the variance.
  # I really need to refactor it since it's too different, coding wise, from the lowest_family_history.
  x <- rnorm(n_samples)
  risk_baseline = integrand_gf(x, T)
  x3 <- risk_baseline$x3
  risk_baseline <- risk_baseline$y
  
  mat <- cbind(1, x, tanh(x))
  
  fit <- .lm.fit(mat, risk_baseline)
  risk_baseline <- fit$coef[1]
  
  # Estimate the std
  std <- sqrt(sum((fit$resid)^2) / (length(x) - fit$rank))
  baseline_std <- sqrt(diag(solve(crossprod(mat)))[1]) * std 
  
  x <- rnorm(n_samples)
  risk_selection = integrand_gf(x, F)
  x1 <- risk_selection$x1
  x2 <- risk_selection$x2
  x3 <- risk_selection$x3
  risk_selection <- risk_selection$y
  
  mat <- cbind(1, x, tanh(x))
  
  fit <- .lm.fit(mat, risk_selection)

  risk_selection <- fit$coef[1]
  
  # Estimate the std
  std <- sqrt(sum((fit$resid)^2) / (length(x) - fit$rank))
  selection_std <- sqrt(diag(solve(crossprod(mat)))[1]) * std 
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction,
           baseline_std, selection_std))
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


