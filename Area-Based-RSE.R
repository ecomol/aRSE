################################################################################################
################################################################################################
#R code for estimating minimal habitat area size 
#for detecting previously-unsampled rare species
#Written by Tsung-Jen Shen and Youhua Chen
#2023-10-30
#Updated: 2024-03-15
################################################################################################
################################################################################################

library(nlme)

## Convert "frequency counts data" to "species counts data" or "abundance data"
f.to.X = function(f) {
  X = c()
  k = 1
  while (k <= length(f)) {
    if (k > 1)
      X = c(X, rep(k, f[k]))
    else {
      X = rep(k, f[k])
    }
    k = k + 1
    
  }
  return(X)
  
}

## Convert "species counts data" to "frequency counts data"
X.to.f = function(X) {
  f = factor(X, levels = 0:max(X))
  f = table(f, exclude = 0) ## frequency counts
  dimnames(f) = NULL
  return(f)
}

#For area-based data
### Poisson model
PoiLam <- function(x, zero = FALSE) {
  x <- unlist(x)
  # print(x)
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (zero == FALSE)
    x <- x[x > 0]
  b.solve <- function(b) {# Find the solution of \beta from Eq. (3a) and (3b) in the main text
    a <- f1 / sum(x * exp(-b * x))# This is an equivalent expression of Eq. (3a)
    obj <-
      sum((x * (1 - a * exp(-b * x))) ^ 2) - sum(x * (x - 1)) + 2 *
      f2
    abs(obj)
  }
  b <-
    tryCatch(
      optimize(b.solve, c(0, 1))$min,
      error = function(a) {
        1
      }
    )
  a <- f1 / sum(x * exp(-b * x))
  
  return(c(a, b))
}


## Chao1 estimator for product Poissons
Chao1.Pois = function(f) {
  est = sum(f)
  if (length(f) > 1) {
    if (f[2] > 0) {
      A = f[1] / f[2]
    } else {
      if (f[1] > 0) {
        A = f[1] - 1
      } else {
        A = 0
      }
    }
  } else
    A = f[1] - 1
  est = est + A * f[1] / 2
  
  return(est)
}


# f is frequency counts data
# a is the surveyed areal size
# h is the areal size of additional sample
# b is the estimated parameters (a, b) in dual equations for estimating true detection intensities of species
# an example: BW.Rh(f = c(10,5,3,0,4), h = 2, b = c(.1,.2))
BW.Rh = function(f, a = 1, h, b, k.show = 5) {
  h = h / a
  
  kmax = length(f)
  cut.pts = k.show
  
  ## lambda i
  d = rep(0, cut.pts)
  ## estimated fk
  est.Rh = rep(0, cut.pts)
  
  ### detection intensities estimation
  ind = 1:cut.pts
  d = ind * (1 - b[1] * exp(-b[2] * ind))
  
  for (i in 1:cut.pts) {
    ## i for tau
    for (j in 1:min(i, kmax)) {
      ## j for q
      est.Rh[i] = est.Rh[i] + 1 / ((1 + 1 / h) ^ i - 1) *
        f[j] * dpois(i - j, h * d[j])
      
    }
  }
  
  ## return a vector of (Rh_[t=1],Rh_[t=2],...)
  return(est.Rh)
}

##### Unweighted estimator
UW.Rh = function(f = NULL,
                 xi = NULL,
                 a = 1,
                 h,
                 b,
                 f0,
                 k.show = 3) {
  h = h / a
  
  if (is.null(f) * is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL)
    
  }
  if (is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if (is.null(xi)) {
    Xi = f.to.X(f)
  }
  
  kmax = length(f)
  
  cut.pts = k.show
  
  est.Rh = rep(0, cut.pts)
  if (length(f) == 1 && f[1] > 0) {
    f = c(f[1], 0)
  }
  if (f[2] == 0) {
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  
  d0 = 0
  if (f0 > 0)
    d0 = f1 / f0
  
  ### detection intensities estimation
  d = Xi * (1 - b[1] * exp(-b[2] * Xi))
  
  for (k in 1:cut.pts) {
    est.Rh[k] = sum(exp(-d) * dpois(k, h * d))
    
    est.Rh[k] = est.Rh[k] + exp(-d0) * f0 * dpois(k, h * d0)
  }# loop: k
  
  return(est.Rh)
}



## Inferring the minimal area size (denoted by h) when the predetermined number of rare species (denoted by sh) is given
## sh (threshold): the predetermined number of rare species
## tt: the maximal abundance of rare species which are asked to be conserved
Est.h = function(f,
                 a = 1,
                 b,
                 tt = 3,
                 sh = 1) {
  ## if f[2] = 0 then f[1] = f[1] - 1; f[2] = 1 
  if(length(f) == 1) f =  c(f[1]-1, 1)  
  ## main constraint
  obj.fun = function(x) {
    val = sum(BW.Rh(
      f = f,
      a = a,
      h = x,
      b = b,
      k.show = tt
    )) - sh
  }
  
  ## find the point (peak.pts) which attains the maximum value of "obj.fun"
  peak.pts <-
    optim(0.1, function(x) {
      -obj.fun(x)
    }, lower = 0.00001, upper = 5000, method = c("L-BFGS-B"))$par
  
  ## Check if obj.fun(peak.pts) is smaller than sh
  AA = try(uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15), silent = T)
  
  if (inherits(AA, "try-error")) {
    ## "sh" is larger than "sum(BW.Rh())"
    h1 = h2 = peak.pts
  } else {
    h2 <- uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15)$root
    # print(h2)
    h1 <- uniroot(obj.fun, c(0.00001, peak.pts), tol = 1e-15)$root
  }
  # print(list(f=f,h=c(h1,h2),AA=peak.pts,tt=tt,sh=sh,a=a))
  c(h1, h2)
}



## Calculate the true value of h whose numerical solution can be solved out by
## the equation 7a or 7b in the main text
Real.h = function(lambda,
                  a = 1,
                  tt = 3,
                  sh = 1) {
  ## main constraint
  obj.fun = function(h) {
    val = 0
    for (i in 1:tt) {
      val = val + sum(exp(-a * lambda) * dpois(i, lambda * h))
    }
    val - sh
  }

  peak.pts <-
    optim(0.1, function(x) {
      -obj.fun(x)
    }, lower = 0.00001, upper = 5000, method = c("L-BFGS-B"))$par
  
  
  AA = try(uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15), silent = T)
  
  
  if (inherits(AA, "try-error")) {
    ## "sh"(or "g" used in the main text) is larger than "sum(BW.Rh())"
    h1 = h2 = peak.pts
  } else {
    h2 <- uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15)$root
    h1 <- uniroot(obj.fun, c(0.00001, peak.pts), tol = 1e-15)$root
  }
  c(h1, h2)
}


## k is the clumping parameter of the NBD
## u stands for the expectation of the NBD
NBD.Real.h = function(k,
                      u,
                      a = 1,
                      tt = 3,
                      sh = 1) {
  ## main constraint
  obj.fun = function(h) {
    val = 0
    for (i in 1:tt) {
      val = val + sum(dnbinom(0, size = k, mu = a * u) * dnbinom(i, size = k, mu = h *
                                                                   u))
    }
    val - sh
  }

  peak.pts <-
    optim(0.1, function(x) {
      -obj.fun(x)
    }, lower = 0.00001, upper = 5000, method = c("L-BFGS-B"))$par
  
  AA = try(uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15), silent = T)
  
  if (inherits(AA, "try-error")) {
    ## "sh"(or the notation "g" used in the main text) is larger than "sum(BW.Rh())"
    h1 = h2 = peak.pts
  } else {
    h2 <- uniroot(obj.fun, c(peak.pts, 5000), tol = 1e-15)$root
    h1 <- uniroot(obj.fun, c(0.00001, peak.pts), tol = 1e-15)$root
  }
  c(h1, h2)
}


### Generating species frequency counts (fk's) using the bootstrapping method
# f is frequency counts
# b is the estimated parameters (a, b) in dual equations for estimating true detection intensities of species
boot.abundance.Pois.fun = function(S.hat, f, b) {
  D = sum(f)
  ind = 1:length(f)
  kmax = length(f)
  d = rep(0, kmax)
  n = sum(ind * f)
  S.hat = ceiling(S.hat)
  
  ## Estimated lambda of a species with Xi > 0
  boot.lambda = ind * (1 - b[1] * exp(-b[2] * ind))
  ## Generate a bootstrap sample
  boot.Xi = unlist(sapply(1:kmax,
                          function(x) {
                            rpois(f[x], boot.lambda[x])
                          }))
  
  
  if (S.hat > D) {
    ## Number of unseen species in the local sample
    f0 = S.hat - D
    ## Estimated lambda of a species with Xi = 0
    g0 = f[1] / f0
    
    boot.Xi = c(rpois(f0, g0), boot.Xi)
  }
  
  ## frequency counts
  f.count = X.to.f(boot.Xi)
  
  if(length(f.count) == 1) f.count =  c(f.count[1]-1, 1)
  
  return(f.count)
}



##### main function for estimating number of rare species with a specific abundance in an additional sample of area size 'h'
## f: species frequency counts data
## xi: species abundance data
Pred.abundance.Pois.rare = function(boot.rep = 50,
                                    a = 1,
                                    f = NULL,
                                    xi = NULL,
                                    h,
                                    k.show = 3,
                                    b.seed = 1234) {
  if (is.null(f) * is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL)
    
  }
  if (is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if (is.null(xi)) {
    Xi = f.to.X(f)
  }
  # #### Point estimates
  a.p.1 = PoiLam(Xi, zero = FALSE)
  est.R = BW.Rh(
    f = f,
    a = a,
    h = h,
    b = a.p.1,
    k.show = k.show
  )
  
  est.f0 = Chao1.Pois(f) - sum(f)
  est.R2 = UW.Rh(
    f = f,
    a = a,
    h = h,
    b = a.p.1,
    f0 = est.f0,
    k.show = k.show
  )
  
  #### Calculating bootstrap SEs and CIs
  #### Create a space for storing bootstrap samples
  boot.output = NULL
  set.seed(b.seed)
  
  for (i in 1:boot.rep) {
    b.f = boot.abundance.Pois.fun(S.hat = est.f0 + sum(f),
                                  f = f,
                                  b = a.p.1)
    
    b.Xi = f.to.X(b.f)
    b.a.p.1 = PoiLam(b.Xi, zero = FALSE)
    b.est.f0 = Chao1.Pois(b.f) - sum(b.f)
    
    b.pred.fk.BW = BW.Rh(
      f = b.f,
      a = a,
      h = h,
      b = b.a.p.1,
      k.show = k.show
    )
    
    b.pred.fk.unweighted = UW.Rh(
      f = b.f,
      a = a,
      h = h,
      b = b.a.p.1,
      f0 = b.est.f0,
      k.show = k.show
    )
    boot.output = rbind(boot.output,
                        c(proposed = b.pred.fk.BW[1:k.show], unweighted = b.pred.fk.unweighted[1:k.show]))
  }### loop: boot.rep
  
  point.est = cbind(proposed = est.R, unweighted = est.R2)
  boot.sd = apply(boot.output, 2, sd, na.rm = T)
  boot.sd = matrix(boot.sd, ncol = 2, byrow = F)
  boot.ci = apply(boot.output, 2, quantile, probs = c(0.025, 0.975))
  
  Normal.ci = cbind(point.est - 1.96 * boot.sd, point.est + 1.96 * boot.sd)
  for (i in 1:nrow(Normal.ci)) {
    for (j in 1:ncol(Normal.ci)) {
      if (Normal.ci[i, j] < 0)
        Normal.ci[i, j] = 0
    }
  }

  output = list()
  output[["Data information"]] = as.matrix(c(a, h),ncol = 1)

  rownames(output[["Data information"]]) = c(
    "  Area size of the original sample (a):    ",
    "  Area size of an additional sample (h):   "  )
  colnames(output[["Data information"]]) = c(" ")
  
  output[["Bayesian-weight method"]] = round(cbind(1:k.show, point.est[, 1], boot.sd[, 1], Normal.ci[, c(1, 3)]), 1)
  output[["Unweighted method"]] = round(cbind(1:k.show, point.est[, 2], boot.sd[, 2], Normal.ci[, c(2, 4)]), 1)
  
  colnames(output[["Bayesian-weight method"]]) =
    colnames(output[["Unweighted method"]]) =
    c("  k",
      "  Estimate",
      "  Estimated SE",
      "  95% lower limit",
      "  95% upper limit")
  
  output
} ### end of R function "Pred.abundance.Pois.rare"


Pred.Additional.area.Pois = function(boot.rep = 50,
                                     a = 1,
                                     f = NULL,
                                     xi = NULL,
                                     tt = 3,
                                     sh = 1,
                                     b.seed = 1234) {
  if (is.null(f) * is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL)
    
  }
  if (is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if (is.null(xi)) {
    Xi = f.to.X(f)
  }
  
  # #### Point estimates
  a.p.1 = PoiLam(Xi, zero = FALSE)
  est.h = Est.h(
    f = f,
    a = a,
    b = a.p.1,
    tt = tt,
    sh = sh
  )
  est.f0 = Chao1.Pois(f) - sum(f)

  
  #### Calculating bootstrap SEs and CIs
  #### Create a space for storing bootstrap samples
  boot.output = NULL
  set.seed(b.seed)
  
  for (i in 1:boot.rep) {
    b.f = boot.abundance.Pois.fun(S.hat = est.f0 + sum(f),
                                  f = f,
                                  b = a.p.1)
    b.Xi = f.to.X(b.f)
    
    ## Parameters estimated by using the bootstrapping sample
    b.a.p.1 = PoiLam(b.Xi, zero = FALSE)
    b.est.f0 = Chao1.Pois(b.f) - sum(b.f)
    b.est.h = min(Est.h(
      f = b.f,
      a = a,
      b = b.a.p.1,
      tt = tt,
      sh = sh
    ))

    boot.output = rbind(boot.output, b.est.h)
  }### loop: boot.rep
  

  point.est = c(proposed = min(est.h))
  boot.sd = apply(boot.output, 2, sd, na.rm = T)
  boot.sd = matrix(boot.sd, ncol = 1, byrow = F)
  boot.ci = apply(boot.output, 2, quantile, probs = c(0.025, 0.975))

  Normal.ci = cbind(point.est - 1.96 * boot.sd, point.est + 1.96 * boot.sd)
  for (i in 1:nrow(Normal.ci)) {
    for (j in 1:ncol(Normal.ci)) {
      if (Normal.ci[i, j] < 0)
        Normal.ci[i, j] = 0
    }
  }
  
  #print(list(point.est,boot.sd,Normal.ci))  
  
  output = list()
  output[["Data information"]] = as.matrix(c(a, tt, sh), ncol = 1)
  
  rownames(output[["Data information"]]) = c(
    "  Area size of the original sample (a):   ",
    "  Predetermined maximum abundance of rare species (\tau):   ",
    "  Predetermined number of rare species (g):   "
  )
  colnames(output[["Data information"]]) = c(" ")
  output[["Bayesian-weight method"]] = round(c(point.est, boot.sd[, 1], Normal.ci[, c(1, 2)]), 4)
  
  names(output[["Bayesian-weight method"]]) =
    c("  Estimated 'h'",
      "  Estimated SE",
      "  95% lower limit",
      "  95% upper limit")
  
  output
} ### end of R function "Pred.Additional.area.Pois"


## Bootstrapping SE of an estimated h using the Bayesian-weight estimator
BSE.BW.h = function(boot.rep = 50, a = 1, S.hat, f = NULL,
                    xi = NULL, b, tt, sh) {
  if (is.null(f) * is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL)
  }
  if (is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if (is.null(xi)) {
    Xi = f.to.X(f)
  }
  ## Chao1 estimator
  est.f0 = Chao1.Pois(f) - sum(f)
  
  ##########################
  boot.output = NULL
  
  for (i in 1:boot.rep) {
    b.f = boot.abundance.Pois.fun(S.hat = est.f0 + sum(f),
                                  f = f,
                                  b = a.p.1)
    b.Xi = f.to.X(b.f)
    
    ## Parameters estimated by using the bootstrapping sample
    b.a.p.1 = PoiLam(b.Xi, zero = FALSE)
    b.est.f0 = Chao1.Pois(b.f) - sum(b.f)
    # b.pp = rep(NA, k.show)
    
    b.pp = min(Est.h(
      f = b.f,
      a = a,
      b = b.a.p.1,
      tt = tt,
      sh = sh
    ))
    
    boot.output = c(boot.output, b.pp)
  }### loop: boot.rep
  
  pp.sd = sd(boot.output, na.rm = T)
  ##########################
  return(pp.sd)
}

## Fit NBD based on a conditional likelihood function
ku.cond.logf <- function(par, f, a = 1) {
  k <- par[1]
  u <- par[2]
  pp = k/(k + u*a)
  
  zz = which(f > 0)  
  rhoN = lgamma(zz + k) - lgamma(zz + 1) - lgamma(k) + zz * 
    log1p(-pp) + k * log1p(pp - 1) - log(1 - pp^k)
  
  ### log likelihood
  res <- -sum(f[zz] * rhoN)
  
  res
}



