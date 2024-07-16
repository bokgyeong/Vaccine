# =============================================================================-
# functions ----
# =============================================================================-

## zero-inflated Poisson distribution ----
dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })


rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })



## zero-inflated negative binomial distribution ----
dZINB <- nimbleFunction(
  run = function(x = integer(), prob = double(), size = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    
    if (x != 0) { ## First handle non-zero data
      if (log) { ## return the log probability if log = TRUE
        return(dnbinom(x, prob = prob, size = size, log = TRUE) + log(1 - zeroProb))
      } else { ## or the probability if log = FALSE
        return((1 - zeroProb) * dnbinom(x, prob = prob, size = size, log = FALSE))
      }

    } else { ## From here down we know x is 0
      totalProbZero <- zeroProb + (1 - zeroProb) * dnbinom(0, prob = prob, size = size, log = FALSE)
      if (log) {
        return(log(totalProbZero))
      } else {
        return(totalProbZero)
      }
    }
  }
)


rZINB <- nimbleFunction(
  run = function(n = integer(), prob = double(), size = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) {
      return(0)
    } else {
      return(rnbinom(1, prob = prob, size = size))    
    }
  }
)





# =============================================================================-
# models ----
# =============================================================================-

## spatial zero-inflated Poisson-lognormal model ----
s_zi_pois_lognormal <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    y[i] ~ dZIP(mu[i], 1-pii[i])
    logit(pii[i]) <- XB1[i] + MG[i]
    log(mu[i]) ~ dnorm(XB2[i] + MD[i], sigma2)
  }
  
  # Constant
  XB1[1:n] <- X[1:n,1:p] %*% beta1[1:p]
  XB2[1:n] <- X[1:n,1:p] %*% beta2[1:p]
  MG[1:n] <- M[1:n,1:q] %*% (gamma[1:q] * Igamma[1:q])
  MD[1:n] <- M[1:n,1:q] %*% (delta[1:q] * Idelta[1:q])
  preMatG[1:q,1:q] <- kappa * Qs[1:q,1:q]
  preMatD[1:q,1:q] <- tau * Qs[1:q,1:q]
  
  # Process Model
  gamma[1:q] ~ dmnorm(mean = gamma0[1:q], prec = preMatG[1:q,1:q])
  delta[1:q] ~ dmnorm(mean = delta0[1:q], prec = preMatD[1:q,1:q])
  for(j in 1:q){
    Igamma[j] ~ dbern(phi)
    Idelta[j] ~ dbern(phi)
  }
  
  # Parameter Model
  for(j in 1:p){
    beta1[j] ~ dnorm(0, sd = sqrt(100))
    beta2[j] ~ dnorm(0, sd = sqrt(100))
  }
  kappa ~ dgamma(0.001,1000)
  tau ~ dgamma(0.001,1000)
  sigma2 ~ dinvgamma(2, 2)
})


## spatial zero-inflated Poisson regression model ----
s_zi_pois_reg <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    y[i] ~ dZIP(mu[i], 1-pii[i])
    logit(pii[i]) <- XB1[i] + MG[i]
    log(mu[i]) <- XB2[i] + MD[i]
  }
  
  # Constant
  XB1[1:n] <- X[1:n,1:p] %*% beta1[1:p]
  XB2[1:n] <- X[1:n,1:p] %*% beta2[1:p]
  MG[1:n] <- M[1:n,1:q] %*% (gamma[1:q] * Igamma[1:q])
  MD[1:n] <- M[1:n,1:q] %*% (delta[1:q] * Idelta[1:q])
  preMatG[1:q,1:q] <- kappa * Qs[1:q,1:q]
  preMatD[1:q,1:q] <- tau * Qs[1:q,1:q]
  
  # Process Model
  gamma[1:q] ~ dmnorm(mean = gamma0[1:q], prec = preMatG[1:q,1:q])
  delta[1:q] ~ dmnorm(mean = delta0[1:q], prec = preMatD[1:q,1:q])
  for(j in 1:q){
    Igamma[j] ~ dbern(phi)
    Idelta[j] ~ dbern(phi)
  }
  
  # Parameter Model
  for(j in 1:p){
    beta1[j] ~ dnorm(0, sd = sqrt(100))
    beta2[j] ~ dnorm(0, sd = sqrt(100))
  }
  kappa ~ dgamma(0.001,1000)
  tau ~ dgamma(0.001,1000)
})


## spatial zero-inflated negative binomial regression model ----
s_zi_nb_reg <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    y[i] ~ dZINB(prob = prob[i], size = theta, zeroProb = 1-pii[i])
    prob[i] <- 1 - mu[i] / (mu[i] + theta)
    logit(pii[i]) <- XB1[i] + MG[i]
    log(mu[i]) <- XB2[i] + MD[i]
  }
  
  # Constant
  XB1[1:n] <- X[1:n,1:p] %*% beta1[1:p]
  XB2[1:n] <- X[1:n,1:p] %*% beta2[1:p]
  MG[1:n] <- M[1:n,1:q] %*% (gamma[1:q] * Igamma[1:q])
  MD[1:n] <- M[1:n,1:q] %*% (delta[1:q] * Idelta[1:q])
  preMatG[1:q,1:q] <- kappa * Qs[1:q,1:q]
  preMatD[1:q,1:q] <- tau * Qs[1:q,1:q]
  
  # Process Model
  gamma[1:q] ~ dmnorm(mean = gamma0[1:q], prec = preMatG[1:q,1:q])
  delta[1:q] ~ dmnorm(mean = delta0[1:q], prec = preMatD[1:q,1:q])
  for(j in 1:q){
    Igamma[j] ~ dbern(phi)
    Idelta[j] ~ dbern(phi)
  }
  
  # Parameter Model
  for(j in 1:p){
    beta1[j] ~ dnorm(0, sd = sqrt(100))
    beta2[j] ~ dnorm(0, sd = sqrt(100))
  }
  kappa ~ dgamma(0.001,1000)
  tau ~ dgamma(0.001,1000)
  log(theta) ~ dnorm(0, sd = sqrt(100))
})


## nonspatial zero-inflated negative binomial regression model ----
zi_nb_reg <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    y[i] ~ dZINB(prob = prob[i], size = theta, zeroProb = 1-pii[i])
    prob[i] <- 1 - mu[i] / (mu[i] + theta)
    logit(pii[i]) <- XB1[i]
    log(mu[i]) <- XB2[i]
  }
  
  # Constant
  XB1[1:n] <- X[1:n,1:p] %*% beta1[1:p]
  XB2[1:n] <- X[1:n,1:p] %*% beta2[1:p]
  
  # Parameter Model
  for(j in 1:p){
    beta1[j] ~ dnorm(0, sd = sqrt(100))
    beta2[j] ~ dnorm(0, sd = sqrt(100))
  }
  log(theta) ~ dnorm(0, sd = sqrt(100))
})
