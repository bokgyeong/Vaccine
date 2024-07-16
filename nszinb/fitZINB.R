rm(list=ls())
require(nimble)
source('src/nimbleFtns.R')

# directory ----
dirn = 'nszinb'

filepath.fit = paste0(dirn, '/fit/')
ifelse(!dir.exists(filepath.fit), dir.create(filepath.fit, recursive = T), FALSE)


# load data ----
load(paste0(dirn, '/data/sim.RData'))


# preliminaries for nimble ----
p = ncol(X)
consts = list(n = nrow(X), p = p)
data = list(y = y, X = X)
inits = list(
  beta1 = rnorm(p),
  beta2 = rnorm(p),
  theta = 1
)


# run the MCMC algorithm ----
ptm = proc.time()[3]
mcmc.out = nimbleMCMC(
  code = zi_nb_reg, data = data, constants = consts, inits = inits,
  niter = 50000, nburnin = 1000, nchains = 1, WAIC = T, summary = T,
  monitors = c('beta1', 'beta2', 'theta')
)
rtime = proc.time()[3] - ptm

save(mcmc.out, rtime, file = paste0(filepath.fit, 'simZINB.RData'))
