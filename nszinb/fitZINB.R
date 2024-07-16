rm(list=ls())
require(nimble)

source('src/nimbleFtns.R')

# runID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# ==============================================================================
# load data ----
# ==============================================================================
dataname = 'refuse4y'

load(paste0('data/', dataname, '.RData'))



# ==============================================================================
# Fit the model ----
# ==============================================================================

filename = paste0('fit/', dataname, 'ZINB.RData')


## preliminaries for nimble ----
p = ncol(X)
consts = list(n = nrow(X), p = p)
data = list(y = y, X = X)
inits = list(
  beta1 = rnorm(p),
  beta2 = rnorm(p),
  theta = 1
)




### run the MCMC algorithm ----
ptm = proc.time()[3]
mcmc.out = nimbleMCMC(code = zi_nb_reg, data = data, constants = consts, inits = inits,
                      niter = 50000, nburnin = 1000, nchains = 1, WAIC = T, summary = T,
                      # niter = 2000, nburnin = 1000, nchains = 1, WAIC = T, summary = T,
                      monitors = c('beta1', 'beta2', 'theta'))
rtime = proc.time()[3] - ptm

# mcmc.out$samples
# mcmc.out$summary
# mcmc.out$WAIC

save(mcmc.out, file = filename)
