rm(list = ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
options(bitmapType='cairo')
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
everyfifth <- function(x){
  x[-seq(5, length(x), 5)] <- ""
  return(x)
}


# =============================================================================-
# Load data and results ----
# =============================================================================-
load("realnu/data/refuse6m.RData")

q = 100
load(paste0('realnu/fit/fitCOMP.RData'))

M = M0[,1:q]
Qs = t(M) %*% Q %*% M

beta1Ind = 1:p1
beta2Ind = (p1+1):(p1+p2)
alphaInd = p1+p2+1
zetaInd = (p1+p2+1+1):(p1+p2+1+q)
gammaInd = (p1+p2+1+q+1):(p1+p2+1+q+q)
IzetaInd = (p1+p2+1+q+q+1):(p1+p2+1+q+q+q)
IgammaInd = (p1+p2+1+q+q+q+1):(p1+p2+1+q+q+q+q)
omegaInd = p1+p2+1+q+q+q+q+1
kappaInd = p1+p2+1+q+q+q+q+1+1

ts.plot(postSamples[,omegaInd])
ts.plot(postSamples[,kappaInd])

which(postSamples[,omegaInd] < 1)[1]
which(postSamples[,kappaInd] < 1)[1]

burn = 20000
burn*thin

niter = burn + 20000
niter * thin

posterior = postSamples[1:niter,]

rtime / nrow(postSamples) * niter / 60 / 60 / 24


dummy = apply(postSamples[seq(updateInterval/thin, nrow(postSamples), updateInterval/thin),IzetaInd], 2, diff)
dim(dummy)
head(dummy)
dummy[dummy == -1] = 1
mean(apply(dummy, 2, mean))

dummy = apply(postSamples[seq(updateInterval/thin, nrow(postSamples), updateInterval/thin),IgammaInd], 2, diff)
dim(dummy)
head(dummy)
dummy[dummy == -1] = 1
mean(apply(dummy, 2, mean))


# =============================================================================-
# Summary ----
# =============================================================================-
bm_est = function(x){ bm(x)$est }
bm_se = function(x){ bm(x)$se }
hpd = function(x){ paste0('(', round(HPDinterval(as.mcmc(x))[1], 2), ',', round(HPDinterval(as.mcmc(x))[2], 2), ')') }
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }

parsInd = c(beta1Ind, beta2Ind, alphaInd, omegaInd, kappaInd)

df.summary = data.frame(
  mean = round(apply(posterior[-(1:burn),parsInd], 2, bm_est), 2),
  median = round(apply(posterior[-(1:burn),parsInd], 2, median), 2),
  hpd = apply(posterior[-(1:burn),parsInd], 2, hpd),
  hpd1 = apply(posterior[-(1:burn),parsInd], 2, hpd1),
  hpd2 = apply(posterior[-(1:burn),parsInd], 2, hpd2),
  mcse = round(apply(posterior[-(1:burn),parsInd], 2, bm_se), 4),
  ess = round(apply(posterior[-(1:burn),parsInd], 2, ess)),
  acc = round(c(rep(Accprob[1], p1), rep(Accprob[2], p2), Accprob[3], Accprob[6:7]), 2), 
  row.names = c(paste0('beta1', 1:p1), paste0('beta2', 1:p2), 'alpha', 'omega', 'kappa'))

df.summary.final = df.summary %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  dplyr::select(-c(hpd1, hpd2))

df.summary.final
df.summary.final %>% dplyr::select(c(mean, median, hpd, ess, sig))




# =============================================================================-
# Convergence check ----
# =============================================================================-

df.posterior = data.frame(Iteration = 1:nrow(posterior),
                          beta1 = unname(posterior[,beta1Ind]),
                          beta2 = unname(posterior[,beta2Ind]),
                          alpha = unname(posterior[,alphaInd]),
                          omega = unname(posterior[,omegaInd]),
                          kappa = unname(posterior[,kappaInd]))

df.hpd = data.frame(t(HPDinterval(as.mcmc(df.posterior[-(1:burn),-1]))))
df.ess = data.frame(t(round(apply(df.posterior[,-1], 2, ess))))
bm2 = function(x){ bm(x)$se }
df.mcse = data.frame(t(round(apply(df.posterior[,-1], 2, bm2), 4)))

df.acc = data.frame(t(round(c(rep(Accprob[1], p1), rep(Accprob[2], p2), Accprob[3], Accprob[6:7]), 2)))
colnames(df.acc) = colnames(df.mcse)

df.posterior = df.posterior %>% mutate(est.beta1.1 = cumsum(beta1.1)/seq_along(beta1.1),
                                       est.beta1.2 = cumsum(beta1.2)/seq_along(beta1.2),
                                       est.beta1.3 = cumsum(beta1.3)/seq_along(beta1.3),
                                       est.beta1.4 = cumsum(beta1.4)/seq_along(beta1.4),
                                       est.beta1.5 = cumsum(beta1.5)/seq_along(beta1.5),
                                       est.beta1.6 = cumsum(beta1.6)/seq_along(beta1.6),
                                       est.beta1.7 = cumsum(beta1.7)/seq_along(beta1.7),
                                       est.beta1.8 = cumsum(beta1.8)/seq_along(beta1.8),
                                       est.beta1.9 = cumsum(beta1.9)/seq_along(beta1.9),
                                       est.beta1.10 = cumsum(beta1.10)/seq_along(beta1.10),
                                       est.beta1.11 = cumsum(beta1.11)/seq_along(beta1.11),
                                       est.beta1.12 = cumsum(beta1.12)/seq_along(beta1.12),
                                       est.beta2.1 = cumsum(beta2.1)/seq_along(beta2.1),
                                       est.beta2.2 = cumsum(beta2.2)/seq_along(beta2.2),
                                       est.beta2.3 = cumsum(beta2.3)/seq_along(beta2.3),
                                       est.beta2.4 = cumsum(beta2.4)/seq_along(beta2.4),
                                       est.beta2.5 = cumsum(beta2.5)/seq_along(beta2.5),
                                       est.beta2.6 = cumsum(beta2.6)/seq_along(beta2.6),
                                       est.beta2.7 = cumsum(beta2.7)/seq_along(beta2.7),
                                       est.beta2.8 = cumsum(beta2.8)/seq_along(beta2.8),
                                       est.beta2.9 = cumsum(beta2.9)/seq_along(beta2.9),
                                       est.beta2.10 = cumsum(beta2.10)/seq_along(beta2.10),
                                       est.beta2.11 = cumsum(beta2.11)/seq_along(beta2.11),
                                       est.beta2.12 = cumsum(beta2.12)/seq_along(beta2.12),
                                       est.alpha = cumsum(alpha)/seq_along(alpha),
                                       est.omega = cumsum(omega)/seq_along(omega),
                                       est.kappa = cumsum(kappa)/seq_along(kappa))


# -----------------------------------------------------------------------------=
# beta1 - est ----
# -----------------------------------------------------------------------------=

est.beta1.1 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.1)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[11]), subtitle = paste0('ESS=', df.ess$beta1.1, ', MCSE=', df.mcse$beta1.1, ', ACC=', df.acc$beta1.1)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.2 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.2)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[12]), subtitle = paste0('ESS=', df.ess$beta1.2, ', MCSE=', df.mcse$beta1.2, ', ACC=', df.acc$beta1.2)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.3 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.3)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[13]), subtitle = paste0('ESS=', df.ess$beta1.3, ', MCSE=', df.mcse$beta1.3, ', ACC=', df.acc$beta1.3)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.4 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.4)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[14]), subtitle = paste0('ESS=', df.ess$beta1.4, ', MCSE=', df.mcse$beta1.4, ', ACC=', df.acc$beta1.4)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.5 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.5)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[15]), subtitle = paste0('ESS=', df.ess$beta1.5, ', MCSE=', df.mcse$beta1.5, ', ACC=', df.acc$beta1.5)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.6 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.6)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[16]), subtitle = paste0('ESS=', df.ess$beta1.6, ', MCSE=', df.mcse$beta1.6, ', ACC=', df.acc$beta1.6)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.7 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.7)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[17]), subtitle = paste0('ESS=', df.ess$beta1.7, ', MCSE=', df.mcse$beta1.7, ', ACC=', df.acc$beta1.7)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.8 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.8)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[18]), subtitle = paste0('ESS=', df.ess$beta1.8, ', MCSE=', df.mcse$beta1.8, ', ACC=', df.acc$beta1.8)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.9 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.9)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[19]), subtitle = paste0('ESS=', df.ess$beta1.9, ', MCSE=', df.mcse$beta1.9, ', ACC=', df.acc$beta1.9)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.10 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.10)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[110]), subtitle = paste0('ESS=', df.ess$beta1.10, ', MCSE=', df.mcse$beta1.10, ', ACC=', df.acc$beta1.10)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.11 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.11)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[111]), subtitle = paste0('ESS=', df.ess$beta1.11, ', MCSE=', df.mcse$beta1.11, ', ACC=', df.acc$beta1.11)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta1.12 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta1.12)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[112]), subtitle = paste0('ESS=', df.ess$beta1.12, ', MCSE=', df.mcse$beta1.12, ', ACC=', df.acc$beta1.12)) +
  theme(plot.subtitle = element_text(size = 7.4))


plot.est.beta1.1 = ggarrange(est.beta1.1, est.beta1.2, est.beta1.3, 
                             est.beta1.4, est.beta1.5, est.beta1.6,
                             ncol = 3, byrow = T)

plot.est.beta1.2 = ggarrange(est.beta1.7, est.beta1.8, est.beta1.9, 
                             est.beta1.10, est.beta1.11, est.beta1.12,
                             ncol = 3, byrow = T)



# -----------------------------------------------------------------------------=
# beta1 - hist ----
# -----------------------------------------------------------------------------=

hist.beta1.1 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.1)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.1[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.1[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[11]), y = 'Density')

hist.beta1.2 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.2)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.2[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.2[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[12]), y = 'Density')

hist.beta1.3 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.3)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.3[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.3[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[13]), y = 'Density')

hist.beta1.4 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.4)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.4[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.4[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[14]), y = 'Density')

hist.beta1.5 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.5)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.5[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.5[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[15]), y = 'Density')

hist.beta1.6 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.6)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.6[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.6[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[16]), y = 'Density')

hist.beta1.7 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.7)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.7[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.7[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[17]), y = 'Density')

hist.beta1.8 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.8)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.8[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.8[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[18]), y = 'Density')

hist.beta1.9 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.9)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.9[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.9[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[19]), y = 'Density')

hist.beta1.10 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.10)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.10[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.10[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[110]), y = 'Density')

hist.beta1.11 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.11)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.11[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.11[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[111]), y = 'Density')

hist.beta1.12 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta1.12)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta1.12[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta1.12[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[112]), y = 'Density')


plot.hist.beta1.1 = ggarrange(hist.beta1.1, hist.beta1.2, hist.beta1.3, 
                              hist.beta1.4, hist.beta1.5, hist.beta1.6,
                              ncol = 3, byrow = T)

plot.hist.beta1.2 = ggarrange(hist.beta1.7, hist.beta1.8, hist.beta1.9, 
                              hist.beta1.10, hist.beta1.11, hist.beta1.12,
                              ncol = 3, byrow = T)


# ggsave(plot = est.hist.beta.1, width = 8, height = 3,
#        filename = 'sim_month2/figures/full2_beta1.eps')




# -----------------------------------------------------------------------------=
# beta2 - est ----
# -----------------------------------------------------------------------------=

est.beta2.1 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.1)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[21]), subtitle = paste0('ESS=', df.ess$beta2.1, ', MCSE=', df.mcse$beta2.1, ', ACC=', df.acc$beta2.1)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.2 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.2)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[22]), subtitle = paste0('ESS=', df.ess$beta2.2, ', MCSE=', df.mcse$beta2.2, ', ACC=', df.acc$beta2.2)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.3 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.3)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[23]), subtitle = paste0('ESS=', df.ess$beta2.3, ', MCSE=', df.mcse$beta2.3, ', ACC=', df.acc$beta2.3)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.4 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.4)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[24]), subtitle = paste0('ESS=', df.ess$beta2.4, ', MCSE=', df.mcse$beta2.4, ', ACC=', df.acc$beta2.4)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.5 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.5)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[25]), subtitle = paste0('ESS=', df.ess$beta2.5, ', MCSE=', df.mcse$beta2.5, ', ACC=', df.acc$beta2.5)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.6 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.6)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[26]), subtitle = paste0('ESS=', df.ess$beta2.6, ', MCSE=', df.mcse$beta2.6, ', ACC=', df.acc$beta2.6)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.7 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.7)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[27]), subtitle = paste0('ESS=', df.ess$beta2.7, ', MCSE=', df.mcse$beta2.7, ', ACC=', df.acc$beta2.7)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.8 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.8)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[28]), subtitle = paste0('ESS=', df.ess$beta2.8, ', MCSE=', df.mcse$beta2.8, ', ACC=', df.acc$beta2.8)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.9 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.9)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[29]), subtitle = paste0('ESS=', df.ess$beta2.9, ', MCSE=', df.mcse$beta2.9, ', ACC=', df.acc$beta2.9)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.10 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.10)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[210]), subtitle = paste0('ESS=', df.ess$beta2.10, ', MCSE=', df.mcse$beta2.10, ', ACC=', df.acc$beta2.10)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.11 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.11)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[211]), subtitle = paste0('ESS=', df.ess$beta2.11, ', MCSE=', df.mcse$beta2.11, ', ACC=', df.acc$beta2.11)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.beta2.12 = df.posterior %>% ggplot(aes(x = Iteration, y = est.beta2.12)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*beta[212]), subtitle = paste0('ESS=', df.ess$beta2.12, ', MCSE=', df.mcse$beta2.12, ', ACC=', df.acc$beta2.12)) +
  theme(plot.subtitle = element_text(size = 7.4))


plot.est.beta2.1 = ggarrange(est.beta2.1, est.beta2.2, est.beta2.3, 
                             est.beta2.4, est.beta2.5, est.beta2.6,
                             ncol = 3, byrow = T)

plot.est.beta2.2 = ggarrange(est.beta2.7, est.beta2.8, est.beta2.9, 
                             est.beta2.10, est.beta2.11, est.beta2.12,
                             ncol = 3, byrow = T)



# -----------------------------------------------------------------------------=
# beta2 - hist ----
# -----------------------------------------------------------------------------=

hist.beta2.1 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.1)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.1[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.1[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[21]), y = 'Density')

hist.beta2.2 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.2)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.2[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.2[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[22]), y = 'Density')

hist.beta2.3 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.3)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.3[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.3[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[23]), y = 'Density')

hist.beta2.4 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.4)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.4[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.4[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[24]), y = 'Density')

hist.beta2.5 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.5)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.5[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.5[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[25]), y = 'Density')

hist.beta2.6 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.6)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.6[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.6[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[26]), y = 'Density')

hist.beta2.7 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.7)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.7[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.7[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[27]), y = 'Density')

hist.beta2.8 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.8)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.8[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.8[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[28]), y = 'Density')

hist.beta2.9 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.9)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.9[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.9[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[29]), y = 'Density')

hist.beta2.10 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.10)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.10[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.10[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[210]), y = 'Density')

hist.beta2.11 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.11)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.11[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.11[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[211]), y = 'Density')

hist.beta2.12 = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = beta2.12)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$beta2.12[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$beta2.12[2], color = 'blue', linetype = 2) +
  labs(x = expression(beta[212]), y = 'Density')


plot.hist.beta2.1 = ggarrange(hist.beta2.1, hist.beta2.2, hist.beta2.3, 
                              hist.beta2.4, hist.beta2.5, hist.beta2.6,
                              ncol = 3, byrow = T)

plot.hist.beta2.2 = ggarrange(hist.beta2.7, hist.beta2.8, hist.beta2.9, 
                              hist.beta2.10, hist.beta2.11, hist.beta2.12,
                              ncol = 3, byrow = T)





# -----------------------------------------------------------------------------=
# kappa, omega, alpha ----
# -----------------------------------------------------------------------------=

est.omega = df.posterior %>% ggplot(aes(x = Iteration, y = est.omega)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*omega), subtitle = paste0('ESS=', df.ess$omega, ', MCSE=', df.mcse$omega, ', ACC=', df.acc$omega)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.kappa = df.posterior %>% ggplot(aes(x = Iteration, y = est.kappa)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*kappa), subtitle = paste0('ESS=', df.ess$kappa, ', MCSE=', df.mcse$kappa, ', ACC=', df.acc$kappa)) +
  theme(plot.subtitle = element_text(size = 7.4))

est.alpha = df.posterior %>% ggplot(aes(x = Iteration, y = est.alpha)) +
  geom_line() +
  labs(x = 'Sample size', y = expression('Estimate of '*alpha), subtitle = paste0('ESS=', df.ess$alpha, ', MCSE=', df.mcse$alpha, ', ACC=', df.acc$alpha)) +
  theme(plot.subtitle = element_text(size = 7.4))

hist.omega = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = omega)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$omega[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$omega[2], color = 'blue', linetype = 2) +
  labs(x = expression(omega), y = 'Density')

hist.kappa = df.posterior %>%
  slice((burn+1):n()) %>% 
  ggplot(aes(x = kappa)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$kappa[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$kappa[2], color = 'blue', linetype = 2) +
  labs(x = expression(kappa), y = 'Density')

hist.alpha = df.posterior %>% 
  slice((burn+1):n()) %>% 
  ggplot(aes(x = alpha)) +
  geom_histogram(aes(y = ..density..),color = "black", fill = "white") +
  geom_vline(xintercept = df.hpd$alpha[1], color = 'blue', linetype = 2) +
  geom_vline(xintercept = df.hpd$alpha[2], color = 'blue', linetype = 2) +
  labs(x = expression(alpha), y = 'Density')


plot.est.hist.kap.omega.alpha = ggarrange(est.omega, est.kappa, est.alpha,
                                          hist.omega, hist.kappa, hist.alpha,
                                          ncol = 3, byrow = T)




# =============================================================================-
# Barplots for basis vectors ----
# =============================================================================-


posterior2 = posterior[-(1:burn),]

Izetahat = apply(posterior2[,IzetaInd], 2, mean)
Igammahat = apply(posterior2[,IgammaInd], 2, mean)

dat_ind = data.frame(Basis = as.factor(1:q), 
                     Name = '(a) Detection', Value = Izetahat)
dat_ind = rbind(dat_ind, data.frame(Basis = as.factor(1:q), 
                                    Name = '(b) Refusal cases', Value = Igammahat))

everynth <- function(x, n){
  x[-seq(n, length(x), n)] <- ""
  return(x)
}

plot_ind = dat_ind %>% ggplot(aes(Basis, Value)) +
  geom_bar(stat = 'identity', color = "black", fill = "gray") +
  coord_cartesian(ylim = c(0, 1)) +
  # facet_wrap(~ Name, ncol = 1) +
  facet_wrap(~ Name, nrow = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45), label = everynth(1:q, 5)) +
  labs(x = 'Basis vector', y = 'Posterior probability of inclusion')
plot_ind


dat_ind %>% 
  filter(Name == '(a) Detection') %>% 
  # filter(Name == '(b) Refusal cases') %>% 
  select(Value) %>% 
  filter(Value > 0.5)

ggsave(plot = plot_ind, width = 6, height = 3.5, file = 'realnu/fig/refusebasis.pdf')

