rm(list = ls())
library(Rcpp); library(RcppArmadillo)
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach); library(doParallel); library(parallel)

source('src/RFtns.R')
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


load('data/refuse4y.RData')
load('fit/refuse4yZINB.RData')


# =============================================================================-
# Compute RQRs ----
# =============================================================================-
zicomp_refuse = data %>%
  dplyr::select(year, month, fips, y, cov, under5) %>%
  add_column(pii = piihat, mu = muhat, nu = nuhat, est = muhat * piihat,
             # rqr = stdRQR_modeZICOMP(data$y, modehat, nuhat, piihat, n = 100000),
             V = Vhat, W = What)

save(zicomp_refuse, file = paste0('realnu/summary/refuse6m.RData'))



# -----------------------------------------------------------------------------=
# Residual plots
# -----------------------------------------------------------------------------=
load(paste0('realnu/summary/refuse6m.RData'))

plot.rqr.scatter = zicomp_refuse %>%
  filter(rqr != Inf) %>%
  ggplot(aes(x = log(mu), y = rqr)) +
  geom_point(size = 0.3) +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = -3, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 3, color = 'red', linetype = 'dashed') +
  theme(strip.text = element_text(size = 13)) +
  labs(x = 'Log of mean estimate', y = 'Standardized RQR', title = '(a) Residual plot')

plot.rqr.qq = zicomp_refuse %>%
  filter(rqr != Inf) %>%
  ggplot(aes(sample = rqr)) +
  stat_qq(size = 0.3) +
  stat_qq_line(color = 'red') +
  theme(strip.text = element_text(size = 13)) +
  labs(x = 'Standard normal quantile', y = 'Standardized RQR', title = '(b) Quantile-quantile plot')

rqr.final = ggarrange(plot.rqr.scatter, plot.rqr.qq, nrow = 1)

ggsave(plot = rqr.final, width = 6, height = 3.2, file = 'realnu/figures/refuseRQR.eps')
