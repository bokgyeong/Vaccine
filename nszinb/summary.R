rm(list = ls())
library(coda); library(tidyverse); library(egg); library(grid); library(patchwork)
library(batchmeans); library(foreach); library(doParallel); library(parallel)
library(ggbreak)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# ==============================================================================
# load results ----
# ==============================================================================

load('data/refuse4y.RData')
load('fit/refuse4yZINB.RData')

summat = mcmc.out$summary

beta1Ind = grep('beta1', rownames(summat))
beta2Ind = grep('beta2', rownames(summat))
thetaInd = grep('theta', rownames(summat))

beta1hat = summat[beta1Ind, 1]
beta2hat = summat[beta2Ind, 1]
thetahat = summat[thetaInd, 1]

etahat = as.vector(X %*% beta1hat )
piihat = exp(etahat)/(1 + exp(etahat))
muhat =  exp( as.vector(X %*% beta2hat) )



# ==============================================================================
# mean and HPD for regression coefficients ----
# ==============================================================================

Predictors = c('Intercept', 
               'log(Interaction)',
               'Health insurance',
               'Pediatrician reporting', 'Household size',
               'Religious congregation', 'Limited English',
               'Private school', 'High income', 'Same area', 
               'State law leniency', 'State autism')


summaryCoef = data.frame(
  Process = 'Detection', 
  Median = summat[beta1Ind,2], 
  CIlow = summat[beta1Ind,4],
  CIup = summat[beta1Ind,5],
  Predictor = Predictors)

summaryCoef = rbind(
  summaryCoef,
  data.frame(
    Process = 'Refusal cases', 
    Median = summat[beta2Ind,2], 
    CIlow = summat[beta2Ind,4],
    CIup = summat[beta2Ind,5],
    Predictor = Predictors))

summaryCoef$Predictor = factor(summaryCoef$Predictor, levels = rev(unique(summaryCoef$Predictor)))
summaryCoef = summaryCoef %>% mutate(hasZero = ifelse(CIlow <= 0 & CIup >= 0, 'Include zero', 'Do not include zero'))

pBinary = summaryCoef %>% 
  filter(Process == 'Detection') %>%
  filter(Predictor != 'Intercept') %>% 
  ggplot(aes(x = Median, y = Predictor, xmin = CIlow, xmax = CIup)) +
  geom_errorbar(aes(linetype = hasZero), width = 0.2) +
  geom_point(aes(shape = hasZero)) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  # scale_x_break(c(0.4, 2.62), ticklabels = c(2.7)) +
  labs(title = '(a) Detection', x = '95% HPD interval', y = 'Covariates') +
  guides(linetype = guide_legend(title=""),
         shape = guide_legend(title="")) +
  theme(legend.position = 'none')

pCount = summaryCoef %>%
  filter(Process == 'Refusal cases') %>%
  filter(Predictor != 'Intercept') %>% 
  ggplot(aes(x = Median, y = Predictor, xmin = CIlow, xmax = CIup)) +
  geom_errorbar(aes(linetype = hasZero), width = 0.2) +
  geom_point(aes(shape = hasZero)) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  # scale_x_break(c(0.25, 0.91), ticklabels = c(0.95)) +
  labs(title = '(b) Refusal cases', x = '95% HPD interval', y = 'Covariates')+
  guides(linetype = guide_legend(title=""),
         shape = guide_legend(title="")) +
  theme(legend.position = 'none')


plot_both = pBinary + pCount
plot_both

ggsave(plot = plot_both, width = 9, height = 3.5, file = 'fig/refuseCoef.pdf')




# =============================================================================-
# Compute RQRs ----
# =============================================================================-
source('src/RFtns.R')

dat.est = data %>%
  dplyr::select(year, month, fips, y, cov, under5) %>%
  add_column(pii = piihat, mu = muhat, nu = nuhat, est = muhat * piihat,
             rqr = stdRQR_ZINB(data$y, muhat, thetahat, piihat))


plot.rqr.scatter = dat.est %>%
  filter(rqr != Inf) %>%
  ggplot(aes(x = log(mu), y = rqr)) +
  geom_point(size = 0.3) +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = -3, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 3, color = 'red', linetype = 'dashed') +
  theme(strip.text = element_text(size = 13)) +
  labs(x = 'Log of mean estimate', y = 'Standardized RQR', title = '(a) Residual plot')

plot.rqr.qq = dat.est %>%
  filter(rqr != Inf) %>%
  ggplot(aes(sample = rqr)) +
  stat_qq(size = 0.3) +
  stat_qq_line(color = 'red') +
  theme(strip.text = element_text(size = 13)) +
  labs(x = 'Standard normal quantile', y = 'Standardized RQR', title = '(b) Quantile-quantile plot')

rqr.final = ggarrange(plot.rqr.scatter, plot.rqr.qq, nrow = 1)

ggsave(plot = rqr.final, width = 6, height = 3.2, file = 'fig/refuseRQR.eps')




# =============================================================================-
# Plot maps ----
# =============================================================================-


