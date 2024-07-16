rm(list = ls())
# require(coda); 
require(tidyverse); require(egg); 
# require(grid)
require(patchwork)
require(batchmeans); require(foreach); 
# require(doParallel); require(parallel)
# require(ggbreak)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# directory ----
dirn = 'nszinb/'
dirn.data = paste0(dirn, 'data/')
dirn.fit = paste0(dirn, 'fit/')
dirn.fig = paste0(dirn, 'fig/')
ifelse(!dir.exists(dirn.fig), dir.create(dirn.fig, recursive = T), FALSE)


# summarize results ----
load(paste0(dirn.data, 'sim.RData'))
load(paste0(dirn.fit, 'simZINB.RData'))

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



# figures for coefficients ----
Predictors = paste0('Predictor ', 1:2)

summaryCoef = data.frame(
  Process = 'Detection', 
  Median = summat[beta1Ind,2], 
  CIlow = summat[beta1Ind,4],
  CIup = summat[beta1Ind,5],
  Predictor = Predictors)

summaryCoef = rbind(
  summaryCoef,
  data.frame(
    Process = 'Count', 
    Median = summat[beta2Ind,2], 
    CIlow = summat[beta2Ind,4],
    CIup = summat[beta2Ind,5],
    Predictor = Predictors))

summaryCoef$Predictor = factor(summaryCoef$Predictor, levels = rev(unique(summaryCoef$Predictor)))
summaryCoef = summaryCoef %>% mutate(hasZero = ifelse(CIlow <= 0 & CIup >= 0, 'Include zero', 'Do not include zero'))

pBinary = summaryCoef %>% 
  filter(Process == 'Detection') %>%
  ggplot(aes(x = Median, y = Predictor, xmin = CIlow, xmax = CIup)) +
  geom_errorbar(aes(linetype = hasZero), width = 0.2) +
  geom_point(aes(shape = hasZero)) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  # scale_x_break(c(0.4, 2.62), ticklabels = c(2.7)) +
  labs(title = '(a) Detection', x = '95% HPD interval', y = '') +
  guides(linetype = guide_legend(title=""),
         shape = guide_legend(title="")) +
  theme(legend.position = 'none')

pCount = summaryCoef %>%
  filter(Process == 'Count') %>%
  ggplot(aes(x = Median, y = Predictor, xmin = CIlow, xmax = CIup)) +
  geom_errorbar(aes(linetype = hasZero), width = 0.2) +
  geom_point(aes(shape = hasZero)) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  # scale_x_break(c(0.25, 0.91), ticklabels = c(0.95)) +
  labs(title = '(b) Count', x = '95% HPD interval', y = '')+
  guides(linetype = guide_legend(title=""),
         shape = guide_legend(title="")) +
  theme(legend.position = 'none')


plot_both = pBinary + pCount
plot_both

ggsave(paste0(dirn.fig, 'simCoef.png'), plot_both, width = 6.2, height = 3)




# residual analysis ----
source('src/RFtns.R')

dat.est = data.frame(
  y = y, pii = piihat, mu = muhat, theta = thetahat, est = muhat * piihat
) %>% mutate(
  rqr = stdRQR_ZINB(y, mu, theta, pii)
)

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

ggsave(paste0(dirn.fig, 'simRQR.png'), rqr.final, width = 6, height = 3.2)


