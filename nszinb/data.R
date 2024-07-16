rm(list=ls())
require(tidyverse); require(ngspatial); require(foreach); require(MASS)
require(egg)

# =============================================================================-
# simulate data ----
# =============================================================================-

## directory ----
filepath.data = 'nszinb/data/'
ifelse(!dir.exists(filepath.data), dir.create(filepath.data, recursive = T), FALSE)


## design matrix ----
m = 30
A = adjacency.matrix(m)
x = y = seq(0, 1, length.out = m)
coord = cbind(rep(x, times = m), rep(y, each = m))
t = 2
X = foreach(j = 1:t, .combine = 'rbind') %do% { cbind(coord) }
p = ncol(X); N = nrow(X)
n = rep(nrow(A), t); nt = c(0, cumsum(n))


## true parameters ----
beta1 = c(-2, 1)
beta2 = c(2, 2)

eta = foreach(j = 1:t, .combine = 'c') %do% { as.vector( X[(nt[j]+1):(nt[j+1]),] %*% beta1 ) }
pii = exp(eta)/(1 + exp(eta))
mu = foreach(j = 1:t, .combine = 'c') %do% { exp( as.vector( X[(nt[j]+1):(nt[j+1]),] %*% beta2 ) ) }
theta = 1

y = rbinom(n = N, size = 1, prob = pii)
y[y == 1] = rnegbin(sum(y == 1), mu = mu[y == 1], theta = theta)

trpar = list(
  beta1 = beta1, beta2 = beta2,
  pii = pii, mu = mu, theta = theta
)

save(coord, t, n, nt, N, p, X, trpar, y, file = paste0(filepath.data, 'sim.RData'))



# =============================================================================-
# figures ----
# =============================================================================-

## directory ----
filepath.fig = 'nszinb/fig/'
ifelse(!dir.exists(filepath.fig), dir.create(filepath.fig, recursive = T), FALSE)

dat = data.frame()
load(paste0(filepath.data, 'sim.RData'))
for(i in 1:t) {
  dat = rbind(
    dat,
    data.frame(
      Year = paste0('Year ', i), 
      x1 = coord[,1], x2 = coord[,2], 
      y = y[(nt[i]+1):(nt[i+1])],
      mu = trpar$mu[(nt[i]+1):(nt[i+1])],
      pii = trpar$pii[(nt[i]+1):(nt[i+1])]
    )
  )
}  

plot.y = dat %>%
  filter(Year == 'Year 1') %>% 
  ggplot(aes(x1, x2, fill = y)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'white', high = "black") +
  labs(x = expression(x[1]), y = expression(x[2]), fill = expression(y), title = '(a) Observations for Year 1') +
  theme_bw() +
  theme(legend.position = 'bottom')

plot.pii = dat %>%
  filter(Year == 'Year 1') %>% 
  ggplot(aes(x1, x2, fill = pii)) +
  geom_tile(color = 'white') +
  scale_fill_viridis_c() +
  labs(x = expression(x[1]), y = expression(x[2]), fill = expression(pi), title = '(b) Probabilities') +
  theme_bw() +
  theme(legend.position = 'bottom')

plot.mu = dat %>%
  filter(Year == 'Year 1') %>% 
  ggplot(aes(x1, x2, fill = mu)) +
  geom_tile(color = 'white') +
  scale_fill_viridis_c() +
  labs(x = expression(x[1]), y = expression(x[2]), fill = expression(mu), title = '(c) Means') +
  theme_bw() +
  theme(legend.position = 'bottom')

plot.data = ggarrange(plot.y, plot.pii, plot.mu, nrow = 1)

ggsave(paste0(filepath.fig, 'simData.png'), plot.data, width = 7.8, height = 3.3)

