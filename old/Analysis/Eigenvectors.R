library(tidyverse); library(egg); library(grid)
library(ngspatial)

n = 30
A = adjacency.matrix(n)
x = y = seq(0, 1, length.out = n)
X = cbind(rep(x, times = n), rep(y, each = n))
p = ncol(X)
n = nrow(A)


### Moran basis matrix
ones = rep(1, n)
Imat = ones %*% solve(t(ones) %*% ones) %*% t(ones)
Ic = diag(1, n) - Imat
M = Ic %*% A %*% Ic
eig = eigen(M)
eigen = eig$vectors

data = data.frame(x1 = X[,1], x2 = X[,2], eigen = eigen[,1:10])

pEigen1 = ggplot(data, aes(x1, x2, fill = eigen.1)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 1', x = '', y = '') +
    theme(legend.position = 'none')

pEigen2 = ggplot(data, aes(x1, x2, fill = eigen.2)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 2', x = '', y = '') +
    theme(legend.position = 'none')

pEigen3 = ggplot(data, aes(x1, x2, fill = eigen.3)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 3', x = '', y = '') +
    theme(legend.position = 'none')

pEigen4 = ggplot(data, aes(x1, x2, fill = eigen.4)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 4', x = '', y = '') +
    theme(legend.position = 'none')

pEigen5 = ggplot(data, aes(x1, x2, fill = eigen.5)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 5', x = '', y = '') +
    theme(legend.position = 'none')

pEigen6 = ggplot(data, aes(x1, x2, fill = eigen.6)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 6', x = '', y = '') +
    theme(legend.position = 'none')

pEigen7 = ggplot(data, aes(x1, x2, fill = eigen.7)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 7', x = '', y = '') +
    theme(legend.position = 'none')

pEigen8 = ggplot(data, aes(x1, x2, fill = eigen.8)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'white', mid = 'gray', high = "black") +
    guides(fill = guide_coloursteps()) +
    labs(title = 'Eigenvector # 8', x = '', y = '') +
    theme(legend.position = 'none')


ggarrange(pEigen1, pEigen2, pEigen3, pEigen4, 
          pEigen5, pEigen6, pEigen7, pEigen8, nrow = 2) # 800*420

