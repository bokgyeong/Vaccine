rm(list=ls())
library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
# q1 = 100; q2 = 50
q1 = 25; q2 = 25
# load(paste0('Analysis/real5_yearly_At_BSFZICOMP', q1, '_BSFvst', q2, '_map.RData'))
load(paste0('Analysis/real5_yearly_A_BSFZICOMP', q1, '_BSFvs', q2, '_map.RData'))

claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i)))
}

Year = 2015
claimMap = claimMap %>% filter(year == Year)

# breaks.both = c(round(quantile(c(claimMap$y, claimMap$mode),  probs = c(0, 0.85, 0.92, 0.96, 0.988, 0.99885))), ceiling(max(c(claimMap$y, claimMap$mode))))
# breaks.both = c(round(quantile(c(claimMap$y, claimMap$mode),  probs = c(0, 0.88, 0.92, 0.96, 0.988, 0.99885))), ceiling(max(c(claimMap$y, claimMap$mode))))

breaks.both = c(round(quantile(c(claimMap$y, claimMap$mode, claimMap$muhat),  probs = c(0, 0.6, 0.85, 0.92, 0.96, 0.988, 0.99885))), ceiling(max(c(claimMap$y, claimMap$mode, claimMap$muhat))))

tmap_mode('plot')

# observed count
yearly.obs = tm_shape(claimMap) +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('(a) Observation ('*y*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# detection probability
yearly.pii = tm_shape(claimMap) +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    n = 5,
    title = expression('(b) Probability ('*pi*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# mode
yearly.mode = tm_shape(claimMap) +
  tm_polygons(
    col = 'mode',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('(c) Mode ('*eta*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# fixed
yearly.fixed = tm_shape(claimMap) +
  tm_polygons(
    col = 'fixedeff',
    midpoint = 0,
    border.alpha = 0.2,
    title = expression('(c) Fixed effect ('*X*beta[2]*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)



# random effect in log(mode)
yearly.V = tm_shape(claimMap) +
  tm_polygons(
    col = 'V',
    midpoint = 0,
    border.alpha = 0.2,
    title = expression('(d) Spatial effect ('*V*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# mode
yearly.muhat = tm_shape(claimMap) +
  tm_polygons(
    col = 'muhat',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('(c) Mean ('*mu*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)



# dispersion
# fivenum_nu = c(unique(round(fivenum(claimMap$nu), 2)), round(max(claimMap$nu),1))
fivenum_nu = c(unique(round(fivenum(claimMap$nu), 2)), ceiling(max(claimMap$nu*10))/10)

yearly.nu = tm_shape(claimMap) +
  tm_polygons(
    col = 'nu',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = fivenum_nu,
    # n = 4,
    title = expression('(e) Dispersion ('*nu*')')
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.55,
            frame = FALSE) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)



# tmap_arrange(yearly.mode, yearly.fixed, yearly.V, ncol = 1)

map.final = tmap_arrange(yearly.pii, yearly.muhat, yearly.mode, yearly.V, yearly.nu, ncol = 1)
# map.final = tmap_arrange(yearly.pii, yearly.muhat, yearly.V, yearly.nu, ncol = 1)


tmap_save(tm = yearly.obs, 
          filename = paste0('Analysis/Poster/real5_', Year, 'obs.eps'),
          device = cairo_ps, width = 4, height = 2)

tmap_save(tm = map.final, 
          filename = paste0('Analysis/Poster/real5_', Year, '.eps'),
          device = cairo_ps, width = 4, height = 10)

