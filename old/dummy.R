rm(list=ls())
library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load('Analysis/real_monthly.RData')
claim = data


claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  for(j in unique(claim$month)){
    claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i, month == j)))
  }
}


# mean / variance for positive counts
claimMap.disp = claimMap %>% filter(y > 0) %>% group_by(fips) %>% summarise(mean = mean(y), var = var(y), dispersion = mean(y)/var(y))
claimMap.disp = claimMap.disp %>% filter(!is.infinite(dispersion))
claimMap.disp = claimMap.disp %>% na.omit()
breaks.disp = c(round(quantile(claimMap.disp$dispersion,  probs = c(0, 0.3, 0.485, 0.6, 0.75, 0.85, 0.92, 0.96, 0.99141, 0.9979904)), 2), ceiling(max(claimMap.disp$dispersion)))

obs.disp.monthly = claimMap.disp %>% 
  tm_shape() +
  tm_polygons(
    col = 'dispersion',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.disp,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.text.size = 0.7,
            frame = FALSE,
            main.title = '(b) Ratio of mean to variance for monthly positive counts', 
            main.title.position = 'left', 
            main.title.size = 1) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)





niter = 800000
q1 = q2 = 400; k = 6
load(paste0('Analysis/real5_yearly_A_BSFZICOMP', q1, '_BSFvs', q2, '_n', niter, '_map.RData'))

claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap.yearly = data.frame()

for(i in unique(claim$year)){
  claimMap.yearly = rbind(claimMap.yearly, inner_join(counties_sf, claim %>% filter(year == i)))
}



# mean / variance for positive counts
claimMap.disp.yearly = claimMap.yearly %>% filter(y > 0) %>% group_by(fips) %>% summarise(mean = mean(y), var = var(y), dispersion = mean(y)/var(y))
claimMap.disp.yearly = claimMap.disp.yearly %>% filter(!is.infinite(dispersion))
claimMap.disp.yearly = claimMap.disp.yearly %>% na.omit()
breaks.disp.yearly = c(round(quantile(claimMap.disp.yearly$dispersion,  probs = c(0, 0.3, 0.485, 0.6, 0.75, 0.85, 0.92, 0.96, 0.99141, 0.9979904)), 2), ceiling(max(claimMap.disp.yearly$dispersion)))

obs.disp.yearly = claimMap.disp.yearly %>% 
  tm_shape() +
  tm_polygons(
    col = 'dispersion',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.disp.yearly,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.text.size = 0.7,
            frame = FALSE,
            main.title = '(a) Ratio of mean to variance for yearly positive counts', 
            main.title.position = 'left', 
            main.title.size = 1) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


plot.ratio = tmap_arrange(obs.disp.yearly, obs.disp.monthly, ncol = 2)


tmap_save(plot.ratio,
          filename = 'Analysis/Slides/ratio_yearly_monthly.eps',
          device = cairo_ps, width = 9, height = 3)
