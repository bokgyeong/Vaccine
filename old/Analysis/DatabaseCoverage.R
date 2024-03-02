library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')


# ------------------------------------------------
# Map
# ------------------------------------------------
load('Analysis/new_fit_zinb_bsf_100_map.RData')
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap2012 = inner_join(counties_sf, claim %>% filter(year == 2012))
claimMap2013 = inner_join(counties_sf, claim %>% filter(year == 2013))
claimMap2014 = inner_join(counties_sf, claim %>% filter(year == 2014))
claimMap2015 = inner_join(counties_sf, claim %>% filter(year == 2015))
claimMap = rbind(claimMap2012, claimMap2013, claimMap2014, claimMap2015)

breaks = c(0, summary(claimMap$y[-which(claimMap$y == 0)])[-4])
breaks_cov = c(0, summary(claimMap$cov)[-4])

tmap_mode('plot')

tm_shape(claimMap) +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    breaks = breaks,
    title = 'Number of refusals'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.21)

tm_shape(claimMap) +
  tm_polygons(
    col = 'cov',
    palette = 'Greens', 
    breaks = breaks_cov,
    title = 'Database Coverage'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.21)
