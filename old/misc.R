library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load('Analysis/old_fit_zinb_bsf_100_map.RData')
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))
claimMap2012 = inner_join(counties_sf, claim %>% filter(year == 2012))
claimMap2013 = inner_join(counties_sf, claim %>% filter(year == 2013))
claimMap2014 = inner_join(counties_sf, claim %>% filter(year == 2014))
claimMap2015 = inner_join(counties_sf, claim %>% filter(year == 2015))
claimMap = rbind(claimMap2012, claimMap2013, claimMap2014, claimMap2015)
claimMap = claimMap %>% mutate(obs.rate = y/cov * 100000, est.rate = est/cov * 100000) %>% add_column(Name = '(a) Before')

load('Analysis/new_fit_zinb_bsf_200_map.RData')
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap2012 = inner_join(counties_sf, claim %>% filter(year == 2012))
claimMap2013 = inner_join(counties_sf, claim %>% filter(year == 2013))
claimMap2014 = inner_join(counties_sf, claim %>% filter(year == 2014))
claimMap2015 = inner_join(counties_sf, claim %>% filter(year == 2015))
claimMap_new = rbind(claimMap2012, claimMap2013, claimMap2014, claimMap2015)
claimMap_new = claimMap_new %>% mutate(obs.rate = y/cov * 100000, est.rate = est/cov * 100000) %>% add_column(Name = '(b) After')
claimMap = rbind(claimMap, claimMap_new)


tmap_mode('plot')

yearly.eta2 = tm_shape(claimMap) +
  tm_polygons(
    col = 'eta2',
    border.alpha = 0.2,
    midpoint = 0,
    n = 8,
    title = '\nSpatial\nrandom effect\nin incidence\n'
  ) +
  tm_facets(by = c('year', 'Name'), ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)

tmap_save(tm = yearly.eta2, 
          filename = 'Analysis/figures/yearly_map_count_spatialeffects_before_after_law.eps', 
          device = cairo_ps, width = 7.5, height = 10)

