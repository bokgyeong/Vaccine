library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load('Analysis/real_yearly_mode_zicmp_bsf_100_map.RData')
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i)))
}


### observed vs estimate --------------------------
# claimMap = claimMap %>% pivot_longer(
#   cols = c('y', 'est'),
#   names_to = 'Subject',
#   values_to = 'Incidence'
# )
# breaks.both = round(quantile(claimMap$Incidence,  probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885, 1)))
# 
# claimMap$Subject = factor(claimMap$Subject, 
#                           levels = c('y', 'est'),
#                           labels = c('Observed', 'Estimate'))
# 
# tmap_mode('plot')
# 
# plot.y.est.1 = claimMap %>% 
#   filter(year %in% c(2012, 2013)) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'Incidence',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     breaks = breaks.both,
#     title = 'Iincidence'
#   ) +
#   tm_facets(by = c('year', 'Subject'), ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   # tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)
# 
# tmap_save(tm = plot.y.est.1, 
#           filename = 'Analysis/Figures_zicomp_bsf/yearly_map_obs_est_1.eps', 
#           device = cairo_ps, width = 8, height = 5.4)
# 
# 
# plot.y.est.2 = claimMap %>% 
#   filter(year %in% c(2014, 2015)) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'Incidence',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     breaks = breaks.both,
#     title = 'Iincidence'
#   ) +
#   tm_facets(by = c('year', 'Subject'), ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   # tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)
# 
# tmap_save(tm = plot.y.est.2, 
#           filename = 'Analysis/Figures_zicomp_bsf/yearly_map_obs_est_2.eps', 
#           device = cairo_ps, width = 8, height = 5.4)


### observed vs mode --------------------------
claimMap = claimMap %>% pivot_longer(
  cols = c('y', 'mode'),
  names_to = 'Subject',
  values_to = 'Incidence'
)
breaks.both = round(quantile(claimMap$Incidence,  probs = c(0, 0.80, 0.88, 0.96, 0.988, 0.99885, 1)))

claimMap$Subject = factor(claimMap$Subject, 
                          levels = c('y', 'mode'),
                          labels = c('Observed', 'Mode'))

tmap_mode('plot')

plot.y.mode.1 = claimMap %>% 
  filter(year %in% c(2012, 2013)) %>% 
  tm_shape() +
  tm_polygons(
    col = 'Incidence',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = 'Incidence'
  ) +
  tm_facets(by = c('year', 'Subject'), ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  # tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = plot.y.mode.1, 
          filename = 'Analysis/Figures_zicomp_bsf/yearly_map_obs_mode_1.eps', 
          device = cairo_ps, width = 8, height = 5.4)


plot.y.mode.2 = claimMap %>% 
  filter(year %in% c(2014, 2015)) %>% 
  tm_shape() +
  tm_polygons(
    col = 'Incidence',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = 'Incidence'
  ) +
  tm_facets(by = c('year', 'Subject'), ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  # tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = plot.y.est.2, 
          filename = 'Analysis/Figures_zicomp_bsf/yearly_map_obs_mode_2.eps', 
          device = cairo_ps, width = 8, height = 5.4)

