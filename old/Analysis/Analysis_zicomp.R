rm(list=ls())
library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
# q1 = 25; q2 = 300; k = 6
q1 = 25; q2 = 600; k = 6
load(paste0('Analysis/real5_yearly_A_BSFZICOMP', q1, '_BSFvs', q2, '_map.RData'))

claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  # for(j in unique(claim$month)){
    claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i)))
    # claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i, month == j)))
  # }
}


breaks.both = c(round(quantile(c(claimMap$y, claimMap$mode, claimMap$mu),  probs = c(0, 0.6, 0.85, 0.92, 0.96, 0.988, 0.99885))), ceiling(max(c(claimMap$y, claimMap$mode, claimMap$mu))))


tmap_mode('plot')

# observed count
yearly.obs = tm_shape(claimMap) +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = '\nObserved\nincidence\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.obs, 
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_obs.eps',  
          device = cairo_ps, width = 8, height = 5.4)


# approximate mean
yearly.mu = tm_shape(claimMap) +
  tm_polygons(
    col = 'mu',
    palette = 'Greens',
    border.alpha = 0.2,
    breaks = breaks.both,
    title = '\nModel estimate\nof mean\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.mu,
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_mean.eps',
          device = cairo_ps, width = 8, height = 5.4)


# mode
yearly.mode = tm_shape(claimMap) +
  tm_polygons(
    col = 'mode',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = '\nModel estimate\nof mode\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


tmap_save(tm = yearly.mode, 
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_mode.eps',
          device = cairo_ps, width = 8, height = 5.4)


# spatial effect
yearly.V = tm_shape(claimMap) +
  tm_polygons(
    col = 'V',
    midpoint = 0,
    border.alpha = 0.2,
    title = '\nModel estimate\nof spatial effect\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


tmap_save(tm = yearly.V, 
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_V.eps',
          device = cairo_ps, width = 8, height = 5.4)



# dispersion
fivenum_nu = c(unique(round(fivenum(claimMap$nu), 2)), ceiling(max(claimMap$nu*10))/10)

yearly.nu = tm_shape(claimMap) +
  tm_polygons(
    col = 'nu',
    palette = 'Greens', 
    # breaks = breaks.nu,
    border.alpha = 0.2,
    breaks = fivenum_nu,
    title = '\nModel estimate\nof dispersion\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


tmap_save(tm = yearly.nu, 
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_nu.eps',
          device = cairo_ps, width = 8, height = 5.4)


# detection probability
yearly.pii = tm_shape(claimMap) +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    n = 5,
    title = '\nModel estimate of\ndetection probability\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.pii, 
          filename = 'Analysis/Figures_BSFZICOMP_BSFvs/yearly_map_pii.eps',
          device = cairo_ps, width = 8, height = 5.4)

# month = 1
# month.pii = claimMap %>% 
#   filter(month == month) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'pii',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     n = 5,
#     title = '\nModel estimate of\ndetection probability\n'
#   ) +
#   tm_facets(by = 'year', ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   # tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3)
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)
# 
# tmap_save(tm = month.pii, 
#           filename = paste0('Analysis/Figures_zicomp_bsf/month_', month, '_map_pii.eps'), 
#           device = cairo_ps, width = 8, height = 5.4)



# # lambda with state names
# yearly.lam_state = tm_shape(claimMap) +
#   tm_polygons(
#     col = 'lambda',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     # breaks = breaks.obs,
#     breaks = breaks.lam,
#     title = '\nEstimated latent\nmean incidence\n'
#   ) +
#   tm_facets(by = 'year', ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   tm_shape(states_sf) + tm_borders(col = 'gray50', lwd = 0.3) + 
#   tm_text("NAME", size = 0.3)
# 
# tmap_save(tm = yearly.lam_state, 
#           filename = 'Analysis/figures/high_yearly_map_lam_count_state_name.eps', 
#           device = cairo_ps, width = 8, height = 5.4)
# 
# 
# # spatial effects with state names
# yearly.eta2_state = tm_shape(claimMap) +
#   tm_polygons(
#     col = 'eta2',
#     border.alpha = 0.2,
#     midpoint = 0,
#     n = 8,
#     title = '\nSpatial\nrandom effect\nin incidence\n'
#   ) +
#   tm_facets(by = 'year', ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3) +
#   tm_text("NAME", size = 0.3)
# 
# 
# tmap_save(tm = yearly.eta2_state, 
#           filename = 'Analysis/figures/high_yearly_map_count_spatialeffects_state_name.eps', 
#           device = cairo_ps, width = 8, height = 5.4)
# 
# 
# 
# 
# # ------------------------------------------------
# # income
# # ------------------------------------------------
# load('Analysis/refusal.RData')
# refusal$fips = sapply(1:nrow(refusal), function(i) paste0(digits(refusal$fips[i], n = 5), collapse = ""))
# 
# refusalMap2012 = inner_join(counties_sf, refusal %>% filter(year == 2012))
# refusalMap2013 = inner_join(counties_sf, refusal %>% filter(year == 2013))
# refusalMap2014 = inner_join(counties_sf, refusal %>% filter(year == 2014))
# refusalMap2015 = inner_join(counties_sf, refusal %>% filter(year == 2015))
# refusalMap = rbind(refusalMap2012, refusalMap2013, refusalMap2014, refusalMap2015)
# 
# 
# head(refusalMap)
# unique(refusalMap$NAME)[grep('Anne', unique(refusalMap$NAME))]
# refusalMap %>% filter(NAME == 'New York') %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80) %>% slice(c(1, 13, 25, 37))
# refusalMap %>% filter(NAME == 'Anne Arundel') %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80) %>% slice(c(1, 13, 25, 37))
# 
# refusalMap %>% arrange(desc(pc_pop_highinc80)) %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80, NAME)
# refusal %>% arrange(desc(pc_pop_highinc80)) %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80)
# 
# tmap_mode('plot')
# 
# plot_highinc = tm_shape(refusalMap) +
#   tm_polygons(
#     col = 'pc_pop_highinc80',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
#     title = 'Upper 80%'
#   ) +
#   tm_facets(by = 'year', ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   tm_shape(states_sf) + tm_borders(col = 'gray30', lwd = 0.3)
# 
# plot_highlow = tm_shape(refusalMap) +
#   tm_polygons(
#     col = 'pc_pop_lowinc20',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
#     title = 'Lower 20%'
#   ) +
#   tm_facets(by = 'year', ncol = 2) +
#   tm_layout(legend.outside.size = 0.18,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   tm_shape(states_sf) + tm_borders(col = 'gray30', lwd = 0.3)



