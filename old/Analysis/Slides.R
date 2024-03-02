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
# q1 = 25; q2 = 600; k = 6
# load(paste0('Analysis/real5_yearly_A_BSFZICOMP', q1, '_BSFvs', q2, '_map.RData'))

niter = 800000
q1 = q2 = 400; k = 6
load(paste0('Analysis/real5_yearly_A_BSFZICOMP', q1, '_BSFvs', q2, '_n', niter, '_map.RData'))

# load('Analysis/real_monthly.RData')
# claim = data


claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  # for(j in unique(claim$month)){
  claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i)))
  # claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i, month == j)))
  # }
}

breaks.both = c(round(quantile(c(claimMap$y, claimMap$mode, claimMap$mu),  probs = c(0, 0.6, 0.85, 0.92, 0.96, 0.988, 0.99885))), ceiling(max(c(claimMap$y, claimMap$mode, claimMap$mu))))


year = 2015

### observations in a year and mean/variance -----------------------------
obs.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.25,
            legend.title.size = 0.62,
            legend.text.size = 0.7,
            frame = FALSE,
            main.title = 'Observations in 2015', 
            # main.title = '(a) Observations in 2015', 
            main.title.position = 'left', 
            main.title.size = 1.2) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# histogram
# obs.hist = claimMap %>% 
#   filter(year == year) %>% 
#   ggplot(aes(x = y)) +
#   geom_histogram(color = "black", fill = "white", binwidth = 20) +
#   labs(x = 'Observed incidence in 2015', y = 'Count')

# mean / variance for positive counts
claimMap.disp = claimMap %>% filter(y > 0) %>% group_by(fips) %>% summarise(dispersion = mean(y)/var(y))
claimMap.disp = claimMap.disp %>% filter(!is.infinite(dispersion))
claimMap.disp = claimMap.disp %>% na.omit()
breaks.disp = c(round(quantile(claimMap.disp$dispersion,  probs = c(0, 0.3, 0.485, 0.6, 0.75, 0.85, 0.92, 0.96, 0.99141, 0.9979904)), 2), ceiling(max(claimMap.disp$dispersion)))

obs.disp = claimMap.disp %>% 
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
            main.title = '(b) Ratio of mean to variance for positive counts', 
            main.title.position = 'left', 
            main.title.size = 1.2) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# plot.obs = grid.arrange(tmap_grob(obs.map), obs.hist, ncol = 2)

# ggsave(plot.obs, 
#        filename = 'Analysis/Slides/obs.eps', 
#        device = cairo_ps, width = 7.5, height = 2.8)

plot.obs.disp = tmap_arrange(obs.map, obs.disp, nrow = 1)
tmap_save(tm = plot.obs.disp,
          filename = 'Analysis/JSM2022/obs_and_dispersion.eps',
          device = cairo_ps, width = 9, height = 3)

tmap_save(tm = obs.map,
          filename = 'Analysis/JSM2022/obs.eps',
          device = cairo_ps, width = 4.5, height = 3)




### observations yearly -------------------------------------------
obs.yearly = tm_shape(claimMap) +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression(y)
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.17,
            legend.title.size = 1.2,
            legend.text.size = 0.9,
            panel.label.size = 1.5) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = obs.yearly,
          filename = 'Analysis/JSM2022/obs_yearly.eps',
          device = cairo_ps, width = 8, height = 5)



### Model estimates ----------------------------------------------------------
obs.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression(y)
  ) +
  tm_layout(main.title = '(a) Observations', 
            main.title.position = 'left', 
            main.title.size = 1.4, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.23,
            legend.text.size = 0.58) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# mean
mu.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'mu',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression(mu)
  ) +
  tm_layout(main.title = '(b) Mean', 
            main.title.position = 'left', 
            main.title.size = 1.4, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.23,
            legend.text.size = 0.58) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# lambda
# mode.map = claimMap %>% 
#   filter(year == year) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'mode',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     breaks = breaks.both,
#     title = ''
#   ) +
#   tm_layout(main.title = 'Mode', main.title.position = 'center', 
#             main.title.size = 1.2, legend.text.size = 0.4) +
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# pii
pii.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    title = expression(pi)
  ) +
  tm_layout(main.title = '(c) Probability', 
            main.title.position = 'left', 
            main.title.size = 1.4, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.23,
            legend.text.size = 0.58) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# spatial effect
# V.breaks = seq(-0.6, 0.45, by = 0.15)
V.breaks = round(seq(-0.6, 0.4, by = 0.2), 2)

V.map = claimMap %>%
  filter(year == year) %>%
  tm_shape() +
  tm_polygons(
    col = 'V',
    midpoint = 0,
    border.alpha = 0.2,
    n = 7,
    breaks = V.breaks,
    title = 'V'
  ) +
  tm_layout(main.title = '(e) Spatial effect in refusal',
            main.title.position = 'left', 
            main.title.size = 1.4, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.23,
            legend.text.size = 0.58) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# nu
fivenum_nu = c(unique(round(fivenum(claimMap$nu), 2)), ceiling(max(claimMap$nu*10))/10)

nu.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'nu',
    palette = 'Greens',
    border.alpha = 0.2,
    breaks = fivenum_nu,
    title = expression(nu)
  ) +
  tm_layout(main.title = '(d) Dispersion',
            main.title.position = 'left', 
            main.title.size = 1.4, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.23,
            legend.text.size = 0.58) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)




# bsf zinb
# plot.comp = tmap_arrange(obs.map, lam.map, pii.map, eta2.map, ncol = 2)
# tmap_save(tm = plot.comp,
#           filename = 'Analysis/Slides/comp_bsf_zinb.eps',
#           device = cairo_ps, width = 8, height = 5.4)

# bsf zicomp
# plot.zicomp = tmap_arrange(obs.map, mu.map, mode.map, V.map, pii.map, nu.map, nrow = 2)
# tmap_save(tm = plot.zicomp,
#           # filename = paste0('Analysis/Slides/BSFZICOMP', q1, '_BSFvs', q2, '_map.eps'),
#           filename = paste0('Analysis/Slides/RJMCMCBSFZICOMP', q1, '_BSFvs', q2, '_n', niter, '_map.eps'),
#           device = cairo_ps, width = 9, height = 5)

plot.zicomp = tmap_arrange(obs.map, mu.map, pii.map, nu.map, V.map, nrow = 2)
tmap_save(tm = plot.zicomp,
          filename = 'Analysis/JSM2022/estimates.eps',
          device = cairo_ps, width = 12, height = 4.6)

# constant zicomp
# plot.comp = tmap_arrange(obs.map, mode.map, pii.map, ncol = 2)
# tmap_save(tm = plot.comp,
#           filename = 'Analysis/Slides/comp_con_zicomp.eps',
#           device = cairo_ps, width = 8, height = 5.4)







nu.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'nu',
    palette = 'Greens',
    border.alpha = 0.2,
    breaks = fivenum_nu,
    title = expression(nu)
  ) +
  tm_layout(main.title = 'Dispersion estimates',
            main.title.position = 'left', 
            main.title.size = 1.2, 
            frame = F,
            legend.outside = T,
            legend.outside.position = 'right',
            legend.outside.size = 0.25,
            legend.text.size = 0.8) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = nu.map,
          filename = 'Analysis/Slides/estimates_nu.eps',
          device = cairo_ps, width = 5, height = 3)





# ------------------------------------------------
# Comparision between methods
# ------------------------------------------------

# bsf zinb 
load('Analysis/high_fit_block_zinb_bsf_150_map.RData')
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))
claimMap_bsf_zinb = data.frame()
for(i in unique(claim$year)){
  claimMap_bsf_zinb = rbind(claimMap_bsf_zinb, inner_join(counties_sf, claim %>% filter(year == i)))
}
claimMap_bsf_zinb = claimMap_bsf_zinb %>% add_column(Model = '(a) ZINB') %>% select(pii, lambda, year, Model)

# bsf zicomp
q = 100
k = 6
# load(paste0('Analysis/real_yearly_ZICOMP_BSFvs', q, '_k', k, '_map.RData'))
load(paste0('Analysis/real_yearly_ZICOMP_BSFvs', q, '_k', k, '_long_map.RData'))
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))
claimMap_bsf_zicomp = data.frame()
for(i in unique(claim$year)){
  claimMap_bsf_zicomp = rbind(claimMap_bsf_zicomp, inner_join(counties_sf, claim %>% filter(year == i)))
}
claimMap_bsf_zicomp = claimMap_bsf_zicomp %>% add_column(Model = '(b) ZICOMP') %>% select(pii, mode, year, Model)
colnames(claimMap_bsf_zicomp) = c('pii', 'lambda', 'year', 'Model', 'geometry')


claimMap = rbind(claimMap_bsf_zinb, claimMap_bsf_zicomp)
breaks.loc = round(quantile(claimMap$lambda,  probs = c(0, 0.88, 0.96, 0.988, 0.99885, 1)))


year = 2015


# pii
pii.comp.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    title = 'Probability'
  ) +
  tm_facets(by = 'Model', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = pii.comp.map,
          filename = 'Analysis/Slides/pii_zinb_mode_zicomp.eps',
          device = cairo_ps, width = 8, height = 2.7)

# location
loc.comp.map = claimMap %>% 
  filter(year == year) %>% 
  tm_shape() +
  tm_polygons(
    col = 'lambda',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.loc,
    title = 'Location'
  ) +
  tm_facets(by = 'Model', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = loc.comp.map,
          filename = 'Analysis/Slides/loc_zinb_mode_zicomp.eps',
          device = cairo_ps, width = 8, height = 2.7)
