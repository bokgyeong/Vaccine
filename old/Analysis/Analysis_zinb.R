library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load('Analysis/map_refuse_yearly_A_ZINB_k6.RData')

# q = 25
# q = 250
# q = 300
# load(paste0('Analysis/map_refuse_yearly_A_BSFZINB', q, '_k6.RData'))
# load(paste0('Analysis/map_refuse_yearly_A_BSFZINBmst', q, '_k6.RData'))



claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))


### for relative risk
claim = claim %>% mutate(ExpNum = sum(y) / sum(under5) * under5)


### for standardized residuals
claim = claim %>% mutate(std.rsd = (y - est) / sd(y - est))
# claim = claim %>% mutate(std.rsd = (y - lambda) / sd(y - lambda))


claimMap = data.frame()
for(i in unique(claim$year)){
  claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i)))
}

claimMap = claimMap %>% mutate(lampercov = lambda / cov * 1000, lamperpop = lambda / under5 * 1000, lamrr = lambda / ExpNum,
                               ypercov = y / cov * 1000, yperpop = y / under5 * 1000, yrr = y / ExpNum)

breaks.both = c(round(quantile(c(claimMap$y, claimMap$lambda),  probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                ceiling(max(c(claimMap$y, claimMap$lambda)*10))/10)

breaks.lampercov = c(trunc(min(c(claimMap$lampercov)*10))/10, round(quantile(c(claimMap$lampercov),  probs = c(0.75, 0.88, 0.96, 0.988, 0.99885)), 2), 
                ceiling(max(c(claimMap$lampercov)*10))/10)

breaks.lamperpop = c(round(quantile(c(claimMap$lamperpop), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                ceiling(max(c(claimMap$lamperpop))))
breaks.yperpop = c(round(quantile(c(claimMap$yperpop), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                     ceiling(max(c(claimMap$yperpop))))
breaks.bothperpop = c(round(quantile(c(claimMap$yperpop, claimMap$lamperpop), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                     ceiling(max(c(claimMap$yperpop, claimMap$lamperpop))))

breaks.lamrr = c(round(quantile(c(claimMap$lamrr), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885)), 2), 
                     ceiling(max(c(claimMap$lamrr))))
breaks.yrr = c(round(quantile(c(claimMap$yrr), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                   ceiling(max(c(claimMap$yrr))))
breaks.bothrr = c(round(quantile(c(claimMap$yrr, claimMap$lamrr), probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885))), 
                      ceiling(max(c(claimMap$yrr, claimMap$lamrr))))



# ===================================================-
# counties with the highest refusal rate ----
# ===================================================-
library(readxl)

salmon = read_excel('Analysis/salmon2015.xlsx')
head(salmon)


inner_join(claimMap, salmon, by = c('NAME' = 'County', 'STATE_NAME' = 'State'))


# ===================================================-
# plot ----
# ===================================================-
tmap_mode('plot')

# observed count
yearly.obs.count = tm_shape(claimMap) +
  tm_polygons(
    col = 'y',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('Observation ('*y*')')
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.obs.count, 
          # filename = 'Analysis/Figures_BSFZINB/obs.eps', 
          filename = 'Analysis/Figures_ZINB/ns_obs.eps', 
          device = cairo_ps, width = 8, height = 5.4)


# estimated mean count
yearly.est.count = tm_shape(claimMap) +
  tm_polygons(
    col = 'est',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('Estimation ('*pi*mu*')')
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.est.count, 
          # filename = 'Analysis/Figures_BSFZINB/estimation.eps', 
          filename = 'Analysis/Figures_ZINB/ns_estimation.eps', 
          device = cairo_ps, width = 8, height = 5.4)



# estimated latent mean count
yearly.lam.count = tm_shape(claimMap) +
  tm_polygons(
    col = 'lambda',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = 'Incidence'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.lam.count, 
          # filename = 'Analysis/Figures_BSFZINB/mu.eps', 
          filename = 'Analysis/Figures_ZINB/ns_mu.eps', 
          device = cairo_ps, width = 8, height = 5.4)


# observed count per 1,000 interactions
yearly.ypercov = tm_shape(claimMap) +
  tm_polygons(
    col = 'ypercov',
    palette = 'Greens', 
    border.alpha = 0.2,
    # breaks = breaks.lampercov,
    title = '\nObserved Incidence per\n1,000 interactions\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# estimated latent mean count per 1,000 interactions
yearly.lampercov = tm_shape(claimMap) +
  tm_polygons(
    col = 'lampercov',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.lampercov,
    # breaks = c(0, 0.5, 3, 5, 10, 15, 20, round(max(claimMap$yperpop)+1)),
    title = '\nIncidence per\n1,000 interactions\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.lampercov, 
          # filename = 'Analysis/Figures_ZINB/ns_mupercov.eps', 
          filename = 'Analysis/Figures_ZINB/refuse_zinb_mupercov.eps', 
          device = cairo_ps, width = 8, height = 5.4)



# observed count per 1,000 children
yearly.yperpop = claimMap %>% 
  filter(year == 2015) %>% 
  tm_shape() +
  tm_polygons(
    col = 'yperpop',
    palette = 'Greens', 
    border.alpha = 0.2,
    # breaks = breaks.yperpop,
    # breaks = breaks.bothperpop,
    # breaks = c(0, 0.5, 3, 5, 10, 15, 20, round(max(claimMap$yperpop)+1)),
    # breaks = c(0, 0.5, 3, 5, 10, 15, 20, round(max(claimMap$lamperpop)+1)),
    # breaks = c(0, 0.5, 3, 5, 10, 20, 100),
    breaks = c(0, 1, 2, 3, 10, 20, 100),
    legend.hist.width = 1,
    title = '\nObserved\ncases per\n1,000 children\n'
  ) +
  # tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.17,
            legend.title.size = 1.2,
            legend.hist.width = 1,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.yperpop, 
          filename = 'Analysis/Figures_ZINB/refuse_zinb_yperpop.eps',
          # filename = 'Analysis/Figures_ZINB/underut_zinb_yperpop.eps',
          device = cairo_ps, width = 8, height = 5.4)



# estimated latent mean count per 1,000 children
yearly.lamperpop = tm_shape(claimMap) +
  tm_polygons(
    col = 'lamperpop',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.lamperpop,
    # breaks = breaks.bothperpop,
    # breaks = c(0, 0.5, 3, 5, 10, 15),
    # breaks = c(0, 0.5, 3, 5, 10, 15, 20, round(max(claimMap$yperpop)+1)),
    # breaks = c(0, 0.5, 3, 5, 10, 13),
    title = '\nExpected\ncases per\n1,000 children\nunder perfect\ndetection\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.17,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.lamperpop, 
          # filename = 'Analysis/Figures_ZINB/ns_muperpop.eps', 
          filename = 'Analysis/Figures_ZINB/refuse_zinb_muperpop.eps',
          # filename = 'Analysis/Figures_ZINB/underut_zinb_muperpop.eps',
          device = cairo_ps, width = 8, height = 5.4)



# observed relative risk
yearly.yrr = tm_shape(claimMap) +
  tm_polygons(
    col = 'yrr',
    palette = 'Greens', 
    border.alpha = 0.2,
    # breaks = breaks.yrr,
    breaks = c(0, 0.2, 0.5, 1, 3, 5, 15, 60),
    # breaks = breaks.bothrr,
    title = '\nObserved\nrelative risk\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.17,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.yrr,
          filename = 'Analysis/Figures_ZINB/refuse_zinb_yrr.eps',
          device = cairo_ps, width = 8, height = 5.4)



# relative risk under perfect detection
yearly.lamrr = tm_shape(claimMap) +
  tm_polygons(
    col = 'lamrr',
    palette = 'Greens', 
    border.alpha = 0.2,
    # breaks = breaks.lamrr,
    breaks = c(0, 0.2, 0.5, 1, 1.5, 2, 3, 5, 8),
    # breaks = c(0, 0.2, 0.5, 1, 3, 5, 15),
    # breaks = breaks.bothrr,
    title = '\nExpected\nRelative risk\nunder perfect\ndetection\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.16,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.lamrr,
          filename = 'Analysis/Figures_ZINB/refuse_zinb_murr.eps',
          device = cairo_ps, width = 8, height = 5.4)



# spatial dependence in refusal cases
yearly.W = claimMap %>% 
  filter(year == 2015) %>%
  tm_shape() +
  tm_polygons(
    col = 'W',
    border.alpha = 0.2,
    midpoint = 0,
    n = 8,
    # title = '\nSpatial dependence\ncaused by\nsocial influence\n'
    title = '\nSpatial\nrandom effects\n'
  ) +
  # tm_facets(by = 'year', ncol = 2) +
  # tm_layout(legend.outside.size = 0.16,
  tm_layout(legend.outside.size = 0.26,
            legend.title.size = 1.2,
            legend.text.size = 0.9,
            legend.outside = T,
            frame = F) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.W,
          filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_W.eps'),
          # filename = paste0('Analysis/Figures_BSFZINBmst/refuse_bsfzinbmst', q, '_W.eps'),
          device = cairo_ps, width = 5.5, height = 3)
          # device = cairo_ps, width = 8, height = 5.4)


yearly.pii = tm_shape(claimMap) +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    # n = 5,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    style = "cont",
    # title = expression('Probability ('*pi*')')
    title = '\nProbability\nof detection\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.pii,
          # filename = 'Analysis/Figures_BSFZINB/pi.eps', 
          filename = 'Analysis/Figures_ZINB/ns_pi.eps', 
          device = cairo_ps, width = 8, height = 5.4)


yearly.rqr = tm_shape(claimMap) +
  tm_polygons(
    col = 'std.rqr',
    border.alpha = 0.2,
    midpoint = 0,
    n = 8,
    title = '\nRandomized\nquantile residuals\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.rqr,
          filename = 'Analysis/Figures_ZINB/refuse_zinb_rqr_map.eps',
          # filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_rqr_map.eps'),
          device = cairo_ps, width = 8, height = 5.4)


yearly.rsd = tm_shape(claimMap) +
  tm_polygons(
    col = 'std.rsd',
    border.alpha = 0.2,
    midpoint = 0,
    # n = 8,
    breaks = c(-20, -10, -3, -1, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 3, 10, 20, 40, 50), 
    title = '\nStandardized\nresiduals\n'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.76) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = yearly.rsd,
          filename = 'Analysis/Figures_ZINB/refuse_zinb_rsd_map.eps',
          # filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_rsd_map.eps'),
          device = cairo_ps, width = 8, height = 5.4)






# ------------------------------------------------------------------------------
# Observations vs estimates
# ------------------------------------------------------------------------------

### relative risk
# claimMapComp = claimMap %>% pivot_longer(
#   cols = c('yrr', 'lamrr'),
#   names_to = 'Subject',
#   values_to = 'RelativeRisk'
# )
# breaks.both = round(quantile(claimMapComp$RelativeRisk,  probs = c(0, 0.80, 0.88, 0.96, 0.988, 0.99885, 1)))
# 
# claimMapComp$Subject = factor(claimMapComp$Subject, 
#                           levels = c('yrr', 'lamrr'),
#                           labels = c('Observed relative risk', 'Expected relative risk under perfect detection'))
# 
# plotComp = claimMapComp %>% 
#   filter(year == 2015) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'RelativeRisk',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     # breaks = breaks.both,
#     breaks = c(0, 0.2, 0.5, 1, 1.5, 2, 3, 5, 8, 15, 60),
#     style = "cont",
#     title = 'Relative risk'
#   ) +
#   tm_facets(by = 'Subject', ncol = 2) +
#   tm_layout(legend.outside.size = 0.16,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9) +
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# tmap_save(tm = plotComp, 
#           filename = 'Analysis/Figures_ZINB/refuse_zinb_yrr_murr.eps',
#           device = cairo_ps, width = 8.5, height = 2.8)




### cases per 1,000 kids
claimMapComp = claimMap %>% pivot_longer(
  cols = c('yperpop', 'lamperpop'),
  names_to = 'Subject',
  values_to = 'CasesPerPop'
)
breaks.both = round(quantile(claimMapComp$CasesPerPop,  probs = c(0, 0.80, 0.88, 0.96, 0.988, 0.99885, 1)))

claimMapComp$Subject = factor(claimMapComp$Subject, 
                              levels = c('yperpop', 'lamperpop'),
                              labels = c('Observed cases', 'Expected cases under perfect detection'))

library(RColorBrewer)
plotComp = claimMapComp %>% 
  # filter(year == 2015, Subject == 'Observed cases') %>%
  # filter(year == 2015, Subject == 'Expected cases under perfect detection') %>%
  tm_shape() +
  tm_polygons(
    col = 'CasesPerPop',
    # palette = 'Greens',
    palette = 'seq',
    border.alpha = 0.2,
    # breaks = breaks.both,
    breaks = c(0, 1, 2, 3, 10, 20, 100),
    # breaks = c(0, 0.5, 3, 5, 10, 15, 20, 100),
    style = "cont",
    # legend.hist = TRUE,
    title = '\nCases per\n1,000 children\n'
  ) +
  tm_facets(by = 'Subject', ncol = 2) +
  # tm_facets(by = c('year', 'Subject'), ncol = 2) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.16,
            # legend.outside.size = 0.3,
            legend.title.size = 1.2,
            legend.text.size = 0.9,
            frame = F,
            legend.hist.width = 0.4,
            legend.hist.size = 0.5,
            aes.palette = list(seq = "-viridis")) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# tmap_save(tm = plotComp,
#           # filename = 'Analysis/Figures_ZINB/refuse_zinb_yperpop_muperpop.eps',
#           filename = 'Analysis/Figures_ZINB/refuse_zinb_yperpop_muperpop_yearly.eps',
#           # filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_yperpop_muperpop.eps'),
#           # device = cairo_ps, width = 8.5, height = 2.8)
#           device = cairo_ps, width = 8.5, height = 9)



### refusal cases
# claimMapComp = claimMap %>% pivot_longer(
#   cols = c('y', 'lambda'),
#   names_to = 'Subject',
#   values_to = 'Cases'
# )
# breaks.both = round(quantile(claimMapComp$Cases,  probs = c(0, 0.75, 0.88, 0.96, 0.988, 0.99885, 1)))
# 
# claimMapComp$Subject = factor(claimMapComp$Subject, 
#                               levels = c('y', 'lambda'),
#                               labels = c('Observed cases', 'Expected cases under perfect detection'))
# 
# plotComp = claimMapComp %>% 
#   filter(year == 2015) %>% 
#   tm_shape() +
#   tm_polygons(
#     col = 'Cases',
#     palette = 'Greens', 
#     border.alpha = 0.2,
#     # breaks = breaks.both,
#     breaks = c(0, 1, 5, 10, 15, 50, 100, 200, 500, 1000, 1500, 2126),
#     style = "cont",
#     title = 'Cases'
#   ) +
#   tm_facets(by = 'Subject', ncol = 2) +
#   tm_layout(legend.outside.size = 0.16,
#             legend.title.size = 1.2,
#             legend.text.size = 0.9,
#             frame = F) +
#   tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# tmap_save(tm = plotComp,
#           filename = 'Analysis/Figures_ZINB/refuse_zinb_y_mu.eps',
#           # filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_y_mu.eps'),
#           device = cairo_ps, width = 8.5, height = 2.8)




### detection prob for 2015
plotPi2015 = claimMap %>% 
  filter(year == 2015) %>% 
  tm_shape() +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    # n = 5,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    style = "cont",
    title = '\nProbability\nof detection\n'
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.16,
            legend.title.size = 1.2,
            legend.text.size = 0.9,
            frame = F) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# tmap_save(tm = plotPi2015, 
#           filename = 'Analysis/Figures_ZINB/refuse_zinb_pi_2015.eps',
#           device = cairo_ps, width = 4.5, height = 2.8)


plotAll2015 = tmap_arrange(plotComp, plotPi2015, ncol = 1, heights = c(1, 0.7))

tmap_save(tm = plotAll2015, 
          filename = 'Analysis/Figures_ZINB/refuse_zinb_yperpop_muperpop_pi.eps',
          # filename = paste0('Analysis/Figures_BSFZINBms/refuse_bsfzinb', q, '_yperpop_muperpop_pi.eps'),
          device = cairo_ps, width = 8.5, height = 5.5)


# ------------------------------------------------------------------------------
# misc
# ------------------------------------------------------------------------------

# lambda with state names
yearly.lam.count_state = tm_shape(claimMap) +
  tm_polygons(
    col = 'lambda',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.both,
    title = expression('Mean ('*mu*')')
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3) +
  tm_text("NAME", size = 0.3)

tmap_save(tm = yearly.lam.count_state, 
          filename = 'Analysis/Figures_BSFZINB/mu_statenames.eps', 
          device = cairo_ps, width = 8, height = 5.4)


# spatial effects with state names
yearly.eta2_state = tm_shape(claimMap) +
  tm_polygons(
    col = 'eta2',
    border.alpha = 0.2,
    midpoint = 0,
    n = 8,
    title = expression('Spatial effect ('*W*')')
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3) +
  tm_text("NAME", size = 0.3)


tmap_save(tm = yearly.eta2_state, 
          filename = 'Analysis/Figures_BSFZINB/eta_statenames.eps', 
          device = cairo_ps, width = 8, height = 5.4)




# ------------------------------------------------
# income
# ------------------------------------------------
load('Analysis/refusal.RData')
refusal$fips = sapply(1:nrow(refusal), function(i) paste0(digits(refusal$fips[i], n = 5), collapse = ""))

refusalMap2012 = inner_join(counties_sf, refusal %>% filter(year == 2012))
refusalMap2013 = inner_join(counties_sf, refusal %>% filter(year == 2013))
refusalMap2014 = inner_join(counties_sf, refusal %>% filter(year == 2014))
refusalMap2015 = inner_join(counties_sf, refusal %>% filter(year == 2015))
refusalMap = rbind(refusalMap2012, refusalMap2013, refusalMap2014, refusalMap2015)


head(refusalMap)
unique(refusalMap$NAME)[grep('Anne', unique(refusalMap$NAME))]
refusalMap %>% filter(NAME == 'New York') %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80) %>% slice(c(1, 13, 25, 37))
refusalMap %>% filter(NAME == 'Anne Arundel') %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80) %>% slice(c(1, 13, 25, 37))

refusalMap %>% arrange(desc(pc_pop_highinc80)) %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80, NAME)
refusal %>% arrange(desc(pc_pop_highinc80)) %>% select(year, month, pc_pop_lowinc20, pc_pop_highinc80)

tmap_mode('plot')

plot_highinc = tm_shape(refusalMap) +
  tm_polygons(
    col = 'pc_pop_highinc80',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    title = 'Upper 80%'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'gray30', lwd = 0.3)

plot_highlow = tm_shape(refusalMap) +
  tm_polygons(
    col = 'pc_pop_lowinc20',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    title = 'Lower 20%'
  ) +
  tm_facets(by = 'year', ncol = 2) +
  tm_layout(legend.outside.size = 0.18,
            legend.title.size = 1.2,
            legend.text.size = 0.9) +
  tm_shape(states_sf) + tm_borders(col = 'gray30', lwd = 0.3)



