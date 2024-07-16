rm(list=ls())
options(bitmapType = 'cairo')
library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load('realnu/sum/refuse6m.RData')

claim = zicomp_refuse
claim$fips = sapply(1:nrow(claim), function(i) paste0(digits(claim$fips[i], n = 5), collapse = ""))

claimMap = data.frame()

for(i in unique(claim$year)){
  for(j in unique(claim$month)){
    claimMap = rbind(claimMap, inner_join(counties_sf, claim %>% filter(year == i, month == j)))
  }
}

claimMap = claimMap %>% mutate(yperpop = y / under5 * 1000, muperpop = mu / under5 * 1000)

# hist(claimMap$mu)
# hist(claimMap$y)
# hist(claimMap$mu / claimMap$under5 * 1000) 

breaks.both = unique(ceiling(quantile(c(claimMap$yperpop, claimMap$muperpop),  probs = c(0, 0.88, 0.96, 0.99885, 1))))

month = 6

## map of outcomes ----
obs.map = claimMap %>% 
  filter(month == month) %>% 
  tm_shape() +
  tm_polygons(
    col = 'yperpop',
    palette = 'YlOrRd',
    style = "cont",
    border.alpha = 0.2,
    breaks = breaks.both,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.23,
            legend.text.size = 0.58,
            frame = FALSE,
            main.title = '(a) Observed cases per 1,000 kids', 
            main.title.position = 'left', 
            main.title.size = 1.4) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# approximate mean
mu.map = claimMap %>% 
  filter(month == month) %>% 
  tm_shape() +
  tm_polygons(
    col = 'muperpop',
    palette = 'YlOrRd',
    style = "cont",
    border.alpha = 0.2,
    breaks = breaks.both,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.23,
            legend.text.size = 0.58,
            frame = FALSE,
            main.title = '(b) Expected cases per 1,000 kids', 
            main.title.position = 'left', 
            main.title.size = 1.4) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

# spatial dependence
library(RColorBrewer)
smu.map = claimMap %>% 
  filter(month == month) %>% 
  tm_shape() +
  tm_polygons(
    col = 'W',
    palette = 'seq',
    midpoint = 0, 
    style = "cont",
    border.alpha = 0.2,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.23,
            legend.text.size = 0.58,
            frame = FALSE,
            main.title = '(c) Spatial effects in refusal', 
            main.title.position = 'left', 
            aes.palette = list(seq = "-RdYlGn"),
            main.title.size = 1.4) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


# detection probability
pii.map = claimMap %>% 
  filter(month == month) %>% 
  tm_shape() +
  tm_polygons(
    col = 'pii',
    palette = 'Greens', 
    border.alpha = 0.2,
    style = "cont",
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.23,
            legend.text.size = 0.58,
            frame = FALSE,
            main.title = '(d) Detection probability', 
            main.title.position = 'left', 
            main.title.size = 1.4) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


spiimu.map = claimMap %>% 
  filter(month == month) %>% 
  tm_shape() +
  tm_polygons(
    col = 'V',
    midpoint = 0, 
    style = "cont",
    border.alpha = 0.2,
    title = ''
  ) +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.23,
            legend.text.size = 0.58,
            frame = FALSE,
            main.title = '(e) Spatial effects in detection', 
            main.title.position = 'left',
            main.title.size = 1.4) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)


plot.all = tmap_arrange(obs.map, mu.map, smu.map, pii.map, spiimu.map, ncol = 2)

tmap_save(tm = plot.all,
          filename = 'realnu/fig/refuseMaps.pdf',
          width = 8, height = 9)


