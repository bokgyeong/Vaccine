rm(list=ls())
library(sf); library(tmap); library(tidyverse); library(tigris); library(TeachingDemos); library(egg); library(grid)

counties_sf = counties(cb = TRUE, resolution = "20m")%>%
  shift_geometry()
counties_sf = counties_sf %>% unite('fips', STATEFP, COUNTYFP, sep = '')

states_sf <- states(cb = TRUE)


# ------------------------------------------------
# Map
# ------------------------------------------------
load(paste0('Analysis/real_ratio.RData'))

data.ysim$fips = sapply(1:nrow(data.ysim), function(i) paste0(digits(data.ysim$fips[i], n = 5), collapse = ""))

data.ysim.Map = data.frame()

for(i in unique(data.ysim$year)){
  # for(j in unique(data.ysim$month)){
  data.ysim.Map = rbind(data.ysim.Map, inner_join(counties_sf, data.ysim %>% filter(year == i)))
  # data.ysim.Map = rbind(data.ysim.Map, inner_join(counties_sf, data.ysim %>% filter(year == i, month == j)))
  # }
}



# mean / variance for positive counts
data.ysim.Map.disp = data.ysim.Map %>% filter(ysim > 0) %>% group_by(ID, fips) %>% summarise(dispersion = mean(ysim)/var(ysim))
data.ysim.Map.disp = data.ysim.Map.disp %>% filter(!is.infinite(dispersion))
data.ysim.Map.disp = data.ysim.Map.disp %>% na.omit()
# breaks.disp = c(round(quantile(data.ysim.Map.disp$dispersion,  probs = c(0, 0.3, 0.485, 0.6, 0.75, 0.85, 0.92, 0.96, 0.99141, 0.9979904)), 2), ceiling(max(data.ysim.Map.disp$dispersion)))
breaks.disp = c(0, 0.1, 0.3, 0.5, 1, 2.5, 4, 5, 10, 50, round(max(data.ysim.Map.disp$dispersion)))

ysim.disp = data.ysim.Map.disp %>% 
  tm_shape() +
  tm_polygons(
    col = 'dispersion',
    palette = 'Greens', 
    border.alpha = 0.2,
    breaks = breaks.disp,
    title = ''
  ) +
  tm_facets(by = 'ID') +
  tm_layout(legend.outside = T,
            legend.outside.size = 0.12,
            legend.text.size = 0.57,
            frame = FALSE,
            # main.title = '(b) Ratio of mean to variance for positive counts', 
            main.title.position = 'left', 
            main.title.size = 1.2) +
  tm_shape(states_sf) + tm_borders(col = 'black', lwd = 0.3)

tmap_save(tm = ysim.disp, filename = 'Analysis/ysim_ratio.eps', device = cairo_ps,
          width = 7.5, height = 6
          )


