

data = read.csv('fishing_points.csv')


library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(sf)
library(ggspatial)
library(tidyterra)
library(terra)
library(raster)
library(marmap)

data$Date = data$Date %>% mdy() %>% format("%m-%Y")

data = data[order(as_date(data$Date, format = "%m-%Y")),] %>% as.data.frame()


states = ne_states(country = 'Brazil', return = 'sf')
pe = states %>% subset(name == 'Pernambuco')
rn = states %>% subset(name == 'Rio Grande do Norte')


data2 = data %>% count(Date, Estado)


data2 = data2[order(as_date(data2$Date, format = "%m-%Y")),] %>% as.data.frame()


lims_date = c("01-2010", "12-2022") %>% as_date(format = "%m-%Y")


bathy = getNOAA.bathy(
  lon1 = -39,
  lon2 = -33,
  lat1 = -10.5,
  lat2 = -3,
  resolution = 1,
  keep = F
) %>% as.xyz() %>%
subset(V3 > -60 & V3 < 0) %>%
  rasterFromXYZ() %>%
  rast()

map_br = 
ggplot()+
  geom_sf(data = states, fill = 'linen', color = 'grey60')+
  geom_sf(data = pe, fill = 'turquoise3')+
  geom_sf(data = rn, fill = 'goldenrod1')+
  theme_void()+
  geom_text(aes(x = -50, y = -13, label = 'Brazil'),
            size = 6)+
  annotate('rect', xmin = -37, xmax = -34, ymin = -9.5, ymax = -4,
           color = 'black', fill = 'transparent', size = 1)

 map = 
  ggplot()+
  geom_spatraster_contour_filled(data = bathy,
                                 breaks = seq(-60, 0, by = 10),
                                 alpha = .6)+
  geom_sf(data = states, fill = 'linen')+
  geom_sf(data = pe, fill = 'turquoise3')+
  geom_sf(data = rn, fill = 'goldenrod1')+
  geom_point(data = data,
             aes(x = Lon, y = Lat,
                 color = as_date(Date, format ="%m-%Y")),
             alpha = .7, size = 2)+
  scale_color_viridis_c(trans = 'date', option = 'magma',
                        guide = 'colorsteps')+
  scale_fill_manual(values = c('grey20', 'grey40', 'grey60', 'grey75', 'grey85', 'grey95'))+
  coord_sf(xlim = c(-37, -34), ylim = c(-9.5, -4))+
  guides(color = guide_colorsteps(barwidth = .5, barheight = 5),
         fill = guide_colorsteps(barwidth = .5, barheight = 5,
                                 show.limits = T))+
  annotation_scale(location = 'br')+
  annotation_north_arrow(location = 'tr', style = north_arrow_nautical)+
  geom_text(aes(x = -36.2, y = -5.8, label = 'RN'))+
  geom_text(aes(x = -35.9, y = -8.2, label = 'PE'))+
  scale_x_continuous(breaks = seq(-37, -34, by = 1))+
  labs(color = 'Date of monitoring', fill = 'Bathymetry (m)')+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_text(color = 'black'))

 maps = map_br + map + plot_layout(widths = c(.6, 1)) 
 maps = ggplotify::as.ggplot(maps)
 
plot = 
ggplot()+
  geom_col(data = data2, aes(x = as_date(Date, format = "%m-%Y"), 
          y = n, fill = Estado),
           alpha = .6, position = 'identity')+
  labs(y = "No. of monitored sets", fill = "State")+
  scale_x_date(breaks = '1 year', date_labels = '%Y')+
  scale_fill_manual(values = c('turquoise3', 'goldenrod1'))+
  theme_classic()+
  theme(axis.text = element_text(color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, color = 'black'))

 
maps / plot

ggsave('teste.png', dpi = 300)
