
library(ggplot2)
library(dplyr)
library(lubridate)

ldata = read.csv("length_data_fixed.csv")

ldata$year = ldata$date %>% mdy() %>% year()

rpor_data = ldata %>% subset(species == "Rhizoprionodon porosus")
sbra_data = ldata %>% subset(species == "Scomberomorus brasiliensis")


ggplot(rpor_data)+
  geom_histogram(aes(y = fl, fill = state))+
  scale_fill_manual(values = c('turquoise3', 'goldenrod1'))+
  facet_wrap(~year)+
  labs(y = 'FL (cm)', fill = 'State')+
  ggtitle('Rhizoprionodon porosus')+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = 'italic'))

ggplot(sbra_data)+
  geom_histogram(aes(y = fl, fill = state))+
  scale_fill_manual(values = c('turquoise3', 'goldenrod1'))+
  facet_wrap(~year)+
  labs(y = 'FL (cm)', fill = 'State')+
  ggtitle('Scomberomorus brasiliensis')+
  theme(axis.title.x = element_blank())+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = 'italic'))
