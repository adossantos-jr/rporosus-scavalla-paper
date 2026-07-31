
library(ggplot2)
library(dplyr)
library(lubridate)

ldata = read.csv("length_data_fixed.csv")

ldata$year = ldata$date %>% mdy() %>% year()

rpor_data = ldata %>% subset(species == "Rhizoprionodon porosus")
sbra_data = ldata %>% subset(species == "Scomberomorus brasiliensis")


ggplot(rpor_data)+
  geom_histogram(aes(y = fl, fill = state))+
  facet_wrap(~year)

ggplot(sbra_data)+
  geom_histogram(aes(y = fl, fill = state))+
  facet_wrap(~year)
