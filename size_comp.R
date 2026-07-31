
library(ggplot2)
library(dplyr)
library(lubridate)

ldata = read.csv("length_data_fixed.csv")

ldata$year = dmy(ldata$date) %>% year()

ggplot(ldata)+
  geom_histogram(aes(y = as.numeric(fl)))




teste = data.frame(a = ldata$date, b = dmy(ldata$date))
