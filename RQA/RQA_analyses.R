# wipe work space
rm(list = ls())

# set directory
setwd("/Users/wallot/Desktop/laurens stuff/")

# load data
data <- read.csv('forSeb.csv')
music <- data$music
music <- music[seq(1,(length(music)),by = 3)]
plot(music, type = "l")
avg_pupil <- data$avg_pupil
avg_pupil <- avg_pupil[seq(1,(length(avg_pupil)),by = 3)]
plot(avg_pupil, type = "l")

# load packages
library(crqa)

# plot time series
plot(avg_pupil, type = "l")
plot(music, type = "l")

del <- 10
emb <- 2
rad <- 0.3

pup <- crqa(ts1 = avg_pupil, ts2 = avg_pupil, delay = del, embed = emb, rescale = 0, radius = rad, normalize = 2, method = "rqa", tw = 0)

parC <- list(unit = 200, labelx = "pupil trace", labely = "pupil trace",
             cols = "black", pcex = .2 , pch = 15, las = 0,
             labax = seq(0, 1000, 250),
             labay = seq(0, 1000, 250))

plotRP(pup$RP,parC)

mus <- crqa(ts1 = music, ts2 = music, delay = del, embed = emb, rescale = 0, radius = rad, normalize = 2, method = "rqa", tw = 0)

parC <- list(unit = 200, labelx = "amplitude envelope", labely = "amplitude envelope",
             cols = "black", pcex = .2 , pch = 15, las = 0,
             labax = seq(0, 1000, 250),
             labay = seq(0, 1000, 250))

plotRP(mus$RP,parC)


pup_mus <- crqa(ts1 = avg_pupil, ts2 = music, delay = del, embed = emb, rescale = 0, radius = rad, normalize = 2, method = "crqa", tw = 0)

parC <- list(unit = 200, labelx = "pupil trace", labely = "amplitude envelope",
             cols = "black", pcex = .2 , pch = 15, las = 0,
             labax = seq(0, 1000, 250),
             labay = seq(0, 1000, 250))

plotRP(pup_mus$RP,parC)

win = 10
prof <- drpfromts(ts1 = avg_pupil, ts2 = music, windowsize = win, radius = rad,
          delay = del, embed = emb, normalize = 2, datatype = "continuous")

plot(seq(-win,win,by = 1), prof$profile, type = "o",
     xlab = "lag",
     ylab = "recurrence rate")

