setwd("...")

library(lme4)
library(lattice)
library(emmeans)
library(lattice)  #For fancy multipanel graphs
source("HighstatLibV10.R")
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(MASS)
library(visreg)
library(mgcv)
library(stringi)

theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

theme_set(theme_sleek())

################################# CPCE

cpce <- read.csv("shadingpercent2.csv")

str(cpce)

#library(rcompanion)
#cpce$tukeypercent <- transformTukey(cpce$percent)

hist(cpce$percent)
dotchart(cpce$percent)
boxplot(cpce$percent ~ cpce$time)
xyplot(cpce$percent ~ cpce$time | cpce$plot.id)

MyVar <- c("time", "plot.id")
pairs(cpce[,MyVar], 
      lower.panel = panel.cor)
corvif(cpce[,MyVar])

100 * sum(cpce$percent == 0) / nrow(cpce)

modelcpce <- lme(percent ~ time * plot.id, random = ~1|name/plot.area, data = cpce, method="REML")
modelcpce.fixed <- gls(percent ~ time * plot.id, data = cpce, method = "REML")

anova(modelcpce, modelcpce.fixed)

E1 <- resid(modelcpce, type = "pearson")
N  <- nrow(cpce)
p  <- length(coef(modelcpce))
sum(E1^2) / (N - p)

hist(resid(modelcpce))
plot(resid(modelcpce), fitted(modelcpce))
summary(modelcpce)

options(max.print=10000)
emmeans(modelcpce, list(pairwise ~ time * plot.id), adjust = "tukey")


cpce$time <- factor(cpce$time, levels=c("January","March", "April", "May", "June"), order=TRUE)

ggplot(cpce, aes(x=plot.id, y=percent)) +
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme_classic()

p<-ggplot(cpce, aes(x=plot.id, y=percent, fill=time)) +
  geom_boxplot()
p

cpce.bi <- read.csv("shadingcpcebinary.csv")

cpce.bi$time <- factor(cpce.bi$time, levels=c("January","March", "April", "May", "June"), order=TRUE)

hist(cpce.bi$binary)
dotchart(cpce.bi$binary)

modelcpce.bi <- glmer(binary ~ time * plot.id + (1|name), data = cpce.bi, family = binomial)
modelcpce.bi.fixed <- glm(binary ~ time * plot.id, data = cpce.bi, family = binomial)

anova(modelcpce.bi, modelcpce.bi.fixed)

hist(resid(modelcpce.bi))
hist(resid(modelcpce.bi.fixed))
plot(resid(modelcpce.bi), fitted(modelcpce.bi))
plot(resid(modelcpce.bi.fixed), fitted(modelcpce.bi.fixed))
summary(modelcpce.bi)

plot(modelcpce.bi)

visreg(modelcpce.bi, 'time', by = 'plot.id')


options(max.print=10000)
emmeans(modelcpce.bi, list(pairwise ~ time * plot.id), adjust = "tukey")

##################################### PAM

pam <- read.csv("shadingpam.csv")

str(pam)

plot(pam$yield)
hist(pam$yield)

dotchart(pam$yield)

MyVar <- c("time", "plot.id")
pairs(pam[,MyVar], 
      lower.panel = panel.cor)
corvif(pam[,MyVar])

100 * sum(pam$yield == 0) / nrow(pam)

modelpam <- lmer(yield ~ time * plot.id + (1|name), data = pam)

hist(resid(modelpam))
plot(resid(modelpam), fitted(modelpam))
summary(modelpam)

emmeans(modelpam, list(pairwise ~ time * plot.id), adjust = "tukey")

xyplot(yield ~ time | factor(plot.id), 
       data = pam,
       col = 1,
       pch = 16)

ggplot(pam, aes(x=time, y=yield)) +
  geom_boxplot()+
  theme_classic()

p<-ggplot(pam, aes(x=time, y=yield)) + geom_boxplot() + facet_grid(~plot.id)
p

######################### Bites

bites <- read.csv("shadingbiterates2.csv")
str(bites)
dotchart(bites$bites)
hist(bites$bites)
boxplot(bites$bitesm2 ~ bites$plot.id)
boxplot(bites$bitesm2 ~ bites$month)
bites$logbites <- log10(bites$bitesm2)
boxplot(bites$logbites ~ bites$month)
boxplot(bites$logbites ~ bites$plot.id)
boxplot(bites$logbites ~ bites$fg.fine)
hist(bites$logbites)


modelbites1 <- glmer(bitesm2 ~ month * plot.id *fg.fine + (1|name), family = "Gamma", data = bites)
modelbites2 <- lmer(logbites ~ month * plot.id *fg.fine + (1|name/plot.area), data = bites)

summary(modelbites1)
E1 <- resid(modelbites1)
plot(bites$month, E1)
plot(bites$plot.id, E1)
plot(bites$fg.fine, E1)

hist(resid(modelbites1))
plot(resid(modelbites1), fitted(modelbites1))
hist(resid(modelbites2))
plot(resid(modelbites2), fitted(modelbites2))


summary(modelbites2)
options(max.print=10000)
emmeans(modelbites2, pairwise ~ plot.id * month * fg.fine)

bites$month <- factor(bites$month, levels=c("March","April","May"), order=TRUE)
se<-function(x) sd(x)/sqrt(length(x))

temp <- bites %>% group_by(plot, month, plot.id, fg.fine) %>% summarise(bitesm2 = sum(bitesm2))
temp <- temp %>% group_by(month, plot.id, fg.fine) %>% summarise(bites = mean(bitesm2), se = se(bitesm2)) 

ggplot(temp, aes(x=month, y=bites)) +
  geom_bar(stat = "identity")+
  theme_classic()

temp$size.id<-stringr::str_split_fixed(temp$plot.id, '\\ ', 2)[,2]
temp$plot.id2<-stringr::str_split_fixed(temp$plot.id, '\\ ', 2)[,1]


p4 <- ggplot(data=temp, aes(x=month, y=bites, fill = plot.id2)) +
  geom_bar(position='dodge', stat="identity", width=0.75) +
  geom_errorbar(aes(ymin=bites-se, ymax=bites+se), position = position_dodge(0.75),width = 0.2) +
  facet_grid(~fg.fine) 
p4

##################### Tiles 

tiles <- read.csv("shadingtiles2.csv")

str(tiles)

dotchart(tiles$binary)
hist(tiles$binary)

MyVar <- c("pavement", "cca")
pairs(tiles[,MyVar], 
      lower.panel = panel.cor)
corvif(tiles[,MyVar])

modeltiles <- glmer(binary ~ plot.id + (1|name), family = binomial, data = tiles)

E1 <- resid(modeltiles, type = "pearson")
N  <- nrow(tiles)
p  <- length(coef(modeltiles))
sum(E1^2) / (N - p)

hist(resid(modeltiles))
plot(resid(modeltiles), fitted(modeltiles))
summary(modeltiles)

emmeans(modeltiles, pairwise ~ plot.id)

xyplot(turf ~ plot.id, 
       data = tiles,
       col = 1,
       pch = 16)

ggplot(tiles, aes(x=plot.id, y=turf)) +
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme_classic()

##################### Loggers

loggers <- read.csv("shadingloggers.csv")
str(loggers)

p<-ggplot(loggers, aes(x=time, y=lux, group=plot)) + geom_line(aes(linetype=plot))
p

############################## Final figures

library(wesanderson)

pal <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")
pal
cols<-c(pal[1], pal[2], pal[3], pal[4], pal[5])

p1<-ggplot(cpce, aes(x=time, y=percent)) +
  geom_boxplot(fill = cols[1]) + facet_grid(~plot.id) + 
  labs(x = '', y = 'Macroalgae cover (%)')
p1


p2<- ggplot(tiles, aes(x=plot.id, y=turf)) +
  geom_boxplot(fill=cols[1]) + 
  labs(x = '', y = 'Tiles: turf algae cover (%)')
p2


p3<-ggplot(pam, aes(x=time, y=yield)) + geom_boxplot(fill=cols[1]) + facet_grid(~plot.id) +
  labs(x = '', y = 'Photochemical quantum yield of photosystem II (Y(II))')
p3


p4 <- ggplot(data=temp, aes(x=month, y=bites, fill = plot.id2)) +
  geom_bar(position='dodge', stat="identity", width=0.75) +
  geom_errorbar(aes(ymin=bites-se, ymax=bites+se), position = position_dodge(0.75),width = 0.2) +
  facet_grid(size.id~fg.fine) + scale_color_manual(values = c(cols[1], cols[3])) +
  scale_fill_manual(name = "", values = c(cols[1], cols[3])) +
  labs(x = '', y = expression(paste("Bites ", m^-2)))
p4
