library(googlesheets)
library(raster)
library(rasterVis)
library(spdep)
library(dplyr)
library(maptools)
library(agricolae)
##################################
#Lectura de datos con googlesheets
##################################

gs_auth()
my_sheets <- gs_ls()
gap <- gs_title("LASPAF_Salinidad")
laspaf <- gap %>%
  gs_read(ws = "Datos")
laspaf <- as.data.frame(laspaf)
rm(my_sheets,gap)
laspaf <- laspaf[,c(1:2,8:41)]

laspaf$SAR <- (laspaf$Nas)/sqrt((laspaf$Cas+laspaf$Mgs)/2)
laspaf$PSI <- (laspaf$Na)*100/laspaf$CIC

table(laspaf$Txt)
##################################
#Modelamiento Exploratorio
##################################
str(laspaf)

#dt$Age[dt$Age == 99] <- NA para indicar valores NA

mod1 <- lm(SAR~PSI,na.action = na.omit,data=laspaf)
summary(mod1)
plot(mod1)
plot(x=laspaf$PSI,y=laspaf$SAR)

car::ncvTest(mod1)
laspaf1 <- subset(laspaf,SAR<40)
mod2 <- lm(SAR~PSI,na.action = na.omit,data=laspaf1)
plot(mod2)
car::ncvTest(mod2)

atipico <- lala[c(3917,3927,4899),]
atipico <- laspaf[c(93),]

##################################
#Modelamiento GLM Gamma
##################################
mod3 <- glm(SAR~PSI,family=Gamma(link=identity),data=laspaf1)
summary(mod3)
plot(mod3)
boot:: glm.diag.plots(mod3)

car::residualPlots(mod3)
car::influenceIndexPlot(mod3)

# Utilizar la distancia de cook
numero <- which(cooks.distance(mod3) > 0.01)
plot(cooks.distance(mod3),type="h")

atipico <- laspaf1[c(numero),]


mod4 <- glm(SAR~PSI,family=Gamma(link=log),data=laspaf1)
summary(mod4)
boot:: glm.diag.plots(mod4)
car::residualPlots(mod4)
car::influenceIndexPlot(mod4)

#source('http://www.poleto.com/funcoes/envel.gama.txt')
#source('http://www.poleto.com/funcoes/diag.gama.txt')
#diag.gama(mod3)
#envel.gama(mod3)

##################################
#Modelamiento GLM doble
##################################
require(dglm)
mod5 = dglm(SAR~PSI, ~ PSI,
            family=Gamma(link=identity),dlink = "log",
            data=laspaf1)
summary(mod5)
#AIC = -2*loglikelihood +2k || k: numero de parametos
2*(2+2) - mod5$m2loglik