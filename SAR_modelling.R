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

laspaf$SAR <- (laspaf$`Na+`)/sqrt((laspaf$`Ca2+`+laspaf$`Mg2+`)/2)
laspaf$PSI <- (laspaf$Na)*100/laspaf$CIC

##################################
#Modelamiento Exploratorio
##################################
str(laspaf)

#dt$Age[dt$Age == 99] <- NA para indicar valores NA

mod1 <- lm(laspaf$SAR~laspaf$PSI,na.action = na.omit)
summary(mod1)

plot(x=laspaf$PSI,y=laspaf$SAR)


laspaf1 <- subset(laspaf,PSI > 80)
plot(x=lala$PSI,y=lala$SAR)

lala<-subset(laspaf,SAR < 40)
hist(lala$SAR)
mod2 <- lm(lala$SAR~lala$PSI,na.action = na.omit)
summary(mod2)
atipico <- lala[c(3917,3927,4899),]
atipico <- laspaf[c(93),]

mod3 <- glm(lala$SAR~lala$PSI,family=Gamma(link=identity))
summary(mod3)
plot(mod3)
source('http://www.poleto.com/funcoes/envel.gama.txt')
envel.gama(mod3)
