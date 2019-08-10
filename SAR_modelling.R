library(googlesheets)
library(raster)
library(rasterVis)
library(spdep)
library(dplyr)
library(maptools)
library(agricolae)
library(soiltexture)
library(plyr)
library(MASS)
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

# rm(list=ls())
##################################
#1 Clasificacion Textura del suelo
##################################

#1.1 Separar los tamanos de particula en objeto
textura <- na.omit(laspaf[,c(1,10:12)])

#1.2 Explorar triangulos texturales
TT.plot( class.sys = "USDA.TT",tri.data = textura,main="Clase Textura USDA",
         pch = 20)

TT.plot( class.sys = "HYPRES.TT",tri.data = textura,main="Clase Textura USDA",
         pch = 20)

#1.3 Clase textural por metodo European Soil Map
textura[,"EU"]<-TT.points.in.classes( #European Soil Map
  tri.data = textura,
  class.sys = "HYPRES.TT",
  PiC.type = "t",collapse = "y"
)
textura$EU <- as.factor(textura$EU)
levels(textura$EU) #para ver el orden
levels(textura$EU) <- c("Gruesa","Fina","Fina", "Media", "Medio Fina","Media",
                        "Medio Fina","Muy fina")
textura$EU <- as.character(textura$EU)

#1.4 Clase textural por USDA
textura[,"USDA"] <- TT.points.in.classes(
  tri.data = textura,
  class.sys = "USDA.TT",
  PiC.type = "t", collapse = ";"
)
textura$USDA <- as.factor(textura$USDA)
levels(textura$USDA) #para ver el orden
levels(textura$USDA) <- c("Arcilloso","Arcilloso","Arcilloso","Franco Arcilloso",
                          "Franco","Franco","Arena Franca", "Arena Franca",
                          "Arenoso","Arcillo Arenoso","Franco Arcillo Arenoso","Franco Arcillo Arenoso",
                          "Franco Arenoso", "Franco Arenoso","Limoso","Arcillo Limoso",
                          "Franco Arcillo Limoso","Franco Limoso")
textura$USDA <- as.character(textura$USDA)

#1.5 Agregar clase textural a laspaf
laspaf[textura$ID,"USDA"] <- textura$USDA
laspaf[textura$ID,"EU"] <- textura$EU


##################################
#2 Preparar la base de datos
##################################

#2.1 Crear SAR y PSI
laspaf$SAR <- (laspaf$Nas)/sqrt((laspaf$Cas+laspaf$Mgs)/2)
laspaf$PSI <- (laspaf$Na)*100/laspaf$CIC
laspaf$ESR <- (laspaf$Na)/(laspaf$CIC-laspaf$Na)

#2.2 Filtrar variables de interes
laspaf <- laspaf[,c(1:3,9,12:13,16,19,41:42,62:64)]

#2.3 Extraer suelos minerales
unique(laspaf$Txt)
laspaf <- laspaf[-which(laspaf$Txt=="Muestra materia organica"),]
laspaf <- laspaf[-which(laspaf$Txt=="Suelo Orgánico"),]
laspaf <- laspaf[-which(laspaf$Txt=="Material orgánico"),]
laspaf <- laspaf[,-6]

laspaf$USDA <- as.factor(laspaf$USDA)
laspaf$EU <- as.factor(laspaf$EU)

##################################
#3 Analisis Exploratorio de datos
##################################

#3.1 Resumen estadistico
summary(laspaf)
boxplot(laspaf$PSI)
boxplot(laspaf$SAR)
laspaf<- subset(laspaf,SAR<100 & ESR <= 1)

#laspaf1 <- subset(laspaf,SAR<12 & PSI <15)


#3.2 Correlacion entre variables
cortab <- laspaf[,c(4:9,12)]
cortab <- na.omit(cortab)
agricolae::correlation(cortab)


##################################
#3 Modelamiento Exploratorio
##################################
mod0 <- lm(PSI~.,na.action = na.omit,data=laspaf[,4:9,12])
summary(mod0)
step()

#3.1 Modelo lineal general
mod1 <- lm(SAR~PSI,na.action = na.omit,data=laspaf)
summary(mod1)
plot(x=laspaf$PSI,y=laspaf$SAR)
abline(mod1)
#Supuestos
plot(mod1)
nortest::ad.test(residuals(mod1))
car::ncvTest(mod1)
lmtest::dwtest(mod1,alternative="two.sided")

#3.2 Modelo lineal general + Clase textural
mod2 <- lm(SAR~PSI+EU,na.action = na.omit,data=laspaf)
summary(mod2)
#Supuestos
plot(mod2)
nortest::ad.test(residuals(mod2))
car::ncvTest(mod2)
lmtest::dwtest(mod2,alternative="two.sided")

atipico <- lala[c(3917,3927,4899),]
atipico <- laspaf[c(93),]

##################################
#4 Modelamiento GLM Gamma
##################################

#4.1 Modelo Gamma enlace identidad
mod3 <- glm(SAR~PSI,family=Gamma(link=identity),data=laspaf)
summary(mod3)

#4.2 Prueba de desvio del modelo
#Desvios altos indican falta de ajuste o modelo inadecuado
#El desvio debe ajustarse a una chi2
MASS::gamma.shape(mod3) #Alpha : parametro de precision
Dstar<-mod3$deviance*gamma.shape(mod3)$alpha 
1-pchisq(q = Dstar,df =  mod3$df.residual) #Ho: sigue una distribucion chi2

#4.3 Modelo gamma enlace log
mod4 <- glm(SAR~PSI,family=Gamma(link=log),data=laspaf)
summary(mod4)
Dstar<-mod4$deviance*gamma.shape(mod4)$alpha 
1-pchisq(q = Dstar,df =  mod4$df.residual)

#4.4 Diagnostico de modelo
car::residualPlots(mod3) #para homocedasticidad
car::influenceIndexPlot(mod3) #para valores influenciables
boot:: glm.diag.plots(mod3)

#Evaluar distancia de cook
car::influenceIndexPlot(mod3,vars="Cook")
plot(cooks.distance(mod3),type="h")
numero <- which(cooks.distance(mod3) > 0.01) #0.01 obtenido del grafico
atipico <- laspaf[c(numero),]
write.csv(atipico,file="D:/Documents/GitHub/LASPAF/atipicos.csv")
#Evaluar valores hat
car::influenceIndexPlot(mod3,vars="hat")
numero <- which(hatvalues(mod3) > 0.002) #0.01 obtenido del grafico
atipico <- laspaf[c(numero),]

#Evaluar AIC de los modelos
AIC(mod3)
AIC(mod4)


#source('http://www.poleto.com/funcoes/envel.gama.txt')
#source('http://www.poleto.com/funcoes/diag.gama.txt')
#diag.gama(mod3)
#envel.gama(mod3)

##################################
#Modelamiento Normal Inversa
##################################
mod5 <- glm(SAR~PSI,data=laspaf,
                 family=inverse.gaussian(link=identity))

mod5 <- glm(ESR~SAR,data=laspaf1,
            family=inverse.gaussian(link=identity))

##################################
#Modelamiento GLM doble
##################################
require(dglm)
mod5 <- dglm(SAR~PSI, ~ PSI,
            family=Gamma(link=identity),dlink = "log",
            data=laspaf1)
summary(mod5)
#AIC = -2*loglikelihood +2k || k: numero de parametos
2*(2+2) - mod5$m2loglik

mod6 <- dglm(SAR~PSI, ~ PSI,
             family=gaussian(link=identity),dlink = "log",
             data=laspaf1)
summary(mod6)
