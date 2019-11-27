# Title:  Environment-phenotype associations
# Author: Van Wishingrad
# Date: 09 March, 2017
# Edited: 16 August, 2019
# Finalized: 27 November, 2019
# This script examines relationships between environmental variables and Sceloporus functional traits

# Script formatted for archival

# Load libraries
library(readr)
library(ggplot2)
library(nlme)
library(plyr)
library(gplots)
library(gstat)
library(sp)
require(MuMIn)

# SE function
se<-function(x) sqrt(var(x, na.rm=T)/length(x))

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

Sceloporus_specimens = read_delim("Sceloporus_specimens.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

# Select lizards >= 10 g only
Sceloporus_specimens_Age10g = Sceloporus_specimens[which(Sceloporus_specimens$Age10g=='Adult'),]

# Descriptive stats
# SVL
summary(Sceloporus_specimens_Age10g$SVL_mm)  
# Mass_g
summary(Sceloporus_specimens_Age10g$Mass_g)
# Dorsal scale counts
summary(Sceloporus_specimens_Age10g$Scales_dorsal)
# Dorsal color
summary(Sceloporus_specimens_Age10g$Color_Dorsal)

# data summary function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# Table 1
df3.SVL <- data_summary(Sceloporus_specimens_Age10g, varname="SVL_mm",
                    groupnames=c("Site"))
write.table(df3.SVL, "df3.SVL.txt⁩", sep="\t")
df3.dorsal.scales <- data_summary(Sceloporus_specimens_Age10g, varname="Scales_dorsal", 
                                   groupnames=c("Site"))
write.table(df3.dorsal.scales, "df3.dorsal.scales.txt⁩", sep="\t")
df3.dorsal.color <- data_summary(Sceloporus_specimens_Age10g, varname="Color_Dorsal", 
                                  groupnames=c("Site"))
write.table(df3.dorsal.color, "df3.dorsal.color.txt⁩", sep="\t")


# SVL
par(mfrow=c(1,1))
hist(Sceloporus_specimens_Age10g$SVL_mm, breaks = 20)

# Summer
# Mean Temp of Warmest Quarter

par(mfrow=c(1,1))
par(mar=c(5,5,2,2))

library(lme4)
library(lmerTest)
reg1d.lmer = lmer(SVL_mm ~ BIO10_meanTwarmestQ + (1|BIO10_meanTwarmestQ), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1d.lmer)
anova(reg1d.lmer)
plot(reg1d.lmer)
r.squaredGLMM(reg1d.lmer)

library(ggplot2)
reg1d.1 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO10_meanTwarmestQ, y = SVL_mm)) + 
  geom_point()
reg1d.1 + theme_classic() + labs(x = "Mean Temperature of Warmest Quarter (°C)", y = "Snout-Vent Length (mm)")


# Max Temp Warmest Month

library(lme4)
library(lmerTest)
reg1e.lmer = lmer(SVL_mm ~ BIO5_maxTwarmestMo + (1|BIO5_maxTwarmestMo), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1e.lmer)
anova(reg1e.lmer)
plot(reg1e.lmer)
r.squaredGLMM(reg1e.lmer)

library(ggplot2)
reg1d.1 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO5_maxTwarmestMo, y = SVL_mm)) + 
  geom_point()
reg1d.1 + theme_classic() + labs(x = "Maximum Temperature of Warmest Month (°C)", y = "Snout-Vent Length (mm)")


# Mean Temp of Coldest Quarter

library(lme4)
library(lmerTest)
reg1b.lmer = lmer(SVL_mm ~ BIO11_meanTcoldestQ + (1|BIO11_meanTcoldestQ), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1b.lmer)
anova(reg1b.lmer)
plot(reg1b.lmer)
r.squaredGLMM(reg1b.lmer)

library(ggplot2)
reg1d.1 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO11_meanTcoldestQ, y = SVL_mm)) + 
  geom_point()
reg1d.1 + theme_classic() + labs(x = "Mean Temperature of Coldest Quarter (°C)", y = "Snout-Vent Length (mm)")

# Min temp coldest month

library(lme4)
library(lmerTest)
reg1c.lmer = lmer(SVL_mm ~ BIO6_minTcoldestMo + (1|BIO6_minTcoldestMo), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1c.lmer)
anova(reg1c.lmer)
plot(reg1c.lmer)
r.squaredGLMM(reg1c.lmer)

library(ggplot2)
reg1c.1 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO6_minTcoldestMo, y = SVL_mm)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
reg1c.1 + theme_classic() + labs(x = "Minimum Temperature of Coldest Month (°C)", y = "Snout-Vent Length (mm)")

reg1c <- lm(Sceloporus_specimens_Age10g$SVL_mm~Sceloporus_specimens_Age10g$BIO6_minTcoldestMo)
summary(reg1c)


## Spatially-explicit models

# note: some values my differ slightly from publication values due to variation in random jitter function of gps points

# Load dataset
library(readr)
library(ggplot2)
library(nlme)
library(plyr)
library(gplots)

# SE function
se<-function(x) sqrt(var(x, na.rm=T)/length(x))

Sceloporus_VertNet = read_delim("Sceloporus_occidentalis_VertNet.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

summary(Sceloporus_VertNet)

# Adjust temperature scale
Sceloporus_VertNet$bio_5 = Sceloporus_VertNet$bio_5/10
Sceloporus_VertNet$bio_6 = Sceloporus_VertNet$bio_6/10
Sceloporus_VertNet$bio_10 = Sceloporus_VertNet$bio_10/10
Sceloporus_VertNet$bio_11 = Sceloporus_VertNet$bio_11/10

# Select lizards 64 g and larger only
Sceloporus_VertNet_65mm = Sceloporus_VertNet[which(Sceloporus_VertNet$SVL_mm>63.5),]

summary(Sceloporus_VertNet_65mm)
hist(Sceloporus_VertNet_65mm$SVL_mm)

par(mar=c(5,5,2,2))

# Mean Temperature of Warmest Quarter (bio10)

library(rgdal)
xy = cbind(Sceloporus_VertNet_65mm$decimallongitude, Sceloporus_VertNet_65mm$decimallatitude)
utms = project(xy, "+proj=utm +zone=11 ellps=WGS84")
Sceloporus_VertNet_65mm$northing = utms[,1]/1000
Sceloporus_VertNet_65mm$easting = utms[,2]/1000

library(sp)
Sceloporus_VertNet_65mm.sp = Sceloporus_VertNet_65mm
coordinates(Sceloporus_VertNet_65mm.sp) = c('easting', 'northing')
Sceloporus_VertNet_65mm.sp$SVL_mm.dev = Sceloporus_VertNet_65mm.sp$SVL_mm - mean(Sceloporus_VertNet_65mm.sp$SVL_mm)
bubble(Sceloporus_VertNet_65mm.sp, zcol = 'SVL_mm.dev', scales = list(draw = T))

library(geoR)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data = Sceloporus_VertNet_65mm$SVL_mm)
par(mfrow=c(1,1))
plot(v1)

library(ncf)
ncf.scor <- spline.correlog(Sceloporus_VertNet_65mm$easting, Sceloporus_VertNet_65mm$northing, Sceloporus_VertNet_65mm$SVL_mm,
                            resamp=10, quiet = FALSE)
plot(ncf.scor)

mod.nocorr = lm(SVL_mm ~ bio_10, data = Sceloporus_VertNet_65mm)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data =resid(mod.nocorr))
plot(v1)

# jitter GPS points 
Sceloporus_VertNet_65mm$j.easting = jitter(Sceloporus_VertNet_65mm$easting)
Sceloporus_VertNet_65mm$j.northing = jitter(Sceloporus_VertNet_65mm$northing)

## corExp
SVL_mm.corExp = gls(SVL_mm ~ bio_10, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corExp(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corExp)                            

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

## corGaus
SVL_mm.corGaus = gls(SVL_mm ~ bio_10, 
                     data = Sceloporus_VertNet_65mm, 
                     correlation = corGaus(form = ~j.easting + j.northing, 
                                           nugget = TRUE, 
                                           value = c(50, 0.5)))
summary(SVL_mm.corGaus)      

SVL_mm.corGaus.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corGaus, type = "normalized"))
plot(SVL_mm.corGaus.var)

# corLin
SVL_mm.corLin = gls(SVL_mm ~bio_10, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corLin(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corLin)                            

SVL_mm.corLin.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corLin, type = "normalized"))
plot(SVL_mm.corLin.var)

# corRatio
SVL_mm.corRatio = gls(SVL_mm ~ bio_10, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corRatio(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corRatio)                            

SVL_mm.corRatio.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corRatio, type = "normalized"))
plot(SVL_mm.corRatio.var)

# corSpher
SVL_mm.corSpher = gls(SVL_mm ~ bio_10, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corSpher(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corSpher)                            

SVL_mm.corSpher.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corSpher, type = "normalized"))
plot(SVL_mm.corSpher.var)

# calculate AIC to select best model
AIC(SVL_mm.corExp, SVL_mm.corGaus, SVL_mm.corLin, SVL_mm.corRatio, SVL_mm.corSpher)
# SVL_mm.corExp is best
# SVL_mm.corSpher & SVL_mm.corLin are also very close

# best model
summary(SVL_mm.corExp)

# other good models
summary(SVL_mm.corSpher)
summary(SVL_mm.corLin)

# all give same interpretation

library(ggplot2)
plot3 = ggplot(Sceloporus_VertNet_65mm, aes(x = bio_10, y = SVL_mm)) + 
  geom_point()
plot3 + theme_classic() + labs(x = "Mean Temperature of Warmest Quarter", y = "Snout-Vent Length (mm)")

reg3 <- lm(Sceloporus_VertNet_65mm$SVL_mm~Sceloporus_VertNet_65mm$bio_10)
summary(reg3)

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

#plot

library(ggplot2)
plot3.1 = ggplot(Sceloporus_VertNet_65mm, aes(x = bio_10, y = SVL_mm)) + 
  geom_point()
plot3.1 + theme_classic() + labs(x = "Mean Temperature of Warmest Quarter", y = "Snout-Vent Length (mm)")




### Max Temp of Warmest Month (bio_5)

library(rgdal)
xy = cbind(Sceloporus_VertNet_65mm$decimallongitude, Sceloporus_VertNet_65mm$decimallatitude)
utms = project(xy, "+proj=utm +zone=11 ellps=WGS84")
Sceloporus_VertNet_65mm$northing = utms[,1]/1000
Sceloporus_VertNet_65mm$easting = utms[,2]/1000

library(sp)
Sceloporus_VertNet_65mm.sp = Sceloporus_VertNet_65mm
coordinates(Sceloporus_VertNet_65mm.sp) = c('easting', 'northing')
Sceloporus_VertNet_65mm.sp$SVL_mm.dev = Sceloporus_VertNet_65mm.sp$SVL_mm - mean(Sceloporus_VertNet_65mm.sp$SVL_mm)
bubble(Sceloporus_VertNet_65mm.sp, zcol = 'SVL_mm.dev', scales = list(draw = T))

library(geoR)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data = Sceloporus_VertNet_65mm$SVL_mm)
par(mfrow=c(1,1))
plot(v1)

library(ncf)
ncf.scor <- spline.correlog(Sceloporus_VertNet_65mm$easting, Sceloporus_VertNet_65mm$northing, Sceloporus_VertNet_65mm$SVL_mm,
                            resamp=10, quiet = FALSE)
plot(ncf.scor)

mod.nocorr = lm(SVL_mm ~ bio_5, data = Sceloporus_VertNet_65mm)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data =resid(mod.nocorr))
plot(v1)

# jitter GPS points 
Sceloporus_VertNet_65mm$j.easting = jitter(Sceloporus_VertNet_65mm$easting)
Sceloporus_VertNet_65mm$j.northing = jitter(Sceloporus_VertNet_65mm$northing)

## corExp
SVL_mm.corExp = gls(SVL_mm ~ bio_5, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corExp(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corExp)                            

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

## corGaus
SVL_mm.corGaus = gls(SVL_mm ~ bio_5, 
                     data = Sceloporus_VertNet_65mm, 
                     correlation = corGaus(form = ~j.easting + j.northing, 
                                           nugget = TRUE, 
                                           value = c(50, 0.5)))
summary(SVL_mm.corGaus)      

SVL_mm.corGaus.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corGaus, type = "normalized"))
plot(SVL_mm.corGaus.var)

# corLin
SVL_mm.corLin = gls(SVL_mm ~bio_5, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corLin(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corLin)                            

SVL_mm.corLin.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corLin, type = "normalized"))
plot(SVL_mm.corLin.var)

# corRatio
SVL_mm.corRatio = gls(SVL_mm ~ bio_5, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corRatio(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corRatio)                            

SVL_mm.corRatio.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corRatio, type = "normalized"))
plot(SVL_mm.corRatio.var)

# corSpher
SVL_mm.corSpher = gls(SVL_mm ~ bio_5, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corSpher(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corSpher)                            

SVL_mm.corSpher.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corSpher, type = "normalized"))
plot(SVL_mm.corSpher.var)

# calculate AIC to select best model
AIC(SVL_mm.corExp, SVL_mm.corGaus, SVL_mm.corLin, SVL_mm.corRatio, SVL_mm.corSpher)
# SVL_mm.corExp

# best model
summary(SVL_mm.corExp)

# other good models
summary(SVL_mm.corLin)
summary(SVL_mm.corSpher)
# same interpretation

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

library(ggplot2)
plot4 = ggplot(Sceloporus_VertNet_65mm, aes(x = bio_5, y = SVL_mm)) + 
  geom_point()
plot4 + theme_classic() + labs(x = "Maximum Temperature of Warmest Month", y = "Snout-Vent Length (mm)")

reg4 <- lm(Sceloporus_VertNet_65mm$SVL_mm~Sceloporus_VertNet_65mm$bio_5)
summary(reg4)


# Winter
# Mean Temp of Coldest Quarter (bio_11)

library(rgdal)
xy = cbind(Sceloporus_VertNet_65mm$decimallongitude, Sceloporus_VertNet_65mm$decimallatitude)
utms = project(xy, "+proj=utm +zone=11 ellps=WGS84")
Sceloporus_VertNet_65mm$northing = utms[,1]/1000
Sceloporus_VertNet_65mm$easting = utms[,2]/1000

library(sp)
Sceloporus_VertNet_65mm.sp = Sceloporus_VertNet_65mm
coordinates(Sceloporus_VertNet_65mm.sp) = c('easting', 'northing')
Sceloporus_VertNet_65mm.sp$SVL_mm.dev = Sceloporus_VertNet_65mm.sp$SVL_mm - mean(Sceloporus_VertNet_65mm.sp$SVL_mm)
bubble(Sceloporus_VertNet_65mm.sp, zcol = 'SVL_mm.dev', scales = list(draw = T))

library(geoR)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data = Sceloporus_VertNet_65mm$SVL_mm)
par(mfrow=c(1,1))
plot(v1)

library(ncf)
ncf.scor <- spline.correlog(Sceloporus_VertNet_65mm$easting, Sceloporus_VertNet_65mm$northing, Sceloporus_VertNet_65mm$SVL_mm,
                            resamp=10, quiet = FALSE)
plot(ncf.scor)

mod.nocorr = lm(SVL_mm ~ bio_11, data = Sceloporus_VertNet_65mm)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data =resid(mod.nocorr))
plot(v1)

# jitter GPS points 
Sceloporus_VertNet_65mm$j.easting = jitter(Sceloporus_VertNet_65mm$easting)
Sceloporus_VertNet_65mm$j.northing = jitter(Sceloporus_VertNet_65mm$northing)

## corExp
SVL_mm.corExp = gls(SVL_mm ~ bio_11, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corExp(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corExp)                            

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

## corGaus
SVL_mm.corGaus = gls(SVL_mm ~ bio_11, 
                     data = Sceloporus_VertNet_65mm, 
                     correlation = corGaus(form = ~j.easting + j.northing, 
                                           nugget = TRUE, 
                                           value = c(50, 0.5)))
summary(SVL_mm.corGaus)      

SVL_mm.corGaus.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corGaus, type = "normalized"))
plot(SVL_mm.corGaus.var)

# corLin
SVL_mm.corLin = gls(SVL_mm ~bio_11, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corLin(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corLin)                            

SVL_mm.corLin.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corLin, type = "normalized"))
plot(SVL_mm.corLin.var)

# corRatio
SVL_mm.corRatio = gls(SVL_mm ~ bio_11, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corRatio(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corRatio)                            

SVL_mm.corRatio.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corRatio, type = "normalized"))
plot(SVL_mm.corRatio.var)

# corSpher
SVL_mm.corSpher = gls(SVL_mm ~ bio_11, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corSpher(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corSpher)                            

SVL_mm.corSpher.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corSpher, type = "normalized"))
plot(SVL_mm.corSpher.var)

# calculate AIC to select best model
AIC(SVL_mm.corExp, SVL_mm.corGaus, SVL_mm.corLin, SVL_mm.corRatio, SVL_mm.corSpher)
#SVL_mm.corLin

summary(SVL_mm.corLin)
#best support

# also look at these with approx. equal support
summary(SVL_mm.corExp)
summary(SVL_mm.corSpher)
# same result

reg1 <- lm(Sceloporus_VertNet_65mm$SVL_mm~Sceloporus_VertNet_65mm$bio_11)
summary(reg1)
abline(reg1, col = "black")

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)


#plot
library(ggplot2)
plot1 = ggplot(Sceloporus_VertNet_65mm, aes(x = bio_11, y = SVL_mm)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
plot1 + theme_classic() + labs(x = "Mean Temperature of Coldest Quarter", y = "Snout-Vent Length (mm)")






### Min Temp of Coldest Month (bio_6)

library(rgdal)
xy = cbind(Sceloporus_VertNet_65mm$decimallongitude, Sceloporus_VertNet_65mm$decimallatitude)
utms = project(xy, "+proj=utm +zone=11 ellps=WGS84")
Sceloporus_VertNet_65mm$northing = utms[,1]/1000
Sceloporus_VertNet_65mm$easting = utms[,2]/1000

library(sp)
Sceloporus_VertNet_65mm.sp = Sceloporus_VertNet_65mm
coordinates(Sceloporus_VertNet_65mm.sp) = c('easting', 'northing')
Sceloporus_VertNet_65mm.sp$SVL_mm.dev = Sceloporus_VertNet_65mm.sp$SVL_mm - mean(Sceloporus_VertNet_65mm.sp$SVL_mm)
bubble(Sceloporus_VertNet_65mm.sp, zcol = 'SVL_mm.dev', scales = list(draw = T))

library(geoR)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data = Sceloporus_VertNet_65mm$SVL_mm)
par(mfrow=c(1,1))
plot(v1)

library(ncf)
ncf.scor <- spline.correlog(Sceloporus_VertNet_65mm$easting, Sceloporus_VertNet_65mm$northing, Sceloporus_VertNet_65mm$SVL_mm,
                            resamp=10, quiet = FALSE)
plot(ncf.scor)

mod.nocorr = lm(SVL_mm ~ bio_6, data = Sceloporus_VertNet_65mm)
v1 <- variog(coords = Sceloporus_VertNet_65mm[,c('easting', 'northing')], data =resid(mod.nocorr))
plot(v1)

# jitter GPS points 
Sceloporus_VertNet_65mm$j.easting = jitter(Sceloporus_VertNet_65mm$easting)
Sceloporus_VertNet_65mm$j.northing = jitter(Sceloporus_VertNet_65mm$northing)

## corExp
SVL_mm.corExp = gls(SVL_mm ~ bio_6, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corExp(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corExp)                            

#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

## corGaus
SVL_mm.corGaus = gls(SVL_mm ~ bio_6, 
                     data = Sceloporus_VertNet_65mm, 
                     correlation = corGaus(form = ~j.easting + j.northing, 
                                           nugget = TRUE, 
                                           value = c(50, 0.5)))
summary(SVL_mm.corGaus)      

SVL_mm.corGaus.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corGaus, type = "normalized"))
plot(SVL_mm.corGaus.var)

# corLin
SVL_mm.corLin = gls(SVL_mm ~bio_6, 
                    data = Sceloporus_VertNet_65mm, 
                    correlation = corLin(form = ~j.easting + j.northing, 
                                         nugget = TRUE, 
                                         value = c(50, 0.5)))
summary(SVL_mm.corLin)                            

SVL_mm.corLin.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corLin, type = "normalized"))
plot(SVL_mm.corLin.var)

# corRatio
SVL_mm.corRatio = gls(SVL_mm ~ bio_6, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corRatio(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corRatio)                            

SVL_mm.corRatio.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corRatio, type = "normalized"))
plot(SVL_mm.corRatio.var)

# corSpher
SVL_mm.corSpher = gls(SVL_mm ~ bio_6, 
                      data = Sceloporus_VertNet_65mm, 
                      correlation = corSpher(form = ~j.easting + j.northing, 
                                             nugget = TRUE, 
                                             value = c(50, 0.5)))
summary(SVL_mm.corSpher)                            

SVL_mm.corSpher.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corSpher, type = "normalized"))
plot(SVL_mm.corSpher.var)

# calculate AIC to select best model
AIC(SVL_mm.corExp, SVL_mm.corGaus, SVL_mm.corLin, SVL_mm.corRatio, SVL_mm.corSpher)
# SVL_mm.corLin

# best support
summary(SVL_mm.corLin)

# also good
summary(SVL_mm.corExp)
summary(SVL_mm.corSpher)
# same result


reg2 <- lm(Sceloporus_VertNet_65mm$SVL_mm~Sceloporus_VertNet_65mm$bio_6)
summary(reg2)


#variogram
SVL_mm.corExp.var <- variog(coords = Sceloporus_VertNet_65mm[,c('j.easting', 'j.northing')], data = resid(SVL_mm.corExp, type = "normalized"))
plot(SVL_mm.corExp.var)

# plot
library(ggplot2)
plot2 = ggplot(Sceloporus_VertNet_65mm, aes(x = bio_6, y = SVL_mm)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
plot2 + theme_classic() + labs(x = "Minimum Temperature of Coldest Month", y = "Snout-Vent Length (mm)")








# Scales

# Load libraries
library(readr)
library(ggplot2)
library(nlme)
library(plyr)
library(gplots)
library(gstat)
library(sp)

# SE function
se<-function(x) sqrt(var(x, na.rm=T)/length(x))

Sceloporus_specimens = read_delim("Sceloporus_specimens.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
summary(Sceloporus_specimens)

# Select lizards > 10 g only
Sceloporus_specimens_Age10g = Sceloporus_specimens[which(Sceloporus_specimens$Age10g=='Adult'),]

# select lizards with scale data only
df.complete = subset(Sceloporus_specimens_Age10g, !(is.na(Scales_dorsal)))

plot(df.complete$Scales_dorsal, df.complete$SVL_mm)
SVL.Scales = lm(df.complete$SVL_mm~df.complete$Scales_dorsal)
abline(SVL.Scales, col = "red")
SVL.Scales.resid = resid(SVL.Scales)
SVL.Scales.predict = predict(SVL.Scales)
length(resid(SVL.Scales))
length(predict(SVL.Scales))

par(mfrow=c(1,1))
plot(predict(SVL.Scales), resid(SVL.Scales))
abline(0,0)
# POSITIVE residuals are larger scales, NEGATIVE resuduals are smaller scales

resid.Scales_dorsal = resid(SVL.Scales)
df.complete$resid.Scales_dorsal = resid.Scales_dorsal

# Correlations among scale traits

df.complete$Scales_size = (df.complete$SVL_mm/df.complete$Scales_dorsal)

df.complete$Scale.count = df.complete$Scales_dorsal
df.complete$Scale.size = df.complete$Scales_size

# size vs number of scales
cor.test(df.complete$SVL_mm, df.complete$Scales_dorsal)

# size vs mean scale size
cor.test(df.complete$SVL_mm, df.complete$Scale.size)


### SIMPLE CORRELATION FIGURES ###

library(ggplot2)
modelR1 = ggplot(df.complete, aes(x = SVL_mm, y = Scales_dorsal)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
modelR1 + theme_classic() + labs(x = "Snout-Vent Length (mm)", y = "Number of Dorsal Scales")

library(ggplot2)
modelR2 = ggplot(df.complete, aes(x = SVL_mm, y = Scale.size)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
modelR2 + theme_classic() + labs(x = "Snout-Vent Length (mm)", y = "Mean Dorsal Scale Size (mm)")



# Mean Temperature of Warmest Quarter

reg3d2.lmer.BIO10 = lmer(resid.Scales_dorsal ~ BIO10_meanTwarmestQ + (1|BIO10_meanTwarmestQ), REML = F, na.action=na.omit, data=df.complete)
#summary(reg3d2.lmer.BIO10)
anova(reg3d2.lmer.BIO10)
plot(reg3d2.lmer.BIO10)
r.squaredGLMM(reg3d2.lmer.BIO10)

library(ggplot2)
plot2 = ggplot(df.complete, aes(x = BIO10_meanTwarmestQ, y = SVL_mm)) + 
  geom_point()
plot2 + theme_classic() + labs(x = "Mean Temperature of Warmest Quarter (°C)", y = "Snout-Vent Length (mm)")


# Maximum Temperature of Warmest Month

reg3d2.lmer.BIO5 = lmer(resid.Scales_dorsal ~ BIO5_maxTwarmestMo + (1|BIO5_maxTwarmestMo), REML = F, na.action=na.omit, data=df.complete)
#summary(reg3d2.lmer.BIO5)
anova(reg3d2.lmer.BIO5)
plot(reg3d2.lmer.BIO5)
r.squaredGLMM(reg3d2.lmer.BIO5)

library(ggplot2)
plot2 = ggplot(df.complete, aes(x = BIO5_maxTwarmestMo, y = SVL_mm)) + 
  geom_point()
plot2 + theme_classic() + labs(x = "Maximum Temperature of Warmest Month (°C)", y = "Snout-Vent Length (mm)")


# Aridity

reg3d2.lmer = lmer(resid.Scales_dorsal ~ Aridity_index + (1|Aridity_index), REML = F, na.action=na.omit, data=df.complete)
summary(reg3d2.lmer)
anova(reg3d2.lmer)
r.squaredGLMM(reg3d2.lmer)

library(ggplot2)
reg3d2 = ggplot(df.complete, aes(x = Aridity_index, y = resid.Scales_dorsal)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")
reg3d2 + theme_classic() + labs(x = "Aridity Index", y = "Dorsal Scale Size (resid.)")






###############################################################################
#################################################################################

# Dorsal color
par(mfrow=c(1,1))
hist(Sceloporus_specimens_Age10g$Color_Dorsal, breaks = 20)

# Summer
# Mean Temp of Warmest Quarter

par(mfrow=c(1,1))
par(mar=c(5,5,2,2))

library(lme4)
library(lmerTest)
reg1d2.lmer = lmer(Color_Dorsal ~ BIO10_meanTwarmestQ + (1|BIO10_meanTwarmestQ), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1d2.lmer)
anova(reg1d2.lmer)
plot(reg1d2.lmer)
r.squaredGLMM(reg1d2.lmer)

library(ggplot2)
reg1d2 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO10_meanTwarmestQ, y = Color_Dorsal)) + 
  geom_point()
reg1d2 + theme_classic() + labs(x = "Mean Temperature of Warmest Quarter (°C)", y = "Dorsal Lightness")

# Max Temp Warmest Month

library(ggplot2)
reg1e2 = ggplot(Sceloporus_specimens_Age10g, aes(x = BIO5_maxTwarmestMo, y = Color_Dorsal)) + 
  geom_point()
reg1e2 + theme_classic() + labs(x = "Maximum Temperature of Warmest Month (°C)", y = "Dorsal Lightness")

library(lme4)
library(lmerTest)
reg1e2.lmer = lmer(Color_Dorsal ~ BIO5_maxTwarmestMo + (1|BIO5_maxTwarmestMo), REML = F, na.action=na.omit, data=Sceloporus_specimens_Age10g)
summary(reg1e2.lmer)
anova(reg1e2.lmer)
plot(reg1e2.lmer)
r.squaredGLMM(reg1e2.lmer)

# end

#
#
#