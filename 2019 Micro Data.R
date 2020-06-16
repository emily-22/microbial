microdf <- read.csv("Emily Data.csv")

spongedf <- subset(microdf, Substrate=="Cellulose Sponge" & Week.Removed==4)
wooddf <- subset(microdf, Substrate=="Wood Veneer" & Week.Removed==6)

library(lme4)
library(ggplot2)
library(ggmap)
library(AICcmodavg)
library(scatterplot3d)
library(devtools)
library(cowplot)
library(MuMIn)
library(usdm)


cs1 <- data.frame(scale(microdf$Conductivity),
                  scale(microdf$H),
                  scale(microdf$DO),
                  scale(microdf$Turbidity),
                  scale(microdf$Temp),
                  scale(microdf$Chloride),
                  scale(microdf$SRP),
                  scale(microdf$Nitrate),
                  scale(microdf$Ammonium),
                  scale(microdf$DIN))
cs2 <- vifstep(cs1, th = 10)
cs2



#Cellulose data - respiration with DIN seperated and no conductivity
resp20 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=spongedf, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod20 <- dredge(resp20)
mod1_20 <- model.avg(get.models(mod20, subset = TRUE))
summary(mod1_20)
confint(mod1_20)
r.squaredGLMM(resp20)




#Wood Veneer - respiration with DIN seperated and no conductivity
resp21 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wooddf, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod21 <- dredge(resp21)
mod1_21 <- model.avg(get.models(mod21, subset = TRUE))
summary(mod1_21)
confint(mod1_21)
r.squaredGLMM(resp21)
















