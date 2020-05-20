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


head(microdf)

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
#looked at with vifstep for collinearity found that conductivity and DIN have
# a collinearity problem



#Cellulose data - respiration with DIN seperated and no conductivity
resp1 <- lmer(RespRateInd ~ scale(DO) + 
               scale(Temp) + 
               scale(Turbidity) +
               scale(Chloride) +
               scale(H) +
               scale(SRP) +
               scale(Nitrate) +
               scale(Ammonium) +
               (1 | Site), data=spongedf, REML = FALSE, na.action = "na.fail",
             control = lmerControl(optCtrl = list(maxfun=20000)))
mod1 <- dredge(resp1)
mod1_1 <- model.avg(get.models(mod1, subset = TRUE))
summary(mod1_1)
confint(mod1_1)
r.squaredGLMM(resp1)

#Cellulose data - respiration with no chloride, pH, or DIN
resp3 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=spongedf, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod3 <- dredge(resp3)
mod1_3 <- model.avg(get.models(mod3, subset = TRUE))
summary(mod1_3)
confint(mod1_3)
r.squaredGLMM(resp3)


#Wood data - respiration with DIN seperated and no conductivity
resp2 <- lmer(RespRateInd ~ scale(DO) + 
               scale(Temp) + 
               scale(Turbidity) +
               scale(Chloride) +
               scale(H) +
               scale(SRP) +
               scale(Nitrate) +
               scale(Ammonium) +
               (1 | Site), data=wooddf, REML = FALSE, na.action = "na.fail",
             control = lmerControl(optCtrl = list(maxfun=20000)))
mod2 <- dredge(resp2)
mod2_2 <- model.avg(get.models(mod2, subset = TRUE))
summary(mod2_2)
confint(mod2_2)
r.squaredGLMM(resp2)

#Wood data - respiration with no chloride, pH, or DIN
resp4 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wooddf, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod4 <- dredge(resp4)
mod2_4 <- model.avg(get.models(mod4, subset = TRUE))
summary(mod2_4)
confint(mod2_4)
r.squaredGLMM(resp4)






###data without NrCr, ChCr01, or CsCr

mdf2 <- read.csv("Emily Data no NR CS CH.csv")

cmdf2 <- subset(mdf2, Substrate=="Cellulose Sponge" & Week.Removed==4)
wmdf2 <- subset(mdf2, Substrate=="Wood Veneer" & Week.Removed==6)


cs3 <- data.frame(scale(mdf2$Conductivity),
                  scale(mdf2$H),
                  scale(mdf2$DO),
                  scale(mdf2$Turbidity),
                  scale(mdf2$Temp),
                  scale(mdf2$Chloride),
                  scale(mdf2$SRP),
                  scale(mdf2$Nitrate),
                  scale(mdf2$Ammonium),
                  scale(mdf2$DIN))
cs4 <- vifstep(cs3, th = 10)
cs4
#found colinearity with DIN, Conductivity, Ammonium


#Cellulose data - respiration with DIN seperated and no conductivity
resp8 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=cmdf2, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod8 <- dredge(resp8)
mod1_8 <- model.avg(get.models(mod8, subset = TRUE))
summary(mod1_8)
confint(mod1_8)
r.squaredGLMM(resp8)

#Cellulose data - respiration with no chloride, pH, or DIN
resp5 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=cmdf2, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod5 <- dredge(resp5)
mod1_5 <- model.avg(get.models(mod5, subset = TRUE))
summary(mod1_5)
confint(mod1_5)
r.squaredGLMM(resp5)

#Cellulose data - respiration with  no conductivity, Ammonium, DIN
resp14 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                (1 | Site), data=cmdf2, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod14 <- dredge(resp14)
mod1_14 <- model.avg(get.models(mod14, subset = TRUE))
summary(mod1_14)
confint(mod1_14)
r.squaredGLMM(resp14)


#Wood data - respiration with DIN seperated and no conductivity
resp6 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wmdf2, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod6 <- dredge(resp6)
mod2_6 <- model.avg(get.models(mod6, subset = TRUE))
summary(mod2_6)
confint(mod2_6)
r.squaredGLMM(resp6)

#Wood data - respiration with no chloride, pH, or DIN
resp9 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wmdf2, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod9 <- dredge(resp9)
mod2_9 <- model.avg(get.models(mod9, subset = TRUE))
summary(mod2_9)
confint(mod2_9)
r.squaredGLMM(resp9)

#Wood data - respiration with  no conductivity, Ammonium, DIN
resp15 <- lmer(RespRateInd ~ scale(DO) + 
                 scale(Temp) + 
                 scale(Turbidity) +
                 scale(Chloride) +
                 scale(H) +
                 scale(SRP) +
                 scale(Nitrate) +
                 (1 | Site), data=wmdf2, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod15 <- dredge(resp15)
mod2_15 <- model.avg(get.models(mod15, subset = TRUE))
summary(mod2_15)
confint(mod2_15)
r.squaredGLMM(resp15)

###data without NrCr, ChCr01

mdf3 <- read.csv("Emily Data no NR CH.csv")

cmdf3 <- subset(mdf3, Substrate=="Cellulose Sponge" & Week.Removed==4)
wmdf3 <- subset(mdf3, Substrate=="Wood Veneer" & Week.Removed==6)


cs5 <- data.frame(scale(mdf2$Conductivity),
                  scale(mdf2$H),
                  scale(mdf2$DO),
                  scale(mdf2$Turbidity),
                  scale(mdf2$Temp),
                  scale(mdf2$Chloride),
                  scale(mdf2$SRP),
                  scale(mdf2$Nitrate),
                  scale(mdf2$Ammonium),
                  scale(mdf2$DIN))
cs6 <- vifstep(cs5, th = 10)
cs6
#colinearity with DIN, Conductivity, Ammonium


#Cellulose data - respiration with DIN seperated and no conductivity
resp10 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=cmdf3, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod10 <- dredge(resp10)
mod1_10 <- model.avg(get.models(mod10, subset = TRUE))
summary(mod1_10)
confint(mod1_10)
r.squaredGLMM(resp10)

#Cellulose data - respiration with no chloride, pH, or DIN
resp11 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=cmdf3, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod11 <- dredge(resp11)
mod1_11 <- model.avg(get.models(mod11, subset = TRUE))
summary(mod1_11)
confint(mod1_11)
r.squaredGLMM(resp11)


#Wood data - respiration with DIN seperated and no conductivity
resp12 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Chloride) +
                scale(H) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wmdf3, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod12 <- dredge(resp12)
mod2_12 <- model.avg(get.models(mod12, subset = TRUE))
summary(mod2_12)
confint(mod2_12)
r.squaredGLMM(resp12)

#Wood data - respiration with no chloride, pH, or DIN
resp13 <- lmer(RespRateInd ~ scale(DO) + 
                scale(Temp) + 
                scale(Turbidity) +
                scale(Conductivity) +
                scale(SRP) +
                scale(Nitrate) +
                scale(Ammonium) +
                (1 | Site), data=wmdf3, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod13 <- dredge(resp13)
mod2_13 <- model.avg(get.models(mod13, subset = TRUE))
summary(mod2_13)
confint(mod2_13)
r.squaredGLMM(resp13)






