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

resp1a <- lmer(RespRateInd ~ scale(Chloride) +
                 scale(H) +
                 scale(Nitrate) +
                 scale(Turbidity) +
                 (1 | Site), data=spongedf, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
r.squaredGLMM(resp1a)



plot(spongedf$Conductivity,spongedf$RespRateAvg)

plot(spongedf$Conductivity, spongedf$Nitrate)




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



#Model Lists

Resp1.models<-list()
Resp1.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
Resp1.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp1 <- paste("model", 1:length(Resp1.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp1.models, modnames = ModnamesResp1, sort = TRUE)

r.squaredGLMM(Resp1.models[[24]])
r.squaredGLMM(Resp1.models[[28]])
r.squaredGLMM(Resp1.models[[17]])
r.squaredGLMM(Resp1.models[[30]])
r.squaredGLMM(Resp1.models[[29]])

#Top 5 models for sponge resp are 24, 28, 17, 30, 29
#model 24: Chloride + Temp + Conductivity mR2 = 0.56, cR2=0.56
#model 28: Chloride + Conductivity + Nitrate mR2 = 0.55, cR2=0.55
#model 17: Chloride + Conductivity mR2 = 0.48, cR2=0.57
#model 30: Chloride + Conductivity + SRP mR2 = 0.49, cR2=0.58
#model 29: Chloride + Conductivity + Ammonium mR2 = 0.48, cR2=0.57


#Wood resp
Resp2.models<-list()
Resp2.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
Resp2.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp2 <- paste("model", 1:length(Resp2.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp2.models, modnames = ModnamesResp2, sort = TRUE)

r.squaredGLMM(Resp2.models[[37]])
r.squaredGLMM(Resp2.models[[31]])
r.squaredGLMM(Resp2.models[[20]])
r.squaredGLMM(Resp2.models[[16]])
r.squaredGLMM(Resp2.models[[34]])

#Top 5 models for wood resp are 37, 31, 20, 16, 34
#model 37: Conductivity + Nitrate + Ammonium mR2 = 0.36, cR2=0.36
#model 31: Chloride + Ammonium + Nitrate mR2 = 0.35, cR2=0.35
#model 20: Chloride + Ammonium mR2 = 0.27, cR2=0.27
#model 16: DIN + Temp mR2 = 0.27, cR2=0.27
#model 34: Conductivity + Ammonium mR2 = 0.30, cR2=0.30


