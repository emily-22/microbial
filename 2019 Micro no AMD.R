micronoamd <- read.csv("Emily Data No AMD.csv")

spongenoamd <- subset(micronoamd, Substrate=="Cellulose Sponge" & Week.Removed==4)
woodnoamd <- subset(micronoamd, Substrate=="Wood Veneer" & Week.Removed==6)

library(lme4)
library(ggplot2)
library(ggmap)
library(AICcmodavg)
library(scatterplot3d)
library(devtools)
library(cowplot)
library(MuMIn)
library(usdm)

plot(micronoamd$DIN, micronoamd$Conductivity)



cs10 <- data.frame(scale(micronoamd$Conductivity),
                  scale(micronoamd$H),
                  scale(micronoamd$DO),
                  scale(micronoamd$Turbidity),
                  scale(micronoamd$Temp),
                  scale(micronoamd$Chloride),
                  scale(micronoamd$SRP),
                  scale(micronoamd$Nitrate),
                  scale(micronoamd$Ammonium),
                  scale(micronoamd$DIN))
cs11 <- vifstep(cs10, th = 10)
cs11

#Cellulose - resp with all
resp50 <- lmer(RespRateInd ~ scale(DO) + 
                 scale(Temp) + 
                 scale(Turbidity) +
                 scale(Chloride) +
                 scale(H) +
                 scale(SRP) +
                 scale(DIN) +
                 scale(Conductivity) +
                 scale(Nitrate) +
                 scale(Ammonium) +
                 (1 | Stream), data=spongenoamd, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod50 <- dredge(resp50)
mod1_50 <- model.avg(get.models(mod50, subset = TRUE))
summary(mod1_50)
confint(mod1_50)
r.squaredGLMM(resp50)



#Cellulose data - respiration with DIN seperated and no conductivity
resp53 <- lmer(RespRateInd ~ scale(DO) + 
                 scale(Temp) + 
                 scale(Turbidity) +
                 scale(Chloride) +
                 scale(H) +
                 scale(SRP) +
                 scale(Nitrate) +
                 scale(Ammonium) +
                 (1 | Stream), data=spongenoamd, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod53 <- dredge(resp53)
mod1_53 <- model.avg(get.models(mod53, subset = TRUE))
summary(mod1_53)
confint(mod1_53)


#Cellulose - resp with no DIN, Conductivity, Temp, DO, or pH
resp51 <- lmer(RespRateInd ~ scale(Turbidity) +
                 scale(Chloride) +
                 scale(SRP) +
                 scale(Nitrate) +
                 scale(Ammonium) +
                 (1 | Stream), data=spongenoamd, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod51 <- dredge(resp51)
mod1_51 <- model.avg(get.models(mod51, subset = TRUE))
summary(mod1_51)
confint(mod1_51)
r.squaredGLMM(resp51)

#Cellulose - resp with no do, nutrient data, temp, or do
resp52 <- lmer(RespRateInd ~ scale(Turbidity) +
                 scale(Chloride) +
                 scale(Conductivity) +
                 scale(H) +
                 (1 | Stream), data=spongenoamd, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod52 <- dredge(resp52)
mod1_52 <- model.avg(get.models(mod52, subset = TRUE))
summary(mod1_52)
confint(mod1_52)
r.squaredGLMM(resp52)





#Wood data - respiration with DIN seperated and no conductivity
resp54 <- lmer(RespRateInd ~ scale(DO) + 
                 scale(Temp) + 
                 scale(Turbidity) +
                 scale(Chloride) +
                 scale(H) +
                 scale(SRP) +
                 scale(Nitrate) +
                 scale(Ammonium) +
                 (1 | Stream), data=woodnoamd, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod54 <- dredge(resp54)
mod1_54 <- model.avg(get.models(mod54, subset = TRUE))
summary(mod1_54)
confint(mod1_54)





############################ 
############################ makeing lists of models for sponge#
Resp50.models<-list()
Resp50.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Resp50.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp50 <- paste("model", 1:length(Resp50.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp50.models, modnames = ModnamesResp50, sort = TRUE)
#############################
r.squaredGLMM(Resp50.models[[7]])
r.squaredGLMM(Resp50.models[[19]])
r.squaredGLMM(Resp50.models[[20]])
r.squaredGLMM(Resp50.models[[12]])
r.squaredGLMM(Resp50.models[[33]])

#Top 5 sponge models were 7, 19, 20, 12, 33

#model 7: Chloride mR2 = 0.57, cR2=0.57 BSE
#model 19: Chloride + Nitrate mR2 = 0.58, cR2=0.58 BSE
#model 20: Chloride + Ammonium mR2 = 0.58, cR2=0.58 BSE
#model 12: Temp + Chloride mR2 = 0.57, cR2=0.57 BSE
#model 33: Chloride + Ammonium + SRP mR2 = 0.597, cR2=0.60 BSE

plot(spongenoamd$Chloride, spongenoamd$RespRateInd)
plot(spongenoamd$Chloride, spongenoamd$Nitrate)
plot(spongenoamd$Chloride, spongenoamd$Ammonium)
plot(spongenoamd$Chloride, spongenoamd$Temp)




#####################################                          ###
##################################### creating models for Wood ###
Resp51.models<-list()
Resp51.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Resp51.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp51 <- paste("model", 1:length(Resp51.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp51.models, modnames = ModnamesResp51, sort = TRUE)
#############################
r.squaredGLMM(Resp51.models[[22]])
r.squaredGLMM(Resp51.models[[28]])
r.squaredGLMM(Resp51.models[[4]])
r.squaredGLMM(Resp51.models[[10]])
r.squaredGLMM(Resp51.models[[31]])
r.squaredGLMM(Resp51.models[[43]])

#Top 5 wood models were 22, 28, 4, 10, 31
#model 22: Nitrate + Ammonium, mR2 = 0.34, cR2=0.34 BSE
#model 28: Chloride + Nitrate + Conductivity mR2 = 0.37, cR2=0.37 BSE
#model 4: Nitrate mR2 = 0.25, cR2=0.25 BSE
#model 10: Temp + Nitrate mR2 = 0.30, cR2=0.30 BSE
#model 31: Chloride + Ammonium + Nitrate mR2 = 0.36, cR2=0.36 BSE

plot(woodnoamd$Nitrate, woodnoamd$Ammonium)
plot(woodnoamd$Nitrate, woodnoamd$RespRateInd)
abline(lm(woodnoamd$RespRateInd~woodnoamd$Nitrate))
plot(woodnoamd$Ammonium, woodnoamd$RespRateInd)
plot(woodnoamd$Ammonium+woodnoamd$Nitrate, woodnoamd$RespRateInd) 
plot(woodnoamd$DIN, woodnoamd$RespRateInd) 
plot(woodnoamd$Nitrate, woodnoamd$Temp)


####BREAKDOWN

##sponge models
Break50.models<-list()
Break50.models[[1]]  <- lmer( Breakdown~ scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[2]]  <- lmer( Breakdown~ scale(DIN) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[3]]  <- lmer( Breakdown~ scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[4]]  <- lmer( Breakdown~ scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[5]]  <- lmer( Breakdown~ scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[6]]  <- lmer( Breakdown~ scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[7]]  <- lmer( Breakdown~ scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[8]]  <- lmer( Breakdown~ scale(Temp) + scale(DIN) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[9]]  <- lmer( Breakdown~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[10]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[11]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[12]]  <- lmer( Breakdown~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[13]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[14]]  <- lmer( Breakdown~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[15]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[16]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[17]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[18]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[19]]  <- lmer( Breakdown~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[20]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[21]]  <- lmer( Breakdown~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[22]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[23]]  <- lmer( Breakdown~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[24]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[25]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[26]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[27]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[28]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[29]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[30]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[31]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[32]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[33]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[34]]  <- lmer( Breakdown~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[35]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[36]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
Break50.models[[37]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)


## Creating a vector of names to trace back models in set
ModnamesBreak50 <- paste("model", 1:length(Break50.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Break50.models, modnames = ModnamesBreak50, sort = TRUE)
#############################
r.squaredGLMM(Break50.models[[29]])
r.squaredGLMM(Break50.models[[35]])
r.squaredGLMM(Break50.models[[31]])
r.squaredGLMM(Break50.models[[8]])
r.squaredGLMM(Break50.models[[10]])

#Top 5 models for sponge breakdown were 29, 35, 31, 8, 10
#model 29: Conductivity + Temp + Nitrate mR2 = 0.87, cR2=1
#model 35: Temp + Nitrate + SRP mR2 = 0.87, cR2=1
#model 31: Conductivity + Temp + SRP mR2 = 0.86, cR2=1
#model 8: Temp + DIN mR2 = 0.42, cR2=1
#model 10: Temp + Nitrate mR2 = 0.55, cR2=1


plot(spongenoamd$Temp, spongenoamd$DIN)
plot(spongenoamd$Temp, spongenoamd$Nitrate)

##Wood Breakdown
Break51.models<-list()
Break51.models[[1]]  <- lmer( Breakdown~ scale(Conductivity) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[2]]  <- lmer( Breakdown~ scale(DIN) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[3]]  <- lmer( Breakdown~ scale(Temp) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[4]]  <- lmer( Breakdown~ scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[5]]  <- lmer( Breakdown~ scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[6]]  <- lmer( Breakdown~ scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[7]]  <- lmer( Breakdown~ scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[8]]  <- lmer( Breakdown~ scale(Temp) + scale(DIN) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[9]]  <- lmer( Breakdown~ scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[10]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[11]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[12]]  <- lmer( Breakdown~ scale(Temp) + scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[13]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[14]]  <- lmer( Breakdown~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[15]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[16]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[17]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[18]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[19]]  <- lmer( Breakdown~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[20]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[21]]  <- lmer( Breakdown~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[22]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[23]]  <- lmer( Breakdown~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[24]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[25]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[26]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[27]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[28]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[29]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[30]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[31]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[32]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[33]]  <- lmer( Breakdown~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[34]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[35]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
Break51.models[[36]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)


## Creating a vector of names to trace back models in set
ModnamesBreak51 <- paste("model", 1:length(Break51.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Break51.models, modnames = ModnamesBreak51, sort = TRUE)
#############################
r.squaredGLMM(Break51.models[[29]])
r.squaredGLMM(Break51.models[[26]])
r.squaredGLMM(Break51.models[[31]])
r.squaredGLMM(Break51.models[[27]])
r.squaredGLMM(Break51.models[[34]])


#Top 5 models for wood breakdown were 29, 26, 31, 27, 34
#model 29: Conductivity + Temp + Nitrate mR2 = 0.77, cR2=1
#model 26: Chloride + Temp + SRP mR2 = 0.81, cR2=1
#model 31: Conductivity + Temp + SRP mR2 = 0.82, cR2=1
#model 27: Chloride + Ammonium + Nitrate mR2 = 0.78, cR2=1
#model 34: Temp + Nitrate + SRP mR2 = 0.86, cR2=1







