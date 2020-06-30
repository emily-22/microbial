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
############################ makeing lists of models#
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


### got this error message when adding these Error in formatCands(cand.set) : Functions do not support mixture of model classes ####
#Resp50.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[35]]  <- lmer( RespRateInd~ scale(Conductivty) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#Resp50.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp50 <- paste("model", 1:length(Resp50.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp50.models, modnames = ModnamesResp50, sort = TRUE)
#############################


#Top 5 models were 7, 19, 20, 12, 33
#model 7: CHloride
#model 19: CHloride + Nitrate
#model 20: CHloride + Ammonium
#model 12: Temp + Chloride
#model 33: Chloride + Ammonium + SRP




