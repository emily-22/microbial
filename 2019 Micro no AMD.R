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

# The VIF score indicates taht the same two variabels are colinear
# Going to use same selected models as we did with the 
# the full model 

#Cellulose - resp with all variabls
# we are not going to use this model, rather we'll use the smaller 
# selection of models 
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


#########################################################################################
#########################################################################################
#####
#####  Using the limited selection of models to predict resp on sponge
##### since we only have ~11 sites, we should be conservative 
##### about how many perdictors in our model
#######################################################################################
#######################################################################################


bResp_NOAMD.models<-list()
bResp_NOAMD.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[12]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[13]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[14]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[15]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[16]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bResp_NOAMD.models[[17]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)

## Creating a vector of names to trace back models in set
ModnamesbResp_NOAMD <- paste("model", 1:length(bResp_NOAMD.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = bResp_NOAMD.models, modnames = ModnamesbResp_NOAMD, sort = TRUE)
#############################
# Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt     LL
# model 7  4 146.25       0.00   0.34   0.34 -68.61
# model 13 5 147.16       0.91   0.21   0.55 -67.79
# model 14 5 147.18       0.93   0.21   0.76 -67.80
# model 11 5 148.39       2.14   0.11   0.87 -68.41
# model 12 5 148.79       2.54   0.09   0.97 -68.60

# The best model is #7 - Chloride, and model 13 and 14 are Chloride plus a form a nitrogen (Nitrate #13) or Ammonium (#14)


r.squaredGLMM(bResp_NOAMD.models[[7]])
r.squaredGLMM(bResp_NOAMD.models[[13]])
r.squaredGLMM(bResp_NOAMD.models[[14]])

# > r.squaredGLMM(bResp_NOAMD.models[[7]])
#             R2m       R2c
# [1,] 0.5661119 0.5661119
# > r.squaredGLMM(bResp_NOAMD.models[[13]])
#             R2m       R2c
# [1,] 0.5821424 0.5821424
# > r.squaredGLMM(bResp_NOAMD.models[[14]])
#            R2m       R2c
# [1,] 0.5819774 0.5819774

# So while we have an issue with the boundary singularity in the top modes
# we have better predictivie ability having removed the AMD sites. 

confint(bResp_NOAMD.models[[7]])
# Computing profile confidence intervals ...
#             2.5 %    97.5 %
#   .sig01          0.0000000 0.7946575
# .sigma          0.9265538 1.4404568
# (Intercept)     1.5346377 2.3361251
# scale(Chloride) 0.9218652 1.6979145

confint(bResp_NOAMD.models[[13]])
# > confint(bResp_NOAMD.models[[13]])
# Computing profile confidence intervals ...
#                      2.5 %    97.5 %
#  .sig01           0.0000000 0.6690573
# .sigma           0.9279010 1.4138576
# (Intercept)      1.5841863 2.2907434
# scale(Chloride)  0.9360866 1.6398929
# scale(Nitrate)  -0.5828259 0.1435785

confint(bResp_NOAMD.models[[14]])
# > confint(bResp_NOAMD.models[[14]])
# Computing profile confidence intervals ...
#                   2.5 %    97.5 %
#   .sig01           0.0000000 0.6811258
# .sigma           0.9231717 1.4141341
# (Intercept)      1.5713178 2.2994147
# scale(Chloride)  1.0281896 1.8410021
# scale(Ammonium) -0.6509414 0.1485567


# Summary of Spong Respiration #########
#########################################
# This these results indicate that something very similar to our model with the
# full set of sites, except that we have a more predictive power (highR2 value)
# This is good, but these results have the CrCs in it. I'm hoping the boundary singularity issues
# might be resolved if we drop the CrCs
#




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



###################################################################################################
#######################################   Respiration Models with Wood ############################
###################################################################################################




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



###############################   Wood Respiration with Restricted Set of Models
###############################
###############################


cResp_NOAMD.models<-list()
cResp_NOAMD.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[12]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[13]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[14]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[15]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[16]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[17]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamd, REML = FALSE)
cResp_NOAMD.models[[18]]  <- lmer( RespRateInd~  (1 | Stream), data=woodnoamd, REML = FALSE)

## Creating a vector of names to trace back models in set
ModnamescResp_NOAMD <- paste("model", 1:length(cResp_NOAMD.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = cResp_NOAMD.models, modnames = ModnamescResp_NOAMD, sort = TRUE)
#############################


# Model selection based on AICc:
  
#  K    AICc Delta_AICc AICcWt Cum.Wt    LL
# model 16 5 -117.44       0.00   0.37   0.37 64.83
# model 4  4 -116.14       1.31   0.19   0.56 62.78
# model 9  5 -115.90       1.55   0.17   0.73 64.06
# model 2  4 -114.57       2.88   0.09   0.82 62.00
# model 15 5 -114.09       3.35   0.07   0.88 63.16


# The top model (#16) is Nitrate + Ammonium indicating that best predictor for resp. on Wood
# using this restricted set of models is nitrogen. 

r.squaredGLMM(cResp_NOAMD.models[[16]])
# > r.squaredGLMM(cResp_NOAMD.models[[16]])
#      R2m       R2c
# [1,] 0.3368773 0.3368773
# Boundary singularity issue. 

confint(cResp_NOAMD.models[[16]])
# > confint(cResp_NOAMD.models[[16]])
# Computing profile confidence intervals ...
#                         2.5 %     97.5 %
#   .sig01          0.0000000000 0.02248582
# .sigma          0.0271406286 0.04410113
# (Intercept)     0.1069357980 0.13096359
# scale(Nitrate)  0.0070088567 0.03683589
# scale(Ammonium) 0.0004316996 0.02480815









###################################################################################################
#######################################   Breakdown Models with Wood ############################
###################################################################################################






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


##############
############## Sponge Breakdown using a more restricted set of models
##############
##############



bBreak_NOAMD.models<-list()
bBreak_NOAMD.models[[1]]  <- lmer( Breakdown~ scale(Conductivity) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[2]]  <- lmer( Breakdown~ scale(DIN) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[3]]  <- lmer( Breakdown~ scale(Temp) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[4]]  <- lmer( Breakdown~ scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[5]]  <- lmer( Breakdown~ scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[6]]  <- lmer( Breakdown~ scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[7]]  <- lmer( Breakdown~ scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[8]]  <- lmer( Breakdown~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[9]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[10]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[11]]  <- lmer( Breakdown~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[12]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[13]]  <- lmer( Breakdown~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[14]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[15]]  <- lmer( Breakdown~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[16]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[17]]  <- lmer( Breakdown~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamd, REML = FALSE)
bBreak_NOAMD.models[[18]]  <- lmer( Breakdown~  (1 | Stream), data=spongenoamd, REML = FALSE)

## Creating a vector of names to trace back models in set
ModnamesbBreak_NOAMD <- paste("model", 1:length(bBreak_NOAMD.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = bBreak_NOAMD.models, modnames = ModnamesbBreak_NOAMD, sort = TRUE)
#############################


# Model selection based on AICc:
  
#   K    AICc Delta_AICc AICcWt Cum.Wt     LL
# model 9  5 -588.27       0.00      1      1 299.92
# model 15 5 -357.18     231.09      0      1 184.38
# model 12 5 -350.24     238.03      0      1 180.91


# So this is pretty cool, the only model that had any preditict power was temp + nitrate.


r.squaredGLMM(bBreak_NOAMD.models[[9]])
# > r.squaredGLMM(bBreak_NOAMD.models[[9]])
# R2m       R2c
# [1,] 0.5468454 0.9999972

# OK, this is ... too good to believe .. and we need to do some investigating. I think is is about how you  have the 
# data entered. Each site will only have one breakdown rate number (beacsue breakdown intergrated overtime
# and comes from estimating a slope to the substrate mass remaining). So we need to get the data such that each site
# has one breakdown number that we are using for the regression stats. 



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
















####
#### 
#### WITHOUT CRAP STREAM ############
####
####
####


micronoamdCS <- read.csv("Emily Data No AMD or CS.csv")

spongenoamdCS <- subset(micronoamdCS, Substrate=="Cellulose Sponge" & Week.Removed==4)
woodnoamdCS <- subset(micronoamdCS, Substrate=="Wood Veneer" & Week.Removed==6)

#Sponge resp
Resp52.models<-list()
Resp52.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Resp52.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp52 <- paste("model", 1:length(Resp52.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp52.models, modnames = ModnamesResp52, sort = TRUE)

r.squaredGLMM(Resp52.models[[7]])
r.squaredGLMM(Resp52.models[[18]])
r.squaredGLMM(Resp52.models[[19]])
r.squaredGLMM(Resp52.models[[12]])
r.squaredGLMM(Resp52.models[[32]])

#Top 5 models for sponge resp are 7, 18, 19, 12, 32
#model 7: Chloride mR2 = 0.56, cR2=0.56
#model 18: Chloride + SRP mR2 = 0.57, cR2=0.58
#model 19: Chloride + Nitrate mR2 = 0.57, cR2=0.57
#model 12: Chloride + Temp mR2 = 0.57, cR2=0.57
#model 32: Chloride + Nitrate + SRP mR2 = 0.6, cR2=0.6


#Wood resp
Resp53.models<-list()
Resp53.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + DIN + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[12]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[18]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[19]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[20]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[21]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[22]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[23]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[24]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[25]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[26]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[27]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[28]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[29]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[30]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[31]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[32]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[33]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[34]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[35]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[36]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[37]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[38]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[39]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[40]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[41]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[42]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Resp53.models[[43]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesResp53 <- paste("model", 1:length(Resp53.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp53.models, modnames = ModnamesResp53, sort = TRUE)

r.squaredGLMM(Resp53.models[[4]])
r.squaredGLMM(Resp53.models[[19]])
r.squaredGLMM(Resp53.models[[41]])
r.squaredGLMM(Resp53.models[[2]])
r.squaredGLMM(Resp53.models[[21]])

#Top 5 models for wood resp are 4, 19, 41, 2, 21
#model 4: Nitrate mR2 = 0.29, cR2=0.29
#model 19: Chloride + Nitrate mR2 = 0.32, cR2=0.32
#model 41: Temp + Nitrate + SRP mR2 = 0.48, cR2=0.59
#model 2: DIN + Temp mR2 = 0.24, cR2=0.24
#model 21: Nitrate + SRP mR2 = 0.31, cR2=0.31








####
####
###
###      Sponge Breakdown
###



Break52.models<-list()
Break52.models[[1]]  <- lmer( Breakdown~ scale(Conductivity) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[2]]  <- lmer( Breakdown~ scale(DIN) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[3]]  <- lmer( Breakdown~ scale(Temp) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[4]]  <- lmer( Breakdown~ scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[5]]  <- lmer( Breakdown~ scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[6]]  <- lmer( Breakdown~ scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[7]]  <- lmer( Breakdown~ scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[8]]  <- lmer( Breakdown~ scale(Temp) + DIN + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[9]]  <- lmer( Breakdown~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[10]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[11]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[12]]  <- lmer( Breakdown~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[13]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[14]]  <- lmer( Breakdown~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[15]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[16]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[17]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[18]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[19]]  <- lmer( Breakdown~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[20]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[21]]  <- lmer( Breakdown~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[22]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[23]]  <- lmer( Breakdown~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[24]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[25]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[26]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[27]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[28]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[29]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[30]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[31]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[32]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[33]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[34]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[35]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[36]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[37]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[38]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[39]]  <- lmer( Breakdown~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[40]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[41]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[42]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
Break52.models[[43]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=spongenoamdCS, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesBreak52 <- paste("model", 1:length(Break52.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Break52.models, modnames = ModnamesBreak52, sort = TRUE)

r.squaredGLMM(Break52.models[[41]])
r.squaredGLMM(Break52.models[[34]])
r.squaredGLMM(Break52.models[[37]])
r.squaredGLMM(Break52.models[[40]])
r.squaredGLMM(Break52.models[[31]])

#Top 5 models for sponge breakdown are 41, 34, 37, 40, 31
#model 41: Temp + Nitrate + SRP mR2 = 0.86, cR2=1
#model 34: Conductivity + Temp + Nitrate mR2 = 0.86, cR2=1
#model 37: Conductivity + Nitrate + Ammonium mR2 = 0.80, cR2=1
#model 40: Temp + Nitrate + Ammonium mR2 = 0.86, cR2=1
#model 31: Chloride + Nitrate + Ammonium mR2 = 0.8, cR2=0.1


#Wood Breakdown
Break53.models<-list()
Break53.models[[1]]  <- lmer( Breakdown~ scale(Conductivity) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[2]]  <- lmer( Breakdown~ scale(DIN) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[3]]  <- lmer( Breakdown~ scale(Temp) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[4]]  <- lmer( Breakdown~ scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[5]]  <- lmer( Breakdown~ scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[6]]  <- lmer( Breakdown~ scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[7]]  <- lmer( Breakdown~ scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[8]]  <- lmer( Breakdown~ scale(Temp) + DIN + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[9]]  <- lmer( Breakdown~ scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[10]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[11]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[12]]  <- lmer( Breakdown~ scale(Temp) + scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[13]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[14]]  <- lmer( Breakdown~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[15]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[16]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[17]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[18]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[19]]  <- lmer( Breakdown~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[20]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[21]]  <- lmer( Breakdown~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[22]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[23]]  <- lmer( Breakdown~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[24]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Conductivity) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[25]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[26]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[27]]  <- lmer( Breakdown~ scale(Chloride) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[28]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[29]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[30]]  <- lmer( Breakdown~ scale(Chloride) + scale(Conductivity) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[31]]  <- lmer( Breakdown~ scale(Chloride) + scale(Ammonium) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[32]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[33]]  <- lmer( Breakdown~ scale(Chloride) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[34]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Nitrate) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[35]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[36]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Temp) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[37]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[38]]  <- lmer( Breakdown~ scale(Conductivity) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[39]]  <- lmer( Breakdown~ scale(Conductivity) + scale(SRP) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[40]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[41]]  <- lmer( Breakdown~ scale(Temp) + scale(Nitrate) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[42]]  <- lmer( Breakdown~ scale(Temp) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
Break53.models[[43]]  <- lmer( Breakdown~ scale(Nitrate) + scale(Ammonium) + scale(SRP) + (1 | Stream), data=woodnoamdCS, REML = FALSE)
#####################

## Creating a vector of names to trace back models in set
ModnamesBreak53 <- paste("model", 1:length(Break53.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Break53.models, modnames = ModnamesBreak53, sort = TRUE)

r.squaredGLMM(Break53.models[[41]])
r.squaredGLMM(Break53.models[[37]])
r.squaredGLMM(Break53.models[[25]])
r.squaredGLMM(Break53.models[[36]])
r.squaredGLMM(Break53.models[[28]])

#Top 5 models for wood resp are 41, 37, 25, 36, 28
#model 41: Temp + Nitrate + SRP mR2 = 0.84, cR2=1
#model 37: Conductivity + Nitrate + Ammonium mR2 = 0.82, cR2=1
#model 25: CHloride + Temp + Nitrate + SRP mR2 = 0.84, cR2=1
#model 36: Conductivity + Temp + Nitrate mR2 = 0.80, cR2=1
#model 28: Chloride + Conductivity + Nitrate mR2 = 0.81, cR2=1










