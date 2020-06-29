micronoamd <- read.csv("Emily Data No AMD.csv")

spongenoamd <- subset(micronoamd, Substrate=="Cellulose Sponge" & Week.Removed==4)
woodnoamd <- subset(micronoamd, Substrate=="Wood Veneer" & Week.Removed==6)


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






