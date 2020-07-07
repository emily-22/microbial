microdf <- read.csv("Emily Data.csv")

spongedf <- subset(microdf, Substrate=="Cellulose Sponge" & Week.Removed==4)
wooddf <- subset(microdf, Substrate=="Wood Veneer" & Week.Removed==6)

library(lme4)
library(ggplot2)
library(ggmap)
library(AICcmodavg)
library(scatterplot3d)
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
                (1 | Stream), data=spongedf, REML = FALSE, na.action = "na.fail",
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
                (1 | Stream), data=wooddf, REML = FALSE, na.action = "na.fail",
              control = lmerControl(optCtrl = list(maxfun=20000)))
mod21 <- dredge(resp21)
mod1_21 <- model.avg(get.models(mod21, subset = TRUE))
summary(mod1_21)
confint(mod1_21)
r.squaredGLMM(resp21)



CellCPResp <- read.csv("2019CellRespCP.csv")

#Cellulose - cat plot with no din or conductivity
CellCPResp$Predictor_Variables<-factor(CellCPResp$Predictor_Variables,
                                       levels=CellCPResp$Predictor_Variables[order((CellCPResp$Estimate_Effects),
                                                                                   decreasing=TRUE)])

CellCPResp.plot<-ggplot(data=CellCPResp,aes(x=Estimate_Effects,y=Predictor_Variables))+
  geom_point(shape=1, size =5)+
  geom_errorbarh(aes(xmin=LCI,xmax=UCI),height=0.0,size=1, colour="black")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous("Mean Breakdown Parameter Estimate") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=16, colour="black"),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=0),
        plot.title = element_text(hjust = 0.5, vjust = -2))+
  theme(plot.title= element_text(size=16, face="bold"))
CellCPResp.plot


WoodCPResp <- read.csv("2019WoodRespCP.csv")

#Wood - cat plot with no din or conductivity
WoodCPResp$Predictor_Variables<-factor(WoodCPResp$Predictor_Variables,
                                       levels=WoodCPResp$Predictor_Variables[order((WoodCPResp$Estimate_Effects),
                                                                                   decreasing=TRUE)])

WoodCPResp.plot<-ggplot(data=WoodCPResp,aes(x=Estimate_Effects,y=Predictor_Variables))+
  geom_point(shape=1, size =5)+
  geom_errorbarh(aes(xmin=LCI,xmax=UCI),height=0.0,size=1, colour="black")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_x_continuous("Mean Breakdown Parameter Estimate") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=16, colour="black"),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=0),
        plot.title = element_text(hjust = 0.5, vjust = -2))+
  theme(plot.title= element_text(size=16, face="bold"))
WoodCPResp.plot


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









