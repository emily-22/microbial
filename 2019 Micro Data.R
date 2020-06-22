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












