emmadf <- read.csv("Emma Data.csv")

emmasponge <- subset(emmadf, Substrate=="Cellulose" & WeekRemoved==4)
emmawood <- subset(emmadf, Substrate=="Wood" & WeekRemoved==4)

cs5 <- data.frame(scale(emmadf$SPC),
                  scale(emmadf$H),
                  scale(emmadf$DO),
                  scale(emmadf$Temp),
                  scale(emmadf$SRP),
                  scale(emmadf$Nitrate),
                  scale(emmadf$Ammonium),
                  scale(emmadf$DIN))
cs6 <- vifstep(cs5, th = 10)
cs6
# Nitrate, DIN, and Conductivity all have collinearity problems


#Cellulose - respiration with Nitrate, DIN, and Conductivity removed
resp22 <- lmer(RespRateIndividual ~ scale(DO) + 
                 scale(Temp) + 
                 scale(H) +
                 scale(SRP) +
                 scale(Ammonium) +
                 (1 | SiteName), data=emmasponge, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod22 <- dredge(resp22)
mod1_22 <- model.avg(get.models(mod22, subset = TRUE))
summary(mod1_22)
confint(mod1_22)
r.squaredGLMM(resp22)

#Wood - respiration with Nitrate, DIN, and Conductivity removed
resp23 <- lmer(RespRateIndividual ~ scale(DO) + 
                 scale(Temp) + 
                 scale(H) +
                 scale(SRP) +
                 scale(Ammonium) +
                 (1 | SiteName), data=emmawood, REML = FALSE, na.action = "na.fail",
               control = lmerControl(optCtrl = list(maxfun=20000)))
mod23 <- dredge(resp23)
mod1_23 <- model.avg(get.models(mod23, subset = TRUE))
summary(mod1_23)
confint(mod1_23)
r.squaredGLMM(resp23)


CellCPResp2 <- read.csv("2018CellRespCP.csv")

#Cellulose - cat plot with no din, nitrate, or conductivity
CellCPResp2$Predictor_Variables<-factor(CellCPResp2$Predictor_Variables,
                                       levels=CellCPResp2$Predictor_Variables[order((CellCPResp2$Estimate_Effects),
                                                                                   decreasing=TRUE)])

CellCPResp2.plot<-ggplot(data=CellCPResp2,aes(x=Estimate_Effects,y=Predictor_Variables))+
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
CellCPResp2.plot




