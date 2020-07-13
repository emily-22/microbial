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

### testing for colinearity in the perdictor variables
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
# VIF score indicates DIN and conductivity has a collinearity problem,
# and need to be removed from the data set.

# This will plot the x~y relationship of all the variables. 
pairs(~Conductivity+    H+     DO+ Turbidity+ Temp+    Chloride+      SRP+      Nitrate+
          +       Ammonium+      DIN, data=spongedf)
#From looking at these gaphs you can see why the conductivity and DIN were removed. DIN appears to be 
# strongly correlated with Nitrate and Ammonia - this is expected, and Conductivity and Chloride appear
# to be correlated. This would be expected also.




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


# top models
# Component models: 
# df logLik    AICc delta weight
# 145       6  77.01 -139.40  0.00   0.15 # Ammonia + pH + temp
# 45        5  74.98 -138.15  1.25   0.08 # ph + temp
# 1245      7  77.47 -137.32  2.08   0.05 # Ammonia + Chloride + Temp + Nitrate
# Term codes: 
#  scale(Ammonium)  scale(Chloride)        scale(DO)         scale(H)   scale(Nitrate) 
#              1                2                3                4                5 
#scale(SRP)      scale(Temp) scale(Turbidity) 
#       6                7                8 

# Model-averaged coefficients:  
#  (full average) 
# Estimate Std. Error Adjusted SE z value Pr(>|z|)    
# (Intercept)       0.1139393  0.0055560   0.0057714  19.742   <2e-16 ***
# scale(Ammonium)   0.0087341  0.0089963   0.0091193   0.958   0.3382    
# scale(H)         -0.0245663  0.0124279   0.0126010   1.950   0.0512 .  
# scale(Nitrate)    0.0170616  0.0095076   0.0096741   1.764   0.0778 .  
# scale(Chloride)  -0.0041334  0.0096223   0.0097416   0.424   0.6713    
# scale(SRP)       -0.0006401  0.0054611   0.0055655   0.115   0.9084    
# scale(Turbidity) -0.0005980  0.0041431   0.0042446   0.141   0.8880    
# scale(DO)        -0.0018846  0.0055695   0.0056635   0.333   0.7393    
# scale(Temp)      -0.0007919  0.0047093   0.0048192   0.164   0.8695    

#(conditional average) 
# Estimate Std. Error Adjusted SE z value Pr(>|z|)    
# (Intercept)       0.113939   0.005556    0.005771  19.742  < 2e-16 ***
#  scale(Ammonium)   0.014060   0.007443    0.007681   1.831  0.06717 .  
#  scale(H)         -0.027820   0.009186    0.009449   2.944  0.00324 ** 
#  scale(Nitrate)    0.019622   0.007329    0.007575   2.590  0.00959 ** 
#  scale(Chloride)  -0.013269   0.013267    0.013544   0.980  0.32721    
#  scale(SRP)       -0.002680   0.010927    0.011146   0.240  0.80998    
#  scale(Turbidity) -0.002950   0.008817    0.009052   0.326  0.74452    
#  scale(DO)        -0.006927   0.008893    0.009109   0.760  0.44696    
#  scale(Temp)      -0.003595   0.009519    0.009765   0.368  0.71274    
# ---
  

# Relative variable importance: 
#   scale(H) scale(Nitrate) scale(Ammonium) scale(Chloride) scale(DO)
# Importance:          0.88     0.87           0.62            0.31            0.27     
# N containing models:  128      128            128             128             128     
# scale(SRP) scale(Temp) scale(Turbidity)
# Importance:          0.24       0.22        0.20            
# N containing models:  128        128         128            

# COnfidence interval for whole model parameters 
# > confint(mod1_21)
#   2.5 %       97.5 %
#  (Intercept)       0.1026276611  0.125250978
# scale(Ammonium)  -0.0009939754  0.029113651
# scale(H)         -0.0463403885 -0.009300549
# scale(Nitrate)    0.0047748006  0.034469866
# scale(Chloride)  -0.0398140527  0.013275740
# scale(SRP)       -0.0245249864  0.019165032
# scale(Turbidity) -0.0206914389  0.014791725
# scale(DO)        -0.0247799653  0.010925773
# scale(Temp)      -0.0227348881  0.015544201
# > r.squaredGLMM(resp21)
# R2m       R2c
# [1,] 0.4981961 0.4981961






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



# Here instead of running all possible combindations of models with function "dredge" 
# we build our candiates models we chose to make combinations of models with no
# more than three variables. This was to in part try to avoid the "boundary singlarity issues
# from over fitting the models
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
#Resp1.models[[13]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Temp) + (1 | Stream), data=spongedf, REML = FALSE)
#Resp1.models[[14]]  <- lmer( RespRateInd~ scale(Conductivity) +scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
#Resp1.models[[15]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
#Resp1.models[[16]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
#Resp1.models[[17]]  <- lmer( RespRateInd~ scale(Conductivity) + scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)

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
# ok, most of the models run just fine, with some exceptions
# model 24 - boundary singularity
# model 28 - boundary singularity
# model 38 - model did not converge - I think this indicates that we shouldn't trust it 
# if that model is shown to be a top model, perhaps some colinearity issues
# we still need to check correlation strength of each model.

## Creating a vector of names to trace back models in set
ModnamesResp1 <- paste("model", 1:length(Resp1.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = Resp1.models, modnames = ModnamesResp1, sort = TRUE)

#Top 5 models for sponge resp are 24, 28, 17, 30, 29
# Top two models have boundary singularity issues - this means our model is overfitted,
# basically we don't need the LMER (ie.random effect... (1| site) part of our model)
# we might just want to run these models over again with just a linear model (lm),
# however that could affect our AIC - I'll look into this. 

# model 24 6 166.84       0.00   0.56   0.56 -76.46
# model 28 6 168.05       1.21   0.31   0.87 -77.07
# model 17 5 171.01       4.18   0.07   0.93 -79.84
# model 30 6 172.49       5.65   0.03   0.97 -79.29
# model 29 6 173.59       6.75   0.02   0.99 -79.84

r.squaredGLMM(Resp1.models[[24]])
r.squaredGLMM(Resp1.models[[28]])
r.squaredGLMM(Resp1.models[[17]])
r.squaredGLMM(Resp1.models[[30]])
r.squaredGLMM(Resp1.models[[29]])

#model 24: Chloride + Temp + Conductivity mR2 = 0.56, cR2=0.56
#model 28: Chloride + Conductivity + Nitrate mR2 = 0.55, cR2=0.55
#model 17: Chloride + Conductivity mR2 = 0.48, cR2=0.57
#model 30: Chloride + Conductivity + SRP mR2 = 0.49, cR2=0.58
#model 29: Chloride + Conductivity + Ammonium mR2 = 0.48, cR2=0.57


summary(lm( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity), data=spongedf, REML = FALSE))
#Coefficients:
#                       Estimate   Std. Error  t value  Pr(>|t|)    
#(Intercept)            1.9070      0.1581  12.064      5.36e-16 ***
#  scale(Chloride)       3.8665     0.5324   7.262     3.28e-09 ***
#  scale(Temp)          -0.6829     0.2185  -3.125     0.00304 ** 
#  scale(Conductivity)  -3.8199     0.5788  -6.600     3.31e-08 ***
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.129 on 47 degrees of freedom
#Multiple R-squared:  0.5571,	Adjusted R-squared:  0.5288 
#F-statistic:  19.7 on 3 and 47 DF,  p-value: 2.06e-08

  
summary(lm( RespRateInd~ scale(Chloride) + scale(Temp) + scale(Conductivity), data=spongedf, REML = FALSE))
  
  


confint((Resp1.models[[24]]))
#                          2.5 %     97.5 %
#  .sig01               0.0000000  0.6559585
#.sigma               0.8881173  1.3335128
#(Intercept)          1.5756617  2.2391412
#scale(Chloride)      2.7878644  4.9515209
#scale(Temp)         -1.1452066 -0.2163937
#scale(Conductivity) -4.9808215 -2.6594951



##### Since Conductivity and DIN had high VIF we should not include these with other variables inthe
##### models, however we can include them on their own. a
aResp1.models<-list()
aResp1.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
aResp1.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)


### limiting the model set to no more than two parameters and only including DIN and Condcutivity on their
### own becasue of VIF scores - i.e. those two variables are correlated with the other variables too much
### we don't want to not look at them on their own incase they are important. 

bResp1.models<-list()
bResp1.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[12]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[13]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[14]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[15]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[16]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)
bResp1.models[[17]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=spongedf, REML = FALSE)

## Creating a vector of names to trace back models in set
bModnamesResp1 <- paste("model", 1:length(bResp1.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = bResp1.models, modnames = bModnamesResp1, sort = TRUE)

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt     LL
#model 13 5 178.95       0.00   0.18   0.18 -83.81
#model 7  4 180.09       1.14   0.10   0.29 -85.61
#model 4  4 180.46       1.51   0.09   0.37 -85.79
#model 12 5 180.47       1.52   0.09   0.46 -84.57
#model 6  4 180.75       1.80   0.07   0.53 -85.94
#model 5  4 181.27       2.32   0.06   0.59 -86.20
#model 3  4 181.30       2.35   0.06   0.65 -86.22


# Top models with restricted model design
# 13 = Chloride + Nitrate
# 7 = Chloride
# 4 = Nitrate
# 12 = Chloride + SRP
# 6 = SRP
# 5 = Ammonium
# 3 = Temp


r.squaredGLMM(bResp1.models[[13]])
r.squaredGLMM(bResp1.models[[7]])
r.squaredGLMM(bResp1.models[[4]])
r.squaredGLMM(bResp1.models[[12]])
r.squaredGLMM(Resp1.models[[6]])

# r.squaredGLMM(bResp1.models[[13]])
#      R2m       R2c
# [1,] 0.248669 0.5604809

# r.squaredGLMM(bResp1.models[[7]])
# R2m       R2c
# [1,] 0.1013571 0.5743215

# r.squaredGLMM(bResp1.models[[4]])
# R2m       R2c
# [1,] 0.06505229 0.5972224

# r.squaredGLMM(bResp1.models[[12]])
# R2m       R2c
# [1,] 0.1881552 0.5793134

# r.squaredGLMM(Resp1.models[[6]])
# R2m       R2c
# [1,] 0.04033858 0.5751821


confint(bResp1.models[[13]])
#> confint(bResp1.models[[13]])
#Computing profile confidence intervals ...
#                        2.5 %    97.5 %
#  .sig01           0.46021129 1.6938726
#.sigma           0.88221961 1.3684866
#(Intercept)      1.26556545 2.7288278
#scale(Chloride)  0.01562696 1.4817865
#scale(Nitrate)  -1.33807857 0.0245815

confint(bResp1.models[[7]])
#> confint(bResp1.models[[7]])
#Computing profile confidence intervals ...
#                      2.5 %   97.5 %
#  .sig01           0.6590273 2.033027
# .sigma           0.8830404 1.367452
# (Intercept)      1.1544630 2.900778
# scale(Chloride) -0.2937579 1.335221


confint(bResp1.models[[4]])
#> confint(bResp1.models[[4]])
#Computing profile confidence intervals ...
#                        2.5 %    97.5 %
#  .sig01          0.7475690 2.1608081
# .sigma          0.8741705 1.3508881
# (Intercept)     1.2451282 3.0633174
# scale(Nitrate) -1.1727697 0.2984845


summary(bResp1.models[[13]])
# > summary(bResp1.models[[13]])
# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: RespRateInd ~ scale(Chloride) + scale(Nitrate) + (1 | Stream)
# Data: spongedf


# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)       1.9839     0.3346   5.929
# scale(Chloride)   0.7885     0.3402   2.318
# scale(Nitrate)   -0.6806     0.3281  -2.075

# Correlation of Fixed Effects:
#   (Intr) scl(C)
# scal(Chlrd) -0.181       
# scale(Ntrt)  0.042 -0.357


##############   Take away from building the model
# When the entire data set is considered our top model (the best thing we have to predicet respiration on sponge
# (a subsitute for leaves and labile carbon) is Cholride concentration + the effects of Nitrate.
# Chloride has an clearly *significant* effect on respiration while the effect of nitrate appears to be largely
# negative, which is different than many other studies have found. Our models are not very strong though, our best
# model "Chloride + Nitrate" to explain the variation in respiration rate on the sponge only explained 25% of the variation observed (the 
# marginal R2 value). 

# so across all our sites and the wide variation they have in chemistry and size, its hard to explain what is
# really controling the speed of microbial activity on labile organic material in streams.










###
###
###
###
###






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





# Using a more limited set of models so that we are not overfitting
# This is the same set of models (i.e. water chemistry predictors) that were
# chosen to use with the sponge data. 


bResp2.models<-list()
bResp2.models[[1]]  <- lmer( RespRateInd~ scale(Conductivity) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[2]]  <- lmer( RespRateInd~ scale(DIN) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[3]]  <- lmer( RespRateInd~ scale(Temp) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[4]]  <- lmer( RespRateInd~ scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[5]]  <- lmer( RespRateInd~ scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[6]]  <- lmer( RespRateInd~ scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[7]]  <- lmer( RespRateInd~ scale(Chloride) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[8]]  <- lmer( RespRateInd~ scale(Temp) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[9]]  <- lmer( RespRateInd~ scale(Temp) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[10]]  <- lmer( RespRateInd~ scale(Temp) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[11]]  <- lmer( RespRateInd~ scale(Temp) + scale(Chloride) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[12]]  <- lmer( RespRateInd~ scale(Chloride) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[13]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Nitrate) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[14]]  <- lmer( RespRateInd~ scale(Chloride) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[15]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(SRP) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[16]]  <- lmer( RespRateInd~ scale(Nitrate) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)
bResp2.models[[17]]  <- lmer( RespRateInd~ scale(SRP) + scale(Ammonium) + (1 | Stream), data=wooddf, REML = FALSE)

## Creating a vector of names to trace back models in set
bModnamesResp2 <- paste("model", 1:length(bResp2.models), sep = " ")

##generate AICc table from candidate models so that you can control the model
aictab(cand.set = bResp2.models, modnames = bModnamesResp2, sort = TRUE)

#Model selection based on AICc:
#  
#  K    AICc Delta_AICc AICcWt Cum.Wt    LL
#model 14 5 -133.61       0.00   0.58   0.58 72.71
#model 13 5 -129.99       3.62   0.09   0.67 70.90
#model 1  4 -129.18       4.43   0.06   0.74 69.18

### So similar to the sponge data, or best model includes Chloride, and a nitrogen, this time Ammonium. The next 
# best model, #13, is not very good compariatively (AICc score great than 2), however that model is 
# Chloride + Nitrate. The boundary singularity error is present indicating that we are overfitting in 
# model 14, which means that we're really artifically increasing the "goodness" of our model.

r.squaredGLMM(bResp2.models[[14]])
#> r.squaredGLMM(bResp2.models[[14]])
#       R2m       R2c
#[1,] 0.2674541 0.2674541

# So the Boundary Singularity error is presented by have the same marginal and conditional R2 value. If we run a 
# basic lm, with the "Random Effect" we will should get a R2 number that is pretty close to 0.26. 

r.squaredGLMM(bResp2.models[[13]])
# > r.squaredGLMM(bResp2.models[[13]])
# R2m      R2c
# [1,] 0.2014021 0.218222

r.squaredGLMM(bResp2.models[[1]])
# > r.squaredGLMM(bResp2.models[[1]])
# R2m       R2c
# [1,] 0.1183853 0.1183853

# only a little of variation in respiration on wood is explained by these models, about 20% by Chloride + Nitrate
# and only 11 percent by Conductivity alone (also BSE error)



confint(bResp2.models[[14]])
# > confint(bResp2.models[[14]])
# Computing profile confidence intervals ...
# 2.5 %       97.5 %
#   .sig01           0.000000000  0.019383829
# .sigma           0.030500391  0.047653842
# (Intercept)      0.101877059  0.126051619
# scale(Chloride) -0.036517018 -0.009215494
# scale(Ammonium)  0.006232609  0.033530257

# Chloride has a clear negative influence while Ammonium has a clear postive influence on respiration rate.



confint(bResp2.models[[13]])
# Computing profile confidence intervals ...
# 2.5 %      97.5 %
#   .sig01           0.0000000000  0.03753218
# .sigma           0.0295620201  0.04991054
# (Intercept)      0.0998754000  0.13181681
# scale(Chloride) -0.0414015023 -0.00561729
# scale(Nitrate)   0.0004584025  0.04115979

# Chloride has a clear negative influence while nitrate has a clear postive influence on respiration rate. 

confint(bResp2.models[[13]])
# Computing profile confidence intervals ...
# 2.5 %      97.5 %
#   .sig01           0.0000000000  0.03753218
# .sigma           0.0295620201  0.04991054
# (Intercept)      0.0998754000  0.13181681
# scale(Chloride) -0.0414015023 -0.00561729
# scale(Nitrate)   0.0004584025  0.04115979


summary(bResp2.models[[14]])
# Fixed effects:
#  Estimate Std. Error t value
# (Intercept)      0.113940   0.006005  18.975
# scale(Chloride) -0.022865   0.006793  -3.366
# scale(Ammonium)  0.019881   0.006793   2.927

# Correlation of Fixed Effects:
#   (Intr) scl(C)
# scal(Chlrd)  0.000       
# scal(Ammnm)  0.000 -0.445

