waterchem2 <- read.csv("Water Chem.csv")

library(ggplot2)

levels(waterchem2$Date)
waterchem2$Date <- factor(waterchem2$Date, levels = c("2/7/18", "2/10/18", "2/20/18", "3/9/18","3/24/18","3/25/18",
                                                      "4/2/18", "4/9/18", "4/12/18","4/13/18","4/16/18", "4/20/18", "5/21/2018",
                                                      "6/6/2018", "6/28/2018", "7/3/2018", "7/9/2018", "7/11/2018", "7/17/2018", 
                                                      "7/24/2018",  "7/26/2018", "8/7/2018", "9/3/2018", "10/3/2018", "10/4/2018",
                                                      "10/10/2018", "10/11/2018", "10/17/2018", "10/18/2018", "10/24/2018",
                                                      "10/30/2018", "10/31/2018", "11/01/2018", "11/05/2018", "11/07/2018", 
                                                      "11/08/2018", "11/13/2018", "11/14/2018", "11/20/2018", "11/29/2018", 
                                                      "12/5/2018", "12/10/2018", "12/13/2018", "12/17/2018", "12/27/2018", 
                                                      "1/3/2019", "1/16/2019", "1/29/2019", "2/15/2019",  "2/19/2019", "2/27/2019", 
                                                      "3/3/2019", "3/6/2019", "3/12/2019", "3/26/2019", "3/29/2019", "4/9/2019", 
                                                      "4/12/2019", "4/16/2019", "5/2/2019", "5/16/2019", "5/22/2019", "5/29/2019",  
                                                      "5/31/2019", "6/3/2019", "6/4/2019", "6/11/2019",  "6/12/2019",  "6/19/2019",
                                                      "6/21/2019",  "6/26/2019", "6/27/2019",  "6/28/2019", "7/1/2019",  "7/2/2019",
                                                      "7/3/2019",  "7/4/2019", "7/8/2019", "7/9/2019", "7/10/2019", "7/11/2019", 
                                                      "7/12/2019",  "7/16/2019", "7/17/2019",  "7/18/2019",  "7/19/2019",   "7/23/2019", 
                                                      "7/24/2019",  "7/25/2019",  "7/26/2019",  "7/29/2019",  
                                                      "7/30/2019",  "7/31/2019", "8/1/2019",   "8/2/2019", "8/5/2019",   "8/6/2019",
                                                      "8/7/2019",   "8/8/2019", "8/20/2019",  "8/21/2019", "8/22/2019",  "8/28/2019",  
                                                      "8/30/2019", "9/8/2019", "9/18/2019", "9/30/2019", "10/2/2019", "10/3/2019",
                                                      "10/21/2019", "10/23/2019", "10/28/2019", "11/1/2019", "11/4/2019", "11/11/2019", 
                                                      "11/13/2019","11/14/2019", "11/16/2019", "11/19/2019","11/20/2019", "11/21/2019", 
                                                      "11/25/2019", "12/5/2019","12/11/2019", "12/18/2019", "12/19/2019", "12/26/2019",   
                                                      "1/8/2020", "1/10/2020", "1/20/2020","1/28/2020", "2/5/2020", 
                                                      "2/10/2020", "2/18/2020", "2/20/2020",  "2/25/2020", "2/28/2020", "3/4/2020",
                                                      "3/9/2020", "3/26/2020", "4/10/2020", "5/1/2020", "5/8/2020", "5/13/2020",
                                                      "5/20/2020", "5/28/2020", "6/5/2020", "6/9/2020", "6/10/2020", "6/17/2020",
                                                      "6/22/2020", "6/23/2020", "6/29/2020", "7/1/2020", "7/8/2020", "7/9/2020",
                                                      "7/29/2020", "8/3/2020", "8/4/2020", "8/11/2020", "8/12/2020", "8/19/2020",
                                                      "8/26/2020", "9/2/2020", "9/9/2020", "9/16/2020", "9/23/2020", "9/30/2020",
                                                      "10/2/2020", "10/7/2020", "10/14/2020", "10/20/2020", "11/18/2020", "11/24/2020",
                                                      "12/3/2020"))

ggplot(waterchem2,aes(x = Date ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.x = element_text(size=5)) +
  labs(y = "Conductivity ms/cm", x = "Date")

ggplot(waterchem2,aes(x = Water.Temp ,y = Conductivity, colour = Site, group = Site)) + 
  geom_point() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Water Temperature C")

ggplot(waterchem2,aes(x = Chloride ,y = Conductivity, colour = Site, group = Site)) + 
  geom_point() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Chloride")

################################################################################
#using julian day


ggplot(waterchem2,aes(x = Julian.Day ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Julian Day")

##Wheeling Creek and tribs
WeCr <- read.csv("WeCr.csv")

ggplot(WeCr,aes(x = Julian.Day ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Julian Day")


##North fork and short creek
NfSh <- read.csv("NfSh.csv")

ggplot(NfSh,aes(x = Julian.Day ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Julian Day")

##Buffalo creek
BfCr <- read.csv("BfCr.csv")

ggplot(BfCr,aes(x = Julian.Day ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Julian Day")

##AMD sites
AMDs <- read.csv("AMDs.csv")

ggplot(AMDs,aes(x = Julian.Day ,y = Conductivity, colour = Site, group = Site)) + 
  geom_line() +
  theme_classic() +
  labs(y = "Conductivity ms/cm", x = "Julian Day")




























