library(ggplot2)
library(readr)
library(gridExtra)
library(splines)

source("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/ColorTheme.R")
setwd("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/NEW PLOTS2/Fig2/Codes for plots")

Info_All <-
  read_csv("MATX_ConfigurationAnalysis_Bistable4R.csv",
           col_names = c("theta3","theta1","theta2","theta4","arc31","arc24","ang_31_24",
                         "areaH1","areaH2","areaH3","areaH4"),
           show_col_types = TRUE)

Plot_theta3_thetas<- ggplot(Info_All,aes(x=theta3))+
  geom_line(aes(y = theta1),size=0.8,color="#58595B")+
  #geom_point(aes(y = theta1), fill="#58595B",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = theta2),size=0.8,color="#FBB040")+
  #geom_point(aes(y = theta2), fill="#FBB040",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = theta4),size=0.8,color="#F15A29")+
  #geom_point(aes(y = theta4), fill="#F15A29",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  theme(
    plot.background = element_rect(fill = "transparent"), 
  )+
  scale_x_continuous(
    breaks = c(0, pi/6, pi/3, pi/2, 2*pi/3),
    labels = c("0", "π/6", "π/3", "π/2", "2π/3"), 
    limits = c(0, pi*2/3)
  )+
  scale_y_continuous(
    breaks = c(0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi),
    labels = c("0", "π/3", "2π/3", "π", "4π/3", "5π/3", "2π"), 
    limits = c(0, 2*pi)
  )

Plot_theta3_thetas


Plot_theta3_arc<- ggplot(Info_All,aes(x=theta3))+
  geom_line(aes(y = arc31),size=0.8,color="#FBB040")+
  #geom_point(aes(y = arc31), fill="#FBB040",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = arc24),size=0.8,color="#F15A29")+
  #geom_point(aes(y = arc31), fill="#F15A29",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  theme(
    plot.background = element_rect(fill = "transparent"), 
  )+
  scale_x_continuous(
    breaks = c(0, pi/6, pi/3, pi/2, 2*pi/3),
    labels = c("0", "π/6", "π/3", "π/2", "2π/3"), 
    limits = c(0, pi*2/3)
  )+
  scale_y_continuous(
    breaks = c(0.1,0.13,0.16,0.19,0.22,0.25, 0.28),
    labels = c("0.1", "0.13", "0.16", "0.19", "0.22", "0.25", "0.28"), 
    limits = c(0.1, 0.28)
  )

Plot_theta3_arc


Plot_theta3_arcang<- ggplot(Info_All,aes(x=theta3))+
  geom_line(aes(y = ang_31_24),size=0.8,color="#58595B")+
  #geom_point(aes(y = ang_31_24), fill="#FBB040",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  theme(
    plot.background = element_rect(fill = "transparent"), 
  )+
  scale_x_continuous(
    breaks = c(0, pi/6, pi/3, pi/2, 2*pi/3),
    labels = c("0", "π/6", "π/3", "π/2", "2π/3"), 
    limits = c(0, pi*2/3)
  )+
  scale_y_continuous(
    breaks = seq(pi/2, 3*pi/4, by = 1/24*pi),
    labels = c("pi/2", "13*pi/24", "7*pi/12", "5*pi/8", "2*pi/3", "17*pi/24", "3*pi/4"), 
    limits = c(pi/2, 3*pi/4)
  )

Plot_theta3_arcang

Plot_theta3_projarea<- ggplot(Info_All,aes(x=theta3))+
  geom_line(aes(y = areaH1),size=0.8,color="#FBB040")+
  #geom_point(aes(y = theta2), fill="#FBB040",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = areaH2),size=0.8,color="#F15A29")+
  #geom_point(aes(y = theta4), fill="#F15A29",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = areaH3),size=0.8,color="#84B7DC")+
  #geom_point(aes(y = theta2), fill="#FBB040",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  geom_line(aes(y = areaH4),size=0.8,color="#5C7AB8")+
  #geom_point(aes(y = theta4), fill="#F15A29",color="transparent",size = 0.3, shape = 21, stroke = 0.3)+
  theme(
    plot.background = element_rect(fill = "transparent"), 
  )+
  scale_x_continuous(
    breaks = c(0, pi/6, pi/3, pi/2, 2*pi/3),
    labels = c("0", "π/6", "π/3", "π/2", "2π/3"), 
    limits = c(0, pi*2/3)
  )+
  scale_y_continuous(
    breaks = seq(0, 0.036, by = 0.006),
    labels = c("0", "0.006", "0.012", "0.018", "0.024", "0.030", "0.036"), 
    limits = c(0, 0.036)
  )

Plot_theta3_projarea

