library(ggplot2)
library(readr)
library(gridExtra)
library(splines)

source("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/ColorTheme.R")
setwd("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/NEW PLOTS2/Fig2/Raw Data")

Info_AllForce <-
  read_csv("ForceCurve_4R.csv",
           col_names = c("Type","DIST_Actu","FORCE"),
           show_col_types = TRUE)

Plot_Force<- ggplot(Info_AllForce,aes(x=DIST_Actu,y=FORCE,group=Type))+
  coord_cartesian(xlim = c(0, 0.18), ylim = c(0, 0.3))+
  geom_line()+
  geom_point(aes(color = Type, fill = Type), size = 1.5, shape = 21, stroke = 0.3)+
  scale_shape_manual(values = c("S1toS2" = 1, "S2toS1" = 4)) +
  scale_color_manual(values = c("S1toS2" = "black", "S2toS1" = "black"))+
  scale_fill_manual(values = c("S1toS2" = "#FDF1AC", "S2toS1" = "#F15A29")) +
  facet_grid(Type ~ .,scales = "fixed", space = "fixed")+
  theme(
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "lightgray"), 
    panel.spacing = unit(0.5, "lines")
    )  
  

Plot_Force
options(units = "mm")
