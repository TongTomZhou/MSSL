library(ggplot2)
library(readr)
library(gridExtra)
library(splines)

source("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/ColorTheme.R")
setwd("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Results")

# source("C:/Users/zh0ut/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/ColorTheme.R")
# setwd("C:/Users/zh0ut/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Results")

Info_AllTolerance <-
  read_csv("GenerateResult_ANG_ID1000_AllTole.csv",
           col_names = c("Tolerance","Theta1_C1","Theta2_C1","Theta3_C1","Theta4_C1",
                         "Theta1_C2","Theta2_C2","Theta3_C2","Theta4_C2",
                         "Deviation_N","Stiffness_C1","Stiffness_C2","Stiffness_Mean",
                         "Flag_Converaged"),
           show_col_types = FALSE)

Info_0_01Tole <- Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,]


Plot_Devia_0_01 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,],
                          aes(x=Theta1_C1,y=Theta1_C2))+
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi))+
    geom_raster(aes(fill=Deviation_N),interpolate = FALSE)+
    scale_fill_gradientn(colors = colors_2,limits = c(0, 0.1))+
    theme_bw()+ theme(panel.grid.major.y=element_blank(),
                    panel.grid.minor.y=element_blank(),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank())+
    coord_fixed(ratio=1)+
    labs(x = element_blank(), y = element_blank()) +
    theme(legend.title=element_blank(),legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            # axis.ticks = element_blank(), 
            # axis.text = element_blank(),
            axis.line = element_blank())

Plot_Stiff1_0_01 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,],
                          aes(x=Theta1_C1,y=Theta1_C2))+
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi))+
    geom_raster(aes(fill=log(Stiffness_C1,base=10)),interpolate = FALSE)+
    scale_fill_gradientn(colors = colors_2,limits = c(-4, 1),oob = scales::squish)+
    theme_bw()+ theme(panel.grid.major.y=element_blank(),
                    panel.grid.minor.y=element_blank(),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank())+
    coord_fixed(ratio=1)+
    labs(x = element_blank(), y = element_blank()) +
    theme(legend.title=element_blank(),legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            # axis.ticks = element_blank(), 
            # axis.text = element_blank(),
            axis.line = element_blank())
Plot_Stiff2_0_01 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,],
                          aes(x=Theta1_C1,y=Theta1_C2))+
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi))+
    geom_raster(aes(fill=log(Stiffness_C2,base=10)),interpolate = FALSE)+
    scale_fill_gradientn(colors = colors_2,limits = c(-4, 1),oob = scales::squish)+
    theme_bw()+ theme(panel.grid.major.y=element_blank(),
                    panel.grid.minor.y=element_blank(),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank())+
    coord_fixed(ratio=1)+
    labs(x = element_blank(), y = element_blank()) +
    theme(legend.title=element_blank(),legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            # axis.ticks = element_blank(), 
            # axis.text = element_blank(),
            axis.line = element_blank())
Plot_Stiff_0_01 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,],
                          aes(x=Theta1_C1,y=Theta1_C2))+
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi))+
    geom_raster(aes(fill=log(Stiffness_Mean,base=10)),interpolate = FALSE)+
    scale_fill_gradientn(colors = colors_2,limits = c(-4, 1),oob = scales::squish)+
    theme_bw()+ theme(panel.grid.major.y=element_blank(),
                    panel.grid.minor.y=element_blank(),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank())+
    coord_fixed(ratio=1)+
    labs(x = element_blank(), y = element_blank()) +
    theme(legend.title=element_blank(),legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            # axis.ticks = element_blank(), 
            # axis.text = element_blank(),
            axis.line = element_blank())


