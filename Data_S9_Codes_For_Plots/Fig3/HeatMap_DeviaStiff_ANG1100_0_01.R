library(ggplot2)
library(readr)
library(gridExtra)
library(splines)

source("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Plots/ColorTheme.R")
setwd("C:/Users/win/Desktop/CurrentProject/MSS_MultistableSpatialStructure/Results")

Info_ANG1100 <-
  read_csv("GenerateResult_ANG_ID1100.csv",
           col_names = c("Theta1_C1","Theta2_C1","Theta3_C1","Theta4_C1",
                         "Theta1_C2","Theta2_C2","Theta3_C2","Theta4_C2",
                         "Deviation_N","Stiffness_C1","Stiffness_C2","Stiffness_Mean",
                         "Flag_Converaged"),
           show_col_types = FALSE)

Plot_DeviaAvg_C1_11 <- ggplot(Info_ANG1100[Info_ANG1100$Theta1_C1==Info_ANG1100$Theta1_C1[1] &
                                           Info_ANG1100$Theta2_C1==Info_ANG1100$Theta2_C1[1] &
                                           Info_ANG1100$Flag_Converaged==1,],
                              aes(x=Theta1_C2,y=Theta2_C2))+
                            geom_raster(aes(fill=Deviation_N),interpolate = FALSE)+
                            scale_fill_gradientn(colors = colors_2,limits = c(0, 0.1),oob = scales::squish)+
                            scale_x_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            scale_y_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            theme_bw()+ 
                            theme(panel.grid.major.y=element_blank(),
                                panel.grid.minor.y=element_blank(),
                                panel.grid.major.x=element_blank(),
                                panel.grid.minor.x=element_blank())+
                            coord_fixed(ratio=1)+
                            labs(x = element_blank(), y = element_blank()) +
                            theme(
                                legend.title=element_blank(),
                                legend.position = "bottom",
                                legend.key.width = unit(3, "cm"),
                                # legend.key.width = unit(3, "cm"),
                                # axis.ticks = element_blank(), 
                                # axis.text = element_blank(),
                                axis.line = element_blank()
                            )

Plot_DeviaAvg_C1_12 <- ggplot(Info_ANG1100[Info_ANG1100$Theta1_C1==Info_ANG1100$Theta1_C1[1] &
                                           Info_ANG1100$Theta2_C1==Info_ANG1100$Theta2_C1[25] &
                                           Info_ANG1100$Flag_Converaged==1,],
                              aes(x=Theta1_C2,y=Theta2_C2))+
                            geom_raster(aes(fill=Deviation_N),interpolate = FALSE)+
                            scale_fill_gradientn(colors = colors_2,limits = c(0, 0.1),oob = scales::squish)+
                            scale_x_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            scale_y_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            theme_bw()+ 
                            theme(panel.grid.major.y=element_blank(),
                                panel.grid.minor.y=element_blank(),
                                panel.grid.major.x=element_blank(),
                                panel.grid.minor.x=element_blank())+
                            coord_fixed(ratio=1)+
                            labs(x = element_blank(), y = element_blank()) +
                            theme(
                                legend.title=element_blank(),
                                legend.position = "bottom",
                                legend.key.width = unit(3, "cm"),
                                # legend.key.width = unit(3, "cm"),
                                # axis.ticks = element_blank(), 
                                # axis.text = element_blank(),
                                axis.line = element_blank()
                            )

Plot_DeviaAvg_C1_13 <- ggplot(Info_ANG1100[Info_ANG1100$Theta1_C1==Info_ANG1100$Theta1_C1[1] &
                                           Info_ANG1100$Theta2_C1==Info_ANG1100$Theta2_C1[48] &
                                           Info_ANG1100$Flag_Converaged==1,],
                              aes(x=Theta1_C2,y=Theta2_C2))+
                            geom_raster(aes(fill=Deviation_N),interpolate = FALSE)+
                            scale_fill_gradientn(colors = colors_2,limits = c(0, 0.1),oob = scales::squish)+
                            scale_x_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            scale_y_continuous(limit=c(-pi,pi),breaks = c(-3, 0, 3))+
                            theme_bw()+ 
                            theme(panel.grid.major.y=element_blank(),
                                panel.grid.minor.y=element_blank(),
                                panel.grid.major.x=element_blank(),
                                panel.grid.minor.x=element_blank())+
                            coord_fixed(ratio=1)+
                            labs(x = element_blank(), y = element_blank()) +
                            theme(
                                legend.title=element_blank(),
                                legend.position = "bottom",
                                legend.key.width = unit(3, "cm"),
                                # legend.key.width = unit(3, "cm"),
                                # axis.ticks = element_blank(), 
                                # axis.text = element_blank(),
                                axis.line = element_blank()
                            )