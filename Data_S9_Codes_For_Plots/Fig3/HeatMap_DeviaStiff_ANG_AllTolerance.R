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

Plot_Devia_0_09 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,],
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
Plot_Stiff_0_09 <- ggplot(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,],
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


Data_DeviaStiffAvg_AllTole <- data.frame(Tolerance=seq(0.01,0.09,0.01),
    Deviation_N=c(colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.02,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.03,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.04,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.05,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.06,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.07,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.08,10]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,10])),
    StiffC1Avg=c(colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.02,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.03,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.04,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.05,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.06,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.07,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.08,11]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,11])), 
    StiffC2Avg=c(colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.02,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.03,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.04,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.05,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.06,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.07,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.08,12]),
                    colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,12])), 
    StiffAvg=c(colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.01,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.02,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.03,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.04,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.05,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.06,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.07,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.08,13]),
                colMeans(Info_AllTolerance[Info_AllTolerance$Tolerance==0.09,13])))

Data_Stiff_AllTole <- data.frame(Config=c(rep("C1",9),rep("C2",9)),
                                 Tolerance=c(seq(0.01,0.09,0.01),seq(0.01,0.09,0.01)),
                                 Stiffness=c(Data_DeviaStiffAvg_AllTole$StiffC1Avg,
                                             Data_DeviaStiffAvg_AllTole$StiffC2Avg))

Plot_Stiff12_AllTole <- ggplot(Data_Stiff_AllTole,aes(x=Tolerance,y=Stiffness,fill=Config))+
                         geom_col(position = "dodge",color="black",width = 0.005)+
                         scale_fill_manual(values = c(colors_2[6],colors_2[9]))+
                         scale_x_continuous(breaks = seq(0.01, 0.09, 0.01))+
                         scale_y_continuous(limit=c(0.00,1.00),breaks = seq(0, 1.0, 0.2))+
                         theme_bw()+ theme(panel.grid.major.y=element_blank(),
                                           panel.grid.minor.y=element_blank(),
                                           panel.grid.major.x=element_blank(),
                                           panel.grid.minor.x=element_blank())+
                                     coord_fixed(ratio=0.05)+
                                     labs(x = element_blank(), y = element_blank()) +
                                     theme(
                                        legend.title=element_blank(),
                                        legend.position = "bottom",
                                        # legend.key.width = unit(3, "cm"),
                                        # axis.ticks = element_blank(), 
                                        # axis.text = element_blank(),
                                        axis.line = element_blank()
                                        )
Plot_StiffAVG_AllTole <- ggplot(Data_DeviaStiffAvg_AllTole,aes(x=Tolerance,y=StiffAvg))+
                         geom_col(position = "dodge",color="black",width = 0.005,fill=colors_2[9])+
                         scale_x_continuous(breaks = seq(0.01, 0.09, 0.01))+
                         scale_y_continuous(limit=c(0.00,1.00),breaks = seq(0, 1.0, 0.2))+
                         theme_bw()+ theme(panel.grid.major.y=element_blank(),
                                           panel.grid.minor.y=element_blank(),
                                           panel.grid.major.x=element_blank(),
                                           panel.grid.minor.x=element_blank())+
                                    #  coord_fixed(ratio=0.05)+
                                     labs(x = element_blank(), y = element_blank()) +
                                     theme(
                                        legend.title=element_blank(),
                                        legend.position = "bottom",
                                        # legend.key.width = unit(3, "cm"),
                                        # axis.ticks = element_blank(), 
                                        # axis.text = element_blank(),
                                        axis.line = element_blank()
                                        )
Plot_DeviaAvg_AllTole <- ggplot(Data_DeviaStiffAvg_AllTole,aes(x=Tolerance))+
                         geom_line(aes(y=Deviation_N),color="black")+
                         geom_point(aes(y=Deviation_N), fill=colors_2[2],shape = 21, color = "black", size = 6)+
                         scale_x_continuous(breaks = seq(0.01, 0.09, 0.01))+
                         scale_y_continuous(limit=c(0.00,0.05),breaks = seq(0, 0.05, 0.01))+
                         theme_bw()+ theme(panel.grid.major.y=element_blank(),
                                           panel.grid.minor.y=element_blank(),
                                           panel.grid.major.x=element_blank(),
                                           panel.grid.minor.x=element_blank())+
                                    #  coord_fixed(ratio=0.05)+
                                     labs(x = element_blank(), y = element_blank()) +
                                     theme(
                                        legend.title=element_blank(),
                                        legend.position = "bottom",
                                        # legend.key.width = unit(3, "cm"),
                                        # axis.ticks = element_blank(), 
                                        # axis.text = element_blank(),
                                        axis.line = element_blank()
                                        )

