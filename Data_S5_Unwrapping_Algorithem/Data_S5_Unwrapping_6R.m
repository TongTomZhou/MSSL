clear;clc;close all;
% load('RUST_SepcialCase_6R_Keshe.mat');
tool=TOOLS_SharedFunction();

MATX_Points_C1=MATX_PointsSBest_NoContact(1:12,:);
MATX_Points_C2=MATX_PointsSBest_NoContact(13:24,:);
MATX_Points_C3=MATX_PointsSBest_NoContact(25:36,:);

figure(1)
subplot(2,3,1); tool.plot6RStruc(MATX_Points_C1,'','c');
subplot(2,3,2); tool.plot6RStruc(MATX_Points_C2,'','m');
subplot(2,3,3); tool.plot6RStruc(MATX_Points_C3,'','y');


TSOR_Points(:,:,1)=MATX_Points_C1;  TSOR_Points(:,:,2)=MATX_Points_C2;  TSOR_Points(:,:,3)=MATX_Points_C3;
TSOR_Points_Slip=tool.SlipHinges_6R(TSOR_Points,[0,20,-20,0,20+20,10]); %[0,0,15,0,38,0]

MATX_Points_C1=TSOR_Points_Slip(:,:,1);  MATX_Points_C2=TSOR_Points_Slip(:,:,2);
MATX_Points_C3=TSOR_Points_Slip(:,:,3);


figure(1)
subplot(2,3,4); tool.plot6RStruc(MATX_Points_C1,'','r');
subplot(2,3,5); tool.plot6RStruc(MATX_Points_C2,'','g');
subplot(2,3,6); tool.plot6RStruc(MATX_Points_C3,'','b');

tool.calStiff6R(MATX_Points_C1,[1,2,7,8],1)
tool.calStiff6R(MATX_Points_C2,[1,2,7,8],1)
tool.calStiff6R(MATX_Points_C3,[1,2,7,8],1)


%% Defination
MATX_Points=MATX_Points_C1;
LISTr_Topo_Tetra61=[6,12,7,1];
LISTr_Topo_Tetra12=[7,1,2,8];    LISTr_Topo_Tetra23=[2,8,9,3];      LISTr_Topo_Tetra34=[9,3,4,10];
LISTr_Topo_Tetra45=[4,10,11,5];  LISTr_Topo_Tetra56=[11,5,6,12];

LISTr_Topo_ExptTetra61=[2,3,4,5,8,9,10,11];
LISTr_Topo_ExptTetra6112=[3,4,5,9,10,11];
LISTr_Topo_ExptTetra611223=[4,5,10,11];

%% Rotation 1
MATX_Points_Tetra61=MATX_Points(LISTr_Topo_Tetra61,:);
MATX_Points_Tetra12=MATX_Points(LISTr_Topo_Tetra12,:);
VALE_Angle_Hinge1=tool.CalProjANG(MATX_Points_Tetra61(1,:)',...
                                  MATX_Points_Tetra12(4,:)',...
                                  MATX_Points_Tetra61(3,:)',MATX_Points_Tetra61(4,:)');
VALE_Angle_Hinge1=pi-VALE_Angle_Hinge1;
MATX_Points_R1=MATX_Points;
MATX_Points_R1(LISTr_Topo_ExptTetra61,:)=tool.RotaPoints(MATX_Points_R1(LISTr_Topo_ExptTetra61,:),...
                                                         MATX_Points_Tetra61(3,:)',MATX_Points_Tetra61(4,:)',...
                                                         VALE_Angle_Hinge1);
MATX_Points_Tetra56_R1=tool.RotaPoints(MATX_Points(LISTr_Topo_Tetra56,:),...
                                       MATX_Points_Tetra61(3,:)',MATX_Points_Tetra61(4,:)',...
                                       VALE_Angle_Hinge1);

%% Rotation 2
MATX_Points_Tetra12=MATX_Points_R1(LISTr_Topo_Tetra12,:);
MATX_Points_Tetra23=MATX_Points_R1(LISTr_Topo_Tetra23,:);
VALE_Angle_Hinge2=tool.CalProjANG(MATX_Points_Tetra12(2,:)',...
                                  MATX_Points_Tetra23(3,:)',...
                                  MATX_Points_Tetra12(3,:)',MATX_Points_Tetra12(4,:)');
VALE_Angle_Hinge2=pi-VALE_Angle_Hinge2;
MATX_Points_R2=MATX_Points_R1;
MATX_Points_R2(LISTr_Topo_ExptTetra6112,:)=tool.RotaPoints(MATX_Points_R2(LISTr_Topo_ExptTetra6112,:),...
                                                           MATX_Points_Tetra12(3,:)',MATX_Points_Tetra12(4,:)',...
                                                           VALE_Angle_Hinge2);
MATX_Points_Tetra56_R2=tool.RotaPoints(MATX_Points_Tetra56_R1,...
                                       MATX_Points_Tetra12(3,:)',MATX_Points_Tetra12(4,:)',...
                                       VALE_Angle_Hinge2);

%% Rotation 3
MATX_Points_Tetra23=MATX_Points_R2(LISTr_Topo_Tetra23,:);
MATX_Points_Tetra34=MATX_Points_R2(LISTr_Topo_Tetra34,:);
VALE_Angle_Hinge3=tool.CalProjANG(MATX_Points_Tetra23(2,:)',...
                                  MATX_Points_Tetra34(3,:)',...
                                  MATX_Points_Tetra23(3,:)',MATX_Points_Tetra23(4,:)');
VALE_Angle_Hinge3=pi-VALE_Angle_Hinge3;
MATX_Points_R3=MATX_Points_R2;
MATX_Points_R3(LISTr_Topo_ExptTetra611223,:)=tool.RotaPoints(MATX_Points_R3(LISTr_Topo_ExptTetra611223,:),...
                                                             MATX_Points_Tetra23(3,:)',MATX_Points_Tetra23(4,:)',...
                                                             VALE_Angle_Hinge3);
MATX_Points_Tetra56_R3=tool.RotaPoints(MATX_Points_Tetra56_R2,...
                                       MATX_Points_Tetra23(3,:)',MATX_Points_Tetra23(4,:)',...
                                       VALE_Angle_Hinge3);

%% Rotation 4
MATX_Points_Tetra34=MATX_Points_R3(LISTr_Topo_Tetra34,:);
MATX_Points_Tetra45=MATX_Points_R3(LISTr_Topo_Tetra45,:);
VALE_Angle_Hinge4=tool.CalProjANG(MATX_Points_Tetra34(2,:)',...
                                  MATX_Points_Tetra45(3,:)',...
                                  MATX_Points_Tetra34(3,:)',MATX_Points_Tetra34(4,:)');
VALE_Angle_Hinge4=pi-VALE_Angle_Hinge4;
MATX_Points_R4=MATX_Points_R3;
MATX_Points_R4(LISTr_Topo_ExptTetra611223,:)=tool.RotaPoints(MATX_Points_R4(LISTr_Topo_ExptTetra611223,:),...
                                                             MATX_Points_Tetra34(3,:)',MATX_Points_Tetra34(4,:)',...
                                                             VALE_Angle_Hinge4);
MATX_Points_Tetra56_R4=tool.RotaPoints(MATX_Points_Tetra56_R3,...
                                       MATX_Points_Tetra34(3,:)',MATX_Points_Tetra34(4,:)',...
                                       VALE_Angle_Hinge4);

%% Rotation 5
MATX_Points_Tetra45=MATX_Points_R4(LISTr_Topo_Tetra45,:);
%MATX_Points_Tetra56=MATX_Points_R4(LISTr_Topo_Tetra56,:);
VALE_Angle_Hinge5=tool.CalProjANG(MATX_Points_Tetra45(1,:)',...
                                  MATX_Points_Tetra56_R4(4,:)',...
                                  MATX_Points_Tetra45(3,:)',MATX_Points_Tetra45(4,:)');
VALE_Angle_Hinge5=pi-VALE_Angle_Hinge5;
MATX_Points_R5=MATX_Points_R4;
MATX_Points_Tetra56_R5=tool.RotaPoints(MATX_Points_Tetra56_R4,...
                                       MATX_Points_Tetra45(3,:)',MATX_Points_Tetra45(4,:)',...
                                       VALE_Angle_Hinge5);

MATX_Points_Tetra61=MATX_Points_R5(LISTr_Topo_Tetra61,:);
MATX_Points_Tetra23=MATX_Points_R5(LISTr_Topo_Tetra23,:);   MATX_Points_Tetra34=MATX_Points_R5(LISTr_Topo_Tetra34,:);
MATX_Points_Tetra34=MATX_Points_R5(LISTr_Topo_Tetra34,:);   MATX_Points_Tetra45=MATX_Points_R5(LISTr_Topo_Tetra45,:);



figure(2);
tool.plot6RStruc(MATX_Points,'','r')
%tool.plot6RStruc(MATX_Points_R1,'','b')

tetramesh([1,2,3,4],MATX_Points_Tetra61,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra12,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra23,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra34,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra45,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra56_R4,...
          'FaceColor','g','FaceAlpha',0.2,'EdgeColor','k');
hold on
tetramesh([1,2,3,4],MATX_Points_Tetra56_R5,...
          'FaceColor','c','FaceAlpha',0.2,'EdgeColor','k');
hold on
axis equal




















