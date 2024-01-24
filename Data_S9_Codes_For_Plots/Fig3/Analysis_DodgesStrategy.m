clear;clc;close all;
tool=TOOLS_SharedFunction;

load('RUST_SpecialCase_4R_Miao.mat');

SlipDist=zeros(1,4);
[pointsDOC1,pointsDPC1,pointsSOC1,pointsSPC1]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                          MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(4,:));
[pointsDOC2,pointsDPC2,pointsSOC2,pointsSPC2]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                          MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(5,:));
MATX_PointsDH_C1=[pointsDOC1;pointsDPC1];   MATX_PointsDH_C2=[pointsDOC2;pointsDPC2];
MATX_Points_C1=[pointsSOC1;pointsSPC1];     MATX_Points_C2=[pointsSOC2;pointsSPC2];
SlipDist_Set=[(MATX_Points_SpecialCase(1,:)-MATX_PointsDH_C1(1,:))/...
                ((MATX_PointsDH_C1(5,:)-MATX_PointsDH_C1(1,:))/norm(MATX_PointsDH_C1(5,:)-MATX_PointsDH_C1(1,:))),...
              (MATX_Points_SpecialCase(2,:)-MATX_PointsDH_C1(2,:))/...
                ((MATX_PointsDH_C1(6,:)-MATX_PointsDH_C1(2,:))/norm(MATX_PointsDH_C1(6,:)-MATX_PointsDH_C1(2,:))),...
              (MATX_Points_SpecialCase(3,:)-MATX_PointsDH_C1(3,:))/...
                ((MATX_PointsDH_C1(7,:)-MATX_PointsDH_C1(3,:))/norm(MATX_PointsDH_C1(7,:)-MATX_PointsDH_C1(3,:))),...
              (MATX_Points_SpecialCase(4,:)-MATX_PointsDH_C1(4,:))/...
                ((MATX_PointsDH_C1(8,:)-MATX_PointsDH_C1(4,:))/norm(MATX_PointsDH_C1(8,:)-MATX_PointsDH_C1(4,:)))]./MATX_ParaStruc(3,:);

StiffC1_Ori=tool.calStiff4R(MATX_PointsDH_C1(1:8,:),[1,4,5,8],1);
StiffC2_Ori=tool.calStiff4R(MATX_PointsDH_C2(1:8,:),[1,4,5,8],1);

p1AVGC1_Ori=1/2*(MATX_PointsDH_C1(1,:)'+MATX_PointsDH_C1(1+4,:)')/100;
p2AVGC1_Ori=1/2*(MATX_PointsDH_C1(2,:)'+MATX_PointsDH_C1(2+4,:)')/100;
p3AVGC1_Ori=1/2*(MATX_PointsDH_C1(3,:)'+MATX_PointsDH_C1(3+4,:)')/100;
p4AVGC1_Ori=1/2*(MATX_PointsDH_C1(4,:)'+MATX_PointsDH_C1(4+4,:)')/100;

p1AVGC2_Ori=1/2*(MATX_PointsDH_C2(1,:)'+MATX_PointsDH_C2(1+4,:)')/100;
p2AVGC2_Ori=1/2*(MATX_PointsDH_C2(2,:)'+MATX_PointsDH_C2(2+4,:)')/100;
p3AVGC2_Ori=1/2*(MATX_PointsDH_C2(3,:)'+MATX_PointsDH_C2(3+4,:)')/100;
p4AVGC2_Ori=1/2*(MATX_PointsDH_C2(4,:)'+MATX_PointsDH_C2(4+4,:)')/100;

arc13C1_Ori=norm(p1AVGC1_Ori-p3AVGC1_Ori);    arc24C1_Ori=norm(p2AVGC1_Ori-p4AVGC1_Ori);
areaH1C1_Ori=calProjectionArea([p1AVGC1_Ori';p2AVGC1_Ori';p3AVGC1_Ori';p4AVGC1_Ori'],MATX_PointsDH_C1(1,:)'-MATX_PointsDH_C1(1+4,:)');
areaH2C1_Ori=calProjectionArea([p1AVGC1_Ori';p2AVGC1_Ori';p3AVGC1_Ori';p4AVGC1_Ori'],MATX_PointsDH_C1(2,:)'-MATX_PointsDH_C1(2+4,:)');
areaH3C1_Ori=calProjectionArea([p1AVGC1_Ori';p2AVGC1_Ori';p3AVGC1_Ori';p4AVGC1_Ori'],MATX_PointsDH_C1(3,:)'-MATX_PointsDH_C1(3+4,:)');
areaH4C1_Ori=calProjectionArea([p1AVGC1_Ori';p2AVGC1_Ori';p3AVGC1_Ori';p4AVGC1_Ori'],MATX_PointsDH_C1(4,:)'-MATX_PointsDH_C1(4+4,:)');

arc13C2_Ori=norm(p1AVGC2_Ori-p3AVGC2_Ori);    arc24C2_Ori=norm(p2AVGC2_Ori-p4AVGC2_Ori);
areaH1C2_Ori=calProjectionArea([p1AVGC2_Ori';p2AVGC2_Ori';p3AVGC2_Ori';p4AVGC2_Ori'],MATX_PointsDH_C2(1,:)'-MATX_PointsDH_C2(1+4,:)');
areaH2C2_Ori=calProjectionArea([p1AVGC2_Ori';p2AVGC2_Ori';p3AVGC2_Ori';p4AVGC2_Ori'],MATX_PointsDH_C2(2,:)'-MATX_PointsDH_C2(2+4,:)');
areaH3C2_Ori=calProjectionArea([p1AVGC2_Ori';p2AVGC2_Ori';p3AVGC2_Ori';p4AVGC2_Ori'],MATX_PointsDH_C2(3,:)'-MATX_PointsDH_C2(3+4,:)');
areaH4C2_Ori=calProjectionArea([p1AVGC2_Ori';p2AVGC2_Ori';p3AVGC2_Ori';p4AVGC2_Ori'],MATX_PointsDH_C2(4,:)'-MATX_PointsDH_C2(4+4,:)');


SlipDist(2:3)=1*[0.3,0.6]; % -1,0,1
NUM_Count=1;
for i=-1.5:0.3:1.5
    for iv=-1.5:0.3:1.5
        SlipDist(1)=i;  SlipDist(4)=iv;

        [pointsDOC1,pointsDPC1,pointsSOC1,pointsSPC1]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                                  MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(4,:));
        [pointsDOC2,pointsDPC2,pointsSOC2,pointsSPC2]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                                  MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(5,:));
%        MATX_PointsDH_C1=[pointsDOC1;pointsDPC1];   MATX_PointsDH_C2=[pointsDOC2;pointsDPC2];
        MATX_Points_C1=[pointsSOC1;pointsSPC1];     MATX_Points_C2=[pointsSOC2;pointsSPC2];

        MATX_FlagContact_C1=ContactDetect4R(MATX_Points_C1);    MATX_FlagContact_C2=ContactDetect4R(MATX_Points_C2);
        FlagContactC1=sum(MATX_FlagContact_C1,"all")/2; FlagContactC2=sum(MATX_FlagContact_C2,"all")/2;

        stiffC1=tool.calStiff4R(MATX_Points_C1,[1,4,5,8],1);    stiffC2=tool.calStiff4R(MATX_Points_C2,[1,4,5,8],1);

        p1AVGC1=1/2*(MATX_Points_C1(1,:)'+MATX_Points_C1(1+4,:)')/100;
        p2AVGC1=1/2*(MATX_Points_C1(2,:)'+MATX_Points_C1(2+4,:)')/100;
        p3AVGC1=1/2*(MATX_Points_C1(3,:)'+MATX_Points_C1(3+4,:)')/100;
        p4AVGC1=1/2*(MATX_Points_C1(4,:)'+MATX_Points_C1(4+4,:)')/100;

        p1AVGC2=1/2*(MATX_Points_C2(1,:)'+MATX_Points_C2(1+4,:)')/100;
        p2AVGC2=1/2*(MATX_Points_C2(2,:)'+MATX_Points_C2(2+4,:)')/100;
        p3AVGC2=1/2*(MATX_Points_C2(3,:)'+MATX_Points_C2(3+4,:)')/100;
        p4AVGC2=1/2*(MATX_Points_C2(4,:)'+MATX_Points_C2(4+4,:)')/100;

        arc13C1=norm(p1AVGC1-p3AVGC1);    arc24C1=norm(p2AVGC1-p4AVGC1);
        areaH1C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(1,:)'-MATX_Points_C1(1+4,:)');
        areaH2C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(2,:)'-MATX_Points_C1(2+4,:)');
        areaH3C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(3,:)'-MATX_Points_C1(3+4,:)');
        areaH4C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(4,:)'-MATX_Points_C1(4+4,:)');

        arc13C2=norm(p1AVGC2-p3AVGC2);    arc24C2=norm(p2AVGC2-p4AVGC2);
        areaH1C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(1,:)'-MATX_Points_C2(1+4,:)');
        areaH2C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(2,:)'-MATX_Points_C2(2+4,:)');
        areaH3C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(3,:)'-MATX_Points_C2(3+4,:)');
        areaH4C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(4,:)'-MATX_Points_C2(4+4,:)');

        MATX_SlipDist(NUM_Count,:)=SlipDist;

        LISTr_ContactC1(NUM_Count)=FlagContactC1;   LISTr_ContactC2(NUM_Count)=FlagContactC2;

        LISTr_StiffC1(NUM_Count)=stiffC1;   LISTr_StiffC2(NUM_Count)=stiffC2;

        LISTr_Arc13C1(NUM_Count)=arc13C1;   LISTr_Arc24C1(NUM_Count)=arc24C1;
        LISTr_AreaH1C1(NUM_Count)=areaH1C1; LISTr_AreaH2C1(NUM_Count)=areaH2C1;
        LISTr_AreaH3C1(NUM_Count)=areaH3C1; LISTr_AreaH4C1(NUM_Count)=areaH4C1;

        LISTr_Arc13C2(NUM_Count)=arc13C2;   LISTr_Arc24C2(NUM_Count)=arc24C2;
        LISTr_AreaH1C2(NUM_Count)=areaH1C2; LISTr_AreaH2C2(NUM_Count)=areaH2C2;
        LISTr_AreaH3C2(NUM_Count)=areaH3C2; LISTr_AreaH4C2(NUM_Count)=areaH4C2;

        NUM_Count=NUM_Count+1;
    end
end

MATX_Info=[MATX_SlipDist,... #4
           LISTr_ContactC1',LISTr_ContactC2',...
           LISTr_StiffC1',LISTr_StiffC2',... #8
           LISTr_Arc13C1',LISTr_Arc24C1',LISTr_Arc13C2',LISTr_Arc24C2',...
           LISTr_AreaH1C1',LISTr_AreaH2C1',LISTr_AreaH3C1',LISTr_AreaH4C1',...#16
           LISTr_AreaH1C2',LISTr_AreaH2C2',LISTr_AreaH3C2',LISTr_AreaH4C2'];
MATX_Info_Compare=[MATX_SlipDist,...
                   LISTr_ContactC1',LISTr_ContactC2',...
                   (LISTr_StiffC1'-StiffC1_Ori)/StiffC1_Ori,(LISTr_StiffC2'-StiffC2_Ori)/StiffC2_Ori,...
                   (LISTr_Arc13C1'-arc13C1_Ori)/arc13C1_Ori,(LISTr_Arc24C1'-arc24C1_Ori)/arc24C1_Ori,...
                   (LISTr_Arc13C2'-arc13C2_Ori)/arc13C2_Ori,(LISTr_Arc24C2'-arc24C2_Ori)/arc24C2_Ori,...
                   (LISTr_AreaH1C1'-areaH1C1_Ori)/areaH1C1_Ori,(LISTr_AreaH2C1'-areaH2C1_Ori)/areaH2C1_Ori,...
                   (LISTr_AreaH3C1'-areaH3C1_Ori)/areaH3C1_Ori,(LISTr_AreaH4C1'-areaH4C1_Ori)/areaH4C1_Ori,...
                   (LISTr_AreaH1C2'-areaH1C2_Ori)/areaH1C2_Ori,(LISTr_AreaH2C2'-areaH2C2_Ori)/areaH2C2_Ori,...
                   (LISTr_AreaH3C2'-areaH3C2_Ori)/areaH3C2_Ori,(LISTr_AreaH4C2'-areaH4C2_Ori)/areaH4C2_Ori];

papercolormap=tool.customcolormap(linspace(0,1,7), {'#AA1E29','#DC4638','#FDB26E','#E7F1D6','#7AB0D6','#3D4DA2','#3E51A1'});

figure(1)
plotStiffC1=scatter(MATX_Info_Compare(:,1),MATX_Info_Compare(:,4),200,MATX_Info_Compare(:,7),'filled');
plotStiffC1.MarkerEdgeColor='k';
plotStiffC1.LineWidth=0.25;
xlim([-2,2]);   ylim([-2,2]);
colormap(papercolormap);
cb=colorbar;
caxis([-1 1]);
cb.Ticks = [-1 -0.5 0 0.5 1];
axis equal
for i=1:size(MATX_Info_Compare,1)
    hold on;
    if MATX_Info_Compare(i,5)+MATX_Info_Compare(i,6)>0
        plotStiffC1_Contact=scatter(MATX_Info_Compare(i,1),MATX_Info_Compare(i,4),500,'marker','s');
        plotStiffC1_Contact.MarkerEdgeColor='none';
        plotStiffC1_Contact.MarkerFaceColor='w'; plotStiffC1_Contact.MarkerFaceAlpha=0.6;
        plotStiffC1.LineWidth=0.5;
        xlim([-2,2]);   ylim([-2,2])
    end
end

figure(2)
plotStiffC2=scatter(MATX_Info_Compare(:,1),MATX_Info_Compare(:,4),200,MATX_Info_Compare(:,8),'filled');
plotStiffC2.MarkerEdgeColor='k';
plotStiffC2.LineWidth=0.25;
xlim([-2,2]);   ylim([-2,2]);
colormap(papercolormap);
cb=colorbar;
caxis([-1 1]);
cb.Ticks = [-1 -0.5 0 0.5 1];
axis equal
for i=1:size(MATX_Info_Compare,1)
    hold on;
    if MATX_Info_Compare(i,5)+MATX_Info_Compare(i,6)>0
        plotStiffC2_Contact=scatter(MATX_Info_Compare(i,1),MATX_Info_Compare(i,4),500,'marker','s');
        plotStiffC2_Contact.MarkerEdgeColor='none';
        plotStiffC2_Contact.MarkerFaceColor='w'; plotStiffC2_Contact.MarkerFaceAlpha=0.6;
        plotStiffC2.LineWidth=0.5;
        xlim([-2,2]);   ylim([-2,2])
    end
end

figure(3)
plotAreaC1=scatter(MATX_Info_Compare(:,1),MATX_Info_Compare(:,4),200,MATX_Info_Compare(:,15),'filled');
plotAreaC1.MarkerEdgeColor='k';
plotAreaC1.LineWidth=0.25;
xlim([-2,2]);   ylim([-2,2]);
colormap(papercolormap);
cb=colorbar;
caxis([-2 2]);
cb.Ticks = [-2 -1 0 1 2];
axis equal
for i=1:size(MATX_Info_Compare,1)
    hold on;
    if MATX_Info_Compare(i,5)+MATX_Info_Compare(i,6)>0
        plotAreaC1_Contact=scatter(MATX_Info_Compare(i,1),MATX_Info_Compare(i,4),500,'marker','s');
        plotAreaC1_Contact.MarkerEdgeColor='none';
        plotAreaC1_Contact.MarkerFaceColor='w'; plotAreaC1_Contact.MarkerFaceAlpha=0.6;
        plotAreaC1.LineWidth=0.5;
        xlim([-2,2]);   ylim([-2,2])
    end
end

figure(4)
plotAreaC2=scatter(MATX_Info_Compare(:,1),MATX_Info_Compare(:,4),200,MATX_Info_Compare(:,19),'filled');
plotAreaC2.MarkerEdgeColor='k';
plotAreaC2.LineWidth=0.25;
xlim([-2,2]);   ylim([-2,2]);
colormap(papercolormap);
cb=colorbar;
caxis([-2 2]);
cb.Ticks = [-2 -1 0 1 2];
axis equal
for i=1:size(MATX_Info_Compare,1)
    hold on;
    if MATX_Info_Compare(i,5)+MATX_Info_Compare(i,6)>0
        plotAreaC2_Contact=scatter(MATX_Info_Compare(i,1),MATX_Info_Compare(i,4),500,'marker','s');
        plotAreaC2_Contact.MarkerEdgeColor='none';
        plotAreaC2_Contact.MarkerFaceColor='w'; plotAreaC2_Contact.MarkerFaceAlpha=0.6;
        plotAreaC2.LineWidth=0.5;
        xlim([-2,2]);   ylim([-2,2])
    end
end

% Pareto analysis
NUM_Count1=1;
NUM_Count2=1;
for i=-1.5:0.5:1.5
    for ii=-1.5:0.5:1.5
        for iii=-1.5:0.5:1.5
            for iv=-1.5:0.5:1.5
                SlipDist(1)=i;  SlipDist(2)=ii;  SlipDist(3)=iii;  SlipDist(4)=iv;

                [pointsDOC1,pointsDPC1,pointsSOC1,pointsSPC1]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                                          MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(4,:));
                [pointsDOC2,pointsDPC2,pointsSOC2,pointsSPC2]=form4RStruc(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                                                                          MATX_ParaStruc(3,:),SlipDist,MATX_ParaStruc(5,:));
                MATX_Points_C1=[pointsSOC1;pointsSPC1];     MATX_Points_C2=[pointsSOC2;pointsSPC2];

                MATX_FlagContact_C1=ContactDetect4R(MATX_Points_C1);    MATX_FlagContact_C2=ContactDetect4R(MATX_Points_C2);
                FlagContactC1=sum(MATX_FlagContact_C1,"all")/2; FlagContactC2=sum(MATX_FlagContact_C2,"all")/2;

                stiffC1=tool.calStiff4R(MATX_Points_C1,[1,4,5,8],1);    stiffC2=tool.calStiff4R(MATX_Points_C2,[1,4,5,8],1);

                p1AVGC1=1/2*(MATX_Points_C1(1,:)'+MATX_Points_C1(1+4,:)')/100;
                p2AVGC1=1/2*(MATX_Points_C1(2,:)'+MATX_Points_C1(2+4,:)')/100;
                p3AVGC1=1/2*(MATX_Points_C1(3,:)'+MATX_Points_C1(3+4,:)')/100;
                p4AVGC1=1/2*(MATX_Points_C1(4,:)'+MATX_Points_C1(4+4,:)')/100;

                p1AVGC2=1/2*(MATX_Points_C2(1,:)'+MATX_Points_C2(1+4,:)')/100;
                p2AVGC2=1/2*(MATX_Points_C2(2,:)'+MATX_Points_C2(2+4,:)')/100;
                p3AVGC2=1/2*(MATX_Points_C2(3,:)'+MATX_Points_C2(3+4,:)')/100;
                p4AVGC2=1/2*(MATX_Points_C2(4,:)'+MATX_Points_C2(4+4,:)')/100;

                arc13C1=norm(p1AVGC1-p3AVGC1);    arc24C1=norm(p2AVGC1-p4AVGC1);
                areaH1C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(1,:)'-MATX_Points_C1(1+4,:)');
                areaH2C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(2,:)'-MATX_Points_C1(2+4,:)');
                areaH3C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(3,:)'-MATX_Points_C1(3+4,:)');
                areaH4C1=calProjectionArea([p1AVGC1';p2AVGC1';p3AVGC1';p4AVGC1'],MATX_Points_C1(4,:)'-MATX_Points_C1(4+4,:)');

                arc13C2=norm(p1AVGC2-p3AVGC2);    arc24C2=norm(p2AVGC2-p4AVGC2);
                areaH1C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(1,:)'-MATX_Points_C2(1+4,:)');
                areaH2C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(2,:)'-MATX_Points_C2(2+4,:)');
                areaH3C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(3,:)'-MATX_Points_C2(3+4,:)');
                areaH4C2=calProjectionArea([p1AVGC2';p2AVGC2';p3AVGC2';p4AVGC2'],MATX_Points_C2(4,:)'-MATX_Points_C2(4+4,:)');

                MATX_SlipDist_all(NUM_Count2,:)=SlipDist;

                LISTr_ContactC1_all(NUM_Count2)=FlagContactC1;   LISTr_ContactC2_all(NUM_Count2)=FlagContactC2;

                LISTr_StiffC1_all(NUM_Count2)=stiffC1;   LISTr_StiffC2_all(NUM_Count2)=stiffC2;

                LISTr_Arc13C1_all(NUM_Count2)=arc13C1;   LISTr_Arc24C1_all(NUM_Count2)=arc24C1;
                LISTr_AreaH1C1_all(NUM_Count2)=areaH1C1; LISTr_AreaH2C1_all(NUM_Count2)=areaH2C1;
                LISTr_AreaH3C1_all(NUM_Count2)=areaH3C1; LISTr_AreaH4C1_all(NUM_Count2)=areaH4C1;

                LISTr_Arc13C2_all(NUM_Count2)=arc13C2;   LISTr_Arc24C2_all(NUM_Count2)=arc24C2;
                LISTr_AreaH1C2_all(NUM_Count2)=areaH1C2; LISTr_AreaH2C2_all(NUM_Count2)=areaH2C2;
                LISTr_AreaH3C2_all(NUM_Count2)=areaH3C2; LISTr_AreaH4C2_all(NUM_Count2)=areaH4C2;

                NUM_Count2=NUM_Count2+1;

                if FlagContactC1+FlagContactC2==0
                    MATX_SlipDist_Pareto(NUM_Count1,:)=SlipDist;

                    LISTr_ContactC1_Pareto(NUM_Count1)=FlagContactC1;   LISTr_ContactC2_Pareto(NUM_Count1)=FlagContactC2;

                    LISTr_StiffC1_Pareto(NUM_Count1)=stiffC1;   LISTr_StiffC2_Pareto(NUM_Count1)=stiffC2;

                    LISTr_Arc13C1_Pareto(NUM_Count1)=arc13C1;   LISTr_Arc24C1_Pareto(NUM_Count1)=arc24C1;
                    LISTr_AreaH1C1_Pareto(NUM_Count1)=areaH1C1; LISTr_AreaH2C1_Pareto(NUM_Count1)=areaH2C1;
                    LISTr_AreaH3C1_Pareto(NUM_Count1)=areaH3C1; LISTr_AreaH4C1_Pareto(NUM_Count1)=areaH4C1;

                    LISTr_Arc13C2_Pareto(NUM_Count1)=arc13C2;   LISTr_Arc24C2_Pareto(NUM_Count1)=arc24C2;
                    LISTr_AreaH1C2_Pareto(NUM_Count1)=areaH1C2; LISTr_AreaH2C2_Pareto(NUM_Count1)=areaH2C2;
                    LISTr_AreaH3C2_Pareto(NUM_Count1)=areaH3C2; LISTr_AreaH4C2_Pareto(NUM_Count1)=areaH4C2;

                    NUM_Count1=NUM_Count1+1;
                end
            end
        end
    end
end

papercolormapRed =tool.customcolormap(linspace(0,1,4),{'#AA1E29','#DC4638','#FDB26E','#E7F1D6'});
papercolormapBlue=tool.customcolormap(linspace(0,1,4),{'#3E51A1','#3D4DA2','#7AB0D6','#E7F1D6'});

MATX_Info_all=[MATX_SlipDist_all,... #4
               LISTr_ContactC1_all',LISTr_ContactC2_all',...
               LISTr_StiffC1_all',LISTr_StiffC2_all',... #8
               LISTr_Arc13C1_all',LISTr_Arc24C1_all',LISTr_Arc13C2_all',LISTr_Arc24C2_all',...
               LISTr_AreaH1C1_all',LISTr_AreaH2C1_all',LISTr_AreaH3C1_all',LISTr_AreaH4C1_all',...#16
               LISTr_AreaH1C2_all',LISTr_AreaH2C2_all',LISTr_AreaH3C2_all',LISTr_AreaH4C2_all'];
MATX_Info_Pareto=[MATX_SlipDist_Pareto,... #4
                  LISTr_ContactC1_Pareto',LISTr_ContactC2_Pareto',...
                  LISTr_StiffC1_Pareto',LISTr_StiffC2_Pareto',... #8
                  LISTr_Arc13C1_Pareto',LISTr_Arc24C1_Pareto',LISTr_Arc13C2_Pareto',LISTr_Arc24C2_Pareto',...
                  LISTr_AreaH1C1_Pareto',LISTr_AreaH2C1_Pareto',LISTr_AreaH3C1_Pareto',LISTr_AreaH4C1_Pareto',...#16
                  LISTr_AreaH1C2_Pareto',LISTr_AreaH2C2_Pareto',LISTr_AreaH3C2_Pareto',LISTr_AreaH4C2_Pareto'];
figure(5)
plotArea=scatter(LISTr_AreaH3C1_all-LISTr_AreaH3C2_all,(LISTr_StiffC1_all'+LISTr_StiffC2_all')/2,...
                  100,sum(abs(MATX_SlipDist_all),2),'filled');hold on
plotArea.MarkerEdgeColor='k'; plotArea.MarkerFaceAlpha=0.8;
plotArea.LineWidth=0.25; plotArea.MarkerEdgeAlpha=0.1;
ylim([0,0.07]);   xlim([-0.02,0.04]);
colormap(papercolormapBlue);
cb=colorbar;
caxis([1 6]);
cb.Ticks = [1,2,3,4,5,6];
hold on

figure(6)
plotPareto=scatter(LISTr_AreaH3C1_Pareto-LISTr_AreaH3C2_Pareto,(LISTr_StiffC1_Pareto'+LISTr_StiffC2_Pareto')/2,...
                    100,sum(abs(MATX_SlipDist_Pareto),2),'filled');
plotPareto.MarkerEdgeColor='k'; plotPareto.MarkerFaceAlpha=0.8;
plotPareto.LineWidth=0.25; plotPareto.MarkerEdgeAlpha=0.1;
ylim([0,0.07]);   xlim([-0.02,0.04]);
colormap(papercolormapRed);
cb=colorbar;
caxis([1 6]);
cb.Ticks = [1,2,3,4,5,6];
hold on




function [pointsDO,pointsDP,pointsSO,pointsSP]=form4RStruc(lVars,alphaVars,offsetVars,...
                                                           SlipDist,thetaVars)

    l12=lVars(1);    l23=lVars(2);    l34=lVars(3);   l41=lVars(4);
    alpha12=alphaVars(1);   alpha23=alphaVars(2);  alpha34=alphaVars(3);
    alpha41=alphaVars(4);
    offset1=offsetVars(1);  offset2=offsetVars(2);   offset3=offsetVars(3);
    offset4=offsetVars(4);
    hingeEdgeO1=SlipDist(1)*offset1;   hingeEdgeP1=SlipDist(1)*offset1;
    hingeEdgeO2=SlipDist(2)*offset2;   hingeEdgeP2=SlipDist(2)*offset2;
    hingeEdgeO3=SlipDist(3)*offset3;   hingeEdgeP3=SlipDist(3)*offset3;
    hingeEdgeO4=SlipDist(4)*offset4;   hingeEdgeP4=SlipDist(4)*offset4;
    theta1=thetaVars(1);  theta2=thetaVars(2);  theta3=thetaVars(3);
    theta4=thetaVars(4);

    % P1
    O1=[0;0;0];
    sys1X=[1;0;0]; sys1Z=[0;0;1];
    p1=O1+offset1*sys1Z;
    vec12=rotawithNormVec(sys1X,sys1Z,theta1);
    O1H=O1+hingeEdgeO1*sys1Z;    p1H=p1+hingeEdgeP1*sys1Z;

    % P2
    O2=p1+l12*vec12;
    sys2X=vec12;   sys2Z=rotawithNormVec(sys1Z,sys2X,alpha12);
    p2=O2+offset2*sys2Z;
    vec23=rotawithNormVec(sys2X,sys2Z,theta2);
    O2H=O2+hingeEdgeO2*sys2Z;    p2H=p2+hingeEdgeP2*sys2Z;

    % P3
    O3=p2+l23*vec23;
    sys3X=vec23;   sys3Z=rotawithNormVec(sys2Z,sys3X,alpha23);
    p3=O3+offset3*sys3Z;
    vec34=rotawithNormVec(sys3X,sys3Z,theta3);
    O3H=O3+hingeEdgeO3*sys3Z;    p3H=p3+hingeEdgeP3*sys3Z;

    % P4
    O4=p3+l34*vec34;
    sys4X=vec34;   sys4Z=rotawithNormVec(sys3Z,sys4X,alpha34);
    p4=O4+offset4*sys4Z;
    vec41=rotawithNormVec(sys4X,sys4Z,theta4);
    O4H=O4+hingeEdgeO4*sys4Z;    p4H=p4+hingeEdgeP4*sys4Z;

    pointsDO=[O1';O2';O3';O4'];
    pointsDP=[p1';p2';p3';p4'];
    pointsSO=[O1H';O2H';O3H';O4H'];
    pointsSP=[p1H';p2H';p3H';p4H'];
end

function vec_=rotawithNormVec(vec,vecN,theta)
    if abs(vec'*vecN)>1e-12
        vec_=[0;0;0];
        error('ERROR in rotawithNormVec: Vectors NOT Orthogonal')
    else
        vecN=vecN/norm(vecN);
        vec_=cos(theta)*vec+(1-cos(theta))*(vecN'*vec)*vecN+sin(theta)*cross(vecN,vec);
        vec_=vec_/norm(vec_);
    end
end

function area=calProjectionArea(NodesList,NormVec)
    n=NormVec/norm(NormVec);
    A=NodesList(1,:)';
    B=NodesList(2,:)'; C=NodesList(3,:)'; D=NodesList(4,:)';

    % Step 2: Project the vertices onto the plane
    proj_A = A - dot(A, n) * n;
    proj_B = B - dot(B, n) * n;
    proj_C = C - dot(C, n) * n;
    proj_D = D - dot(D, n) * n;

    % Step 3: Calculate the projected area using cross product
    % Calculate the vectors between projected points to form two sides of the quadrilateral
    side1 = proj_B - proj_A;
    side2 = proj_D - proj_A;

    % Calculate the cross product to get the area vector
    area = norm(cross(side1, side2));
end

function MATX_FlagContact=ContactDetect4R(MATX_PointsOrigin)
    tool=TOOLS_SharedFunction;
    topoTetra=[1,5,2,6;3,7,2,6;3,7,4,8;1,5,4,8];

    IndexHinge1=[1,2];  IndexHinge2=[3,4];
    for i=1:size(topoTetra,1)
        TSOR_PointsTetra(:,:,i)=[MATX_PointsOrigin(topoTetra(i,1),:);...
                                 MATX_PointsOrigin(topoTetra(i,2),:);...
                                 MATX_PointsOrigin(topoTetra(i,3),:);...
                                 MATX_PointsOrigin(topoTetra(i,4),:)];
        TSOR_PointsTetraStretch(:,:,i)=StretchTetra(TSOR_PointsTetra(:,:,i),IndexHinge1,IndexHinge2);
        CELL_StrucTetraStretch{i}=NodesForm2StrucForm(TSOR_PointsTetraStretch(:,:,i));
    end

    INDX_FlagContact=[1,2;1,3;1,4;2,3;2,4;3,4];
    for i=1:length(INDX_FlagContact)
        FLAG_ContactTetra(i)=tool.GJK(CELL_StrucTetraStretch{INDX_FlagContact(i,1)},...
                                      CELL_StrucTetraStretch{INDX_FlagContact(i,2)},6);
    end

    MATX_FlagContact=zeros([4,4]);
    for i=1:length(FLAG_ContactTetra)
        if FLAG_ContactTetra(i)~=0
            MATX_FlagContact(INDX_FlagContact(i,1),INDX_FlagContact(i,2))=1;
            MATX_FlagContact(INDX_FlagContact(i,2),INDX_FlagContact(i,1))=1;
        end
    end
end

function NodesList_NEW=StretchTetra(NodesList,IndexHinge1,IndexHinge2)
    HingeCenter1=1/2*(NodesList(IndexHinge1(1),:)'+NodesList(IndexHinge1(2),:)');
    HingeCenter2=1/2*(NodesList(IndexHinge2(1),:)'+NodesList(IndexHinge2(2),:)');

    MoveVec12=HingeCenter2-HingeCenter1;    MoveVec21=-MoveVec12;

    MoveDist=1/3000*norm(MoveVec12);

    NodesList_NEW(IndexHinge1,:)=NodesList(IndexHinge1,:)+MoveDist*MoveVec12';
    NodesList_NEW(IndexHinge2,:)=NodesList(IndexHinge2,:)+MoveDist*MoveVec21';

%    NodesList_NEW=NodesList;
end

function StrucForm=NodesForm2StrucForm(NodesList)
    S.XData=NodesList(:,1);
    S.YData=NodesList(:,2);
    S.ZData=NodesList(:,3);
%    S.Faces=topoFacets;
%    StrucForm=patch(S);
    StrucForm=S;
end