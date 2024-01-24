clear;clc;close all;
tool=TOOLS_SharedFunction();

papercolormap=tool.customcolormap([0 0.25 0.5 0.75 1], {'#DC4638','#FDB26E','#E7F1D6','#7AB0D6','#3D4DA2'});

INDX_FixedNodes=[3,4,7,8];
topoBarsStruc=[1,5;2,6;3,7;4,8;1,2;1,6;5,2;5,6;...
               2,3;2,7;3,6;6,7;3,4;3,8;4,7;7,8;...
               1,4;4,5;1,8;5,8];
topoTetra=[1,5,2,6;2,6,3,7;3,7,4,8;4,8,1,5];

MATX_Points_Motion=readmatrix('MATX_Points_Motion_Bistable4R_FixedH3H4.csv');
j=1;
for i=1:8:size(MATX_Points_Motion,1)
    TSOR_Points_Motion(:,:,j)=MATX_Points_Motion(i:i+7,:);
    j=j+1;
end

for i=1:size(TSOR_Points_Motion,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_Points_Motion(:,:,i));
    theta3=-MATX_DHStruc(4,3);
    LISTr_theta3(i)=theta3;
end

LISTr_Lengths_Ori=Lengths_Bars(TSOR_Points_Motion(:,:,1),topoBarsStruc);

INDX_Marker=round(linspace(1,102,5));

for i=1:length(INDX_Marker)
    MATX_Points_Motion=TSOR_Points_Motion(:,:,INDX_Marker(i));
    LISTr_Lengths=Lengths_Bars(MATX_Points_Motion,topoBarsStruc);
    LISTr_DeltaLength=(LISTr_Lengths-LISTr_Lengths_Ori)./LISTr_Lengths_Ori;

    MATX_DeltaLength(i,:)=LISTr_DeltaLength;
end
VALE_MaxElement=max(MATX_DeltaLength(:));   VALE_MinElement=min(MATX_DeltaLength(:));
VALE_Range=0.1;%max([abs(VALE_MaxElement),abs(VALE_MinElement)]);
LISTr_IntervalDeltaLength=linspace(-VALE_Range,VALE_Range,64);

for i=1:length(INDX_Marker)
    MATX_Points_Motion=TSOR_Points_Motion(:,:,INDX_Marker(i));
    LISTr_Lengths=Lengths_Bars(MATX_Points_Motion,topoBarsStruc);
    LISTr_DeltaLength=(LISTr_Lengths-LISTr_Lengths_Ori)./LISTr_Lengths_Ori;

    f(i)=figure(i);clf;
    for ii=1:length(LISTr_DeltaLength)
        VALE_DeltaLength=LISTr_DeltaLength(ii);
        LISTr_IntervalDeltaLength_temp=sort([LISTr_IntervalDeltaLength,VALE_DeltaLength]);
        INDX_DeltaLength=find(LISTr_IntervalDeltaLength_temp==VALE_DeltaLength); INDX_DeltaLength=INDX_DeltaLength(1);
        if (INDX_DeltaLength==1) || (INDX_DeltaLength==65)
            MATX_DeltaLength_ColorIndex(ii,:)=[LISTr_DeltaLength(ii),INDX_DeltaLength,INDX_DeltaLength];
        else
            MATX_DeltaLength_ColorIndex(ii,:)=[LISTr_DeltaLength(ii),INDX_DeltaLength-1,INDX_DeltaLength];
        end

        colorSelec(ii,:)=papercolormap(INDX_DeltaLength,:);
        p(i)=plot3(MATX_Points_Motion(topoBarsStruc(ii,:),1),...
                MATX_Points_Motion(topoBarsStruc(ii,:),2),...
                MATX_Points_Motion(topoBarsStruc(ii,:),3));
        p(i).LineWidth=5; p(i).Color=colorSelec(ii,:);
        hold on
        for iii=1:size(MATX_Points_Motion,1)
            p(i)=plot3(MATX_Points_Motion(iii,1),MATX_Points_Motion(iii,2),MATX_Points_Motion(iii,3),'ko');
            p(i).MarkerFaceColor='#E6E7E8'; p(i).MarkerSize=10;
            hold on
        end
    end
    tetramesh(topoTetra,MATX_Points_Motion,'FaceColor','#E6E7E8','FaceAlpha',0.2,'EdgeColor','none');
    axis equal
    view(40,70)
    colormap(papercolormap);    colorbar;
    caxis([-VALE_Range,VALE_Range]);
    f(i).Renderer='Painters';

    TSOR_DeltaLength_ColorIndex(:,:,i)=[MATX_DeltaLength_ColorIndex(:,1),colorSelec];
    TSOR_Points_ColorIndex(:,:,i)=MATX_Points_Motion;
end



function Lengths=Lengths_Bars(Points,Bars)
    [row_bars,~]=size(Bars);
    Lengths=[];
    for i=1:row_bars
        length{i}=sqrt((Points(Bars(i,1),1)-Points(Bars(i,2),1))^2+...
                       (Points(Bars(i,1),2)-Points(Bars(i,2),2))^2+...
                       (Points(Bars(i,1),3)-Points(Bars(i,2),3))^2);
        Lengths=cat(1,Lengths,length{i});
    end
end