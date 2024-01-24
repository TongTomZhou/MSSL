clear;clc;close all;
tool=TOOLS_SharedFunction();

INDX_FixedNodes=[1,4,5,8];
topoBarsStruc=[1,5;2,6;3,7;4,8;1,2;1,6;5,2;5,6;...
               2,3;2,7;3,6;6,7;3,4;3,8;4,7;7,8;...
               1,4;4,5;1,8;5,8];

MATX_Points_Motion=readmatrix('MATX_Points_Motion_Bistable4R.csv');
j=1;
for i=1:8:size(MATX_Points_Motion,1)
    TSOR_Points_Motion(:,:,j)=MATX_Points_Motion(i:i+7,:);
    j=j+1;
end
LISTv_Lengths_Ori=Lengths_Bars(TSOR_Points_Motion(:,:,1),topoBarsStruc);

for i=1:size(TSOR_Points_Motion,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_Points_Motion(:,:,i));
    theta1=MATX_DHStruc(4,1); theta2=MATX_DHStruc(4,2);theta3=-MATX_DHStruc(4,3); theta4=MATX_DHStruc(4,4);
    if theta1<-0
        theta1=theta1+2*pi;
    elseif theta1>pi
        theta1=theta1-2*pi;
    end
    if theta2<-0
        theta2=theta2+2*pi;
    elseif theta2>pi
        theta2=theta2-2*pi;
    end
%    if theta3<-0
%        theta3=theta3+2*pi;
%    elseif theta3>pi
%        theta3=theta3-2*pi;
%    end
    if theta4<-0
        theta4=theta4+2*pi;
    elseif theta4>pi
        theta4=theta4-2*pi;
    end

    LISTr_Theta1(i)=theta1; LISTr_Theta2(i)=theta2; LISTr_Theta3(i)=theta3; LISTr_Theta4(i)=theta4;

    % Energy analysis
    LISTv_Lengths_i=Lengths_Bars(TSOR_Points_Motion(:,:,i),topoBarsStruc);
    LISTv_DeltaLengths_i=0;
    for ii=1:length(LISTv_Lengths_i)
        LISTv_DeltaLengths_i=LISTv_DeltaLengths_i+(LISTv_Lengths_i(ii)/100-LISTv_Lengths_Ori(ii)/100)^2;
    end
    LISTv_DeltaLengths(i)=LISTv_DeltaLengths_i;

    % Shape analysis
    p1AVG=1/2*(TSOR_Points_Motion(1,:,i)'+TSOR_Points_Motion(1+4,:,i)')/100;
    p2AVG=1/2*(TSOR_Points_Motion(2,:,i)'+TSOR_Points_Motion(2+4,:,i)')/100;
    p3AVG=1/2*(TSOR_Points_Motion(3,:,i)'+TSOR_Points_Motion(3+4,:,i)')/100;
    p4AVG=1/2*(TSOR_Points_Motion(4,:,i)'+TSOR_Points_Motion(4+4,:,i)')/100;

    arc13=norm(p1AVG-p3AVG);    arc24=norm(p2AVG-p4AVG);
    areaH1=calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],TSOR_Points_Motion(1,:,i)'-TSOR_Points_Motion(1+4,:,i)');
    areaH2=calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],TSOR_Points_Motion(2,:,i)'-TSOR_Points_Motion(2+4,:,i)');
    areaH3=calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],TSOR_Points_Motion(3,:,i)'-TSOR_Points_Motion(3+4,:,i)');
    areaH4=calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],TSOR_Points_Motion(4,:,i)'-TSOR_Points_Motion(4+4,:,i)');
    
    vec31=p1AVG-p3AVG; vec24=p4AVG-p2AVG;
    vec31Norm=vec31/norm(vec31);    vec24Norm=vec24/norm(vec24);
    ang_31_24=acos(vec31Norm'*vec24Norm);
    
    LISTv_Arc13(i)=arc13;   LISTv_Arc24(i)=arc24;
    LISTv_AreaH1(i)=areaH1; LISTv_AreaH2(i)=areaH2; LISTv_AreaH3(i)=areaH3; LISTv_AreaH4(i)=areaH4;
    LISTv_Ang_31_24(i)=ang_31_24;
end

for i=1:size(LISTv_DeltaLengths,2)-1
    LISTv_DeltaLengths_Slope(i)=(LISTv_DeltaLengths(i+1)-LISTv_DeltaLengths(i))/(LISTr_Theta3(i+1)-LISTr_Theta3(i));
end

figure(1)
plot(LISTr_Theta3/pi*180,LISTv_DeltaLengths,'k--')
xlim([-20,140]);   ylim([0,0.02]);
grid on

figure(2)
plot(LISTr_Theta3(1:end-1)/pi*180,LISTv_DeltaLengths_Slope,'k--')
xlim([-20,140]);   ylim([-0.05,0.04]);
grid on

figure(3)
plot(LISTr_Theta3/pi*180,LISTr_Theta1'/pi*180,'r-'); hold on
plot(LISTr_Theta3/pi*180,LISTr_Theta2'/pi*180,'g-'); hold on
plot(LISTr_Theta3/pi*180,LISTr_Theta4'/pi*180,'b-'); hold on
xlim([-20,140]);   ylim([0,360]); yticks(0:36:360);

figure(4)
plot(LISTr_Theta3/pi*180,LISTv_Arc13','r-'); hold on
plot(LISTr_Theta3/pi*180,LISTv_Arc24','g-'); hold on
xlim([-20,140]);   ylim([0.1,0.3]); yticks(0.1:0.02:0.3);

figure(5)
plot(LISTr_Theta3/pi*180,LISTv_AreaH1','r-'); hold on
plot(LISTr_Theta3/pi*180,LISTv_AreaH2','g-'); hold on
plot(LISTr_Theta3/pi*180,LISTv_AreaH3','b-'); hold on
plot(LISTr_Theta3/pi*180,LISTv_AreaH4','m-'); hold on
xlim([-20,140]);   ylim([0,0.036]); yticks(0:0.0036:0.036);

figure(6)
plot(LISTr_Theta3/pi*180,LISTv_Ang_31_24,'k--')
xlim([-20,140]);   %ylim([0,0.02]);
grid on

MATX_Info=[LISTr_Theta3',... 
           LISTr_Theta1',LISTr_Theta2',LISTr_Theta4',... 
           LISTv_Arc13',LISTv_Arc24',LISTv_Ang_31_24',... 
           LISTv_AreaH1',LISTv_AreaH2',LISTv_AreaH3',LISTv_AreaH4'];
writematrix(MATX_Info,"MATX_ConfigurationAnalysis_Bistable4R.csv");



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
