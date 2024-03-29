%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Theme:     This is a supplementary code of the paper: Multi-stable spatial linkages.
%%%%            Run it to generate a single-loop bi-stable 4R linkage.
%%%% Author:    Tong Zhou
%%%% Date:      2024.01.01
%%%% Result of the introductory case in Fig. 2:
%%%%            l12=10.399710333559419, l23=10.566539864467707, l34=14.080282366989588, l41=14.260584194007084.
%%%%            alpha12=0.965831304681946, alpha23=1.880977481770469, alpha34=1.281112142627775, alpha41=2.191739778773956.
%%%%            r1=7.231794114984756, r2=6.869831414844362, r3=6.565829382115895, r4=9.109619407442986.
%%%%            d1=-4,      d2=2,        d3=4,        d4=-11.
%%%%            theta1S1=0.016770135538180, theta2S1=3.613980211371010, theta3S1=6.278128970831695, theta4S1=3.625122306986279.
%%%%            theta1S2=3.158099253657525, theta2S2=1.107786207703004, theta3S2=4.227727051568603, theta4S2=5.188510172458066.
%%%% Note:      As the initial guess is randomly selected, the bi-stable 4R linkage generated by this code
%%%%            might be distinct from the result in Fig. 2 (i.e. the result provided above).
%%%%            It is still possible to review the result in Fig. 2 by substituting it into the constraints in the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all; clc;
format short g

VALE_ThetaTargetC1=0*pi; VALE_ThetaInterval=1*pi;
LISTr_ThetaTarget=[VALE_ThetaTargetC1,VALE_ThetaTargetC1+VALE_ThetaInterval];
INDX_FixedNodes=[1,4,5,8]; VALE_InteraMax=100000;

%% Initial Guess
MATX_ParaStrucInit=standardOCM(LISTr_ThetaTarget(1),LISTr_ThetaTarget(2));
MATX_ParaStrucInit=[MATX_ParaStrucInit;zeros(1,4)];

%% Generate 4R Tristable
A=[];b=[];  Aeq=[];beq=[];  lb=[];ub=[];
options=optimoptions('fmincon',... %  'ConstraintTolerance',1e-12,...
                     'MaxFunctionEvaluations',VALE_InteraMax,...
                     'MaxIterations',VALE_InteraMax);
LISTr_ParaStrucInit=reshape(MATX_ParaStrucInit,[],1);

obj4R=@(vars)objFunc4R(vars);
nonlcon4R=@(vars)conFunc4R(vars,LISTr_ThetaTarget);

[LISTv_Vars4R,fval4R,exitflag4R]=...
    fmincon(obj4R,LISTr_ParaStrucInit,A,b,Aeq,beq,lb,ub,nonlcon4R,options);

MATX_ParaStruc=reshape(LISTv_Vars4R,6,4);

%% Post-Analysis
[PointsDOC1,PointsDPC1,PointsSOC1,PointsSPC1]=...
        form4RStruc4(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                          MATX_ParaStruc(3,:),MATX_ParaStruc(6,:),MATX_ParaStruc(4,:));
[PointsDOC2,PointsDPC2,PointsSOC2,PointsSPC2]=...
        form4RStruc4(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),...
                          MATX_ParaStruc(3,:),MATX_ParaStruc(6,:),MATX_ParaStruc(5,:));

MATX_PointsD(:,:,1)=[PointsDOC1;PointsDPC1];    MATX_PointsD(:,:,2)=[PointsDOC2;PointsDPC2];
MATX_PointsS(:,:,1)=[PointsSOC1;PointsSPC1];    MATX_PointsS(:,:,2)=[PointsSOC2;PointsSPC2];

VALE_StiffC1=calStiff4R(MATX_PointsS(:,:,1),INDX_FixedNodes,1);
VALE_StiffC2=calStiff4R(MATX_PointsS(:,:,2),INDX_FixedNodes,1);
VALE_DeviaTarget=abs(MATX_ParaStruc(4:5,1)-reshape(LISTr_ThetaTarget,[],1));

MATX_SelfContact_C1=ContactDetect4R([PointsSOC1;PointsSPC1]);
MATX_SelfContact_C2=ContactDetect4R([PointsSOC2;PointsSPC2]);
NUMB_SelfContact_C1=sum(MATX_SelfContact_C1,"all")/2;
NUMB_SelfContact_C2=sum(MATX_SelfContact_C2,"all")/2;

subplot(2,2,1)
plot4RDH(MATX_PointsD(:,:,1),'','c');
subplot(2,2,2)
plot4RDH(MATX_PointsD(:,:,2),'','m');
subplot(2,2,3)
plot4RStruc(MATX_PointsS(:,:,1),'','c');
subplot(2,2,4)
plot4RStruc(MATX_PointsS(:,:,2),'','m');

[c,ceq]=conFunc4R(LISTv_Vars4R,LISTr_ThetaTarget)


function obj=objFunc4R(vars)
    FixNodes=[1,4,5,8];
    ParaStruc=reshape(vars,6,4);

    [~,~,PointsSOC1,PointsSPC1]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),...
                                                  ParaStruc(6,:),ParaStruc(4,:));
    [~,~,PointsSOC2,PointsSPC2]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),...
                                                  ParaStruc(6,:),ParaStruc(5,:));

    StiffC1Val=calStiff4R([PointsSOC1;PointsSPC1],FixNodes,1);
    StiffC2Val=calStiff4R([PointsSOC2;PointsSPC2],FixNodes,1);
    obj=-(StiffC1Val+StiffC2Val)/2;

%    % No Contact Constraints
%    MATX_SelfContact_C1=ContactDetect4R([PointsSOC1;PointsSPC1]);
%    MATX_SelfContact_C2=ContactDetect4R([PointsSOC2;PointsSPC2]);
%    NUMB_SelfContact_C1=sum(MATX_SelfContact_C1,"all")/2;
%    NUMB_SelfContact_C2=sum(MATX_SelfContact_C2,"all")/2;
%    obj=NUMB_SelfContact_C1+NUMB_SelfContact_C2;
end

function [c,ceq]=conFunc4R(vars,ThetaTarget)
    c=[];   ceq=[];
    ParaStruc=reshape(vars,6,4);

    [cClosC1,ceqClosC1]=con4RClosure(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),ParaStruc(4,:));
    [cClosC2,ceqClosC2]=con4RClosure(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),ParaStruc(5,:));

    [cShape,ceqShape]=con4RShape(ParaStruc);
    [cTarget,ceqTarget]=con4RTarget(ParaStruc,ThetaTarget);
    [cNoC,ceqNoC]=con4RNoConflict(ParaStruc);

    ceq=[ceqClosC1;ceqClosC2;ceqShape;ceqTarget;ceqNoC];
    c=[cClosC1;cClosC2;cShape;cTarget;cNoC];
end

function [c,ceq]=con4RNoConflict(ParaStruc)
    c=[];   ceq=[];

    [~,~,PointsSOC1,PointsSPC1]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),...
                                             ParaStruc(3,:),ParaStruc(6,:),ParaStruc(4,:));
    [~,~,PointsSOC2,PointsSPC2]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),...
                                             ParaStruc(3,:),ParaStruc(6,:),ParaStruc(5,:));

    MATX_SelfContact_C1=ContactDetect4R([PointsSOC1;PointsSPC1]);
    MATX_SelfContact_C2=ContactDetect4R([PointsSOC2;PointsSPC2]);
    NUMB_SelfContact_C1=sum(MATX_SelfContact_C1,"all")/2;
    NUMB_SelfContact_C2=sum(MATX_SelfContact_C2,"all")/2;

    ceq=[ceq;NUMB_SelfContact_C1+NUMB_SelfContact_C2];
end

function [c,ceq]=con4RTarget(ParaStruc,ThetaTarget)
    c=[];   ceq=[];

    [~,~,PointsSOC1,PointsSPC1]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),...
                                                  ParaStruc(3,:),ParaStruc(6,:),ParaStruc(4,:));
    [~,~,PointsSOC2,PointsSPC2]=form4RStruc4(ParaStruc(1,:),ParaStruc(2,:),...
                                                  ParaStruc(3,:),ParaStruc(6,:),ParaStruc(5,:));
    PointsSC1=1/2*(PointsSOC1+PointsSPC1); PointsSC2=1/2*(PointsSOC2+PointsSPC2);

    theta1C1=sign([0,0,1]*cross([1;0;0],PointsSC1(2,:)'-PointsSC1(1,:)'))*...
             acos(((PointsSC1(2,:)-PointsSC1(1,:))*[1;0;0])/norm(PointsSC1(2,:)-PointsSC1(1,:)));
    theta1C2=sign([0,0,1]*cross([1;0;0],PointsSC2(2,:)'-PointsSC2(1,:)'))*...
             acos(((PointsSC2(2,:)-PointsSC2(1,:))*[1;0;0])/norm(PointsSC2(2,:)-PointsSC2(1,:)));
    ceq=[theta1C1-ThetaTarget(1);theta1C2-ThetaTarget(2)];
end

function [c,ceq]=con4RClosure(lenList,alphaList,offsetList,thetaList)
    c=[];   ceq=[];

    % Closure Condition
    trans21=transMatDH(thetaList(1),alphaList(1),lenList(1),offsetList(1),'DES');
    trans32=transMatDH(thetaList(2),alphaList(2),lenList(2),offsetList(2),'DES');
    trans34=transMatDH(thetaList(3),alphaList(3),lenList(3),offsetList(3),'ASC');
    trans41=transMatDH(thetaList(4),alphaList(4),lenList(4),offsetList(4),'ASC');

    conditionClosure4R=trans21*trans32-trans41*trans34;
    conditionClosure4R=reshape(conditionClosure4R,[],1);
    ceq=[ceq;conditionClosure4R];
end

function [c,ceq]=con4RShape(ParaStruc)
    c=[];   ceq=[];

    lenList=ParaStruc(1,:);  alphaList=ParaStruc(2,:);    offsetList=ParaStruc(3,:);
    thetaListC1=ParaStruc(4,:); thetaListC2=ParaStruc(5,:);

    c=[c;-lenList(1);   -lenList(2);    -lenList(3);    -lenList(4);...
          lenList(1)-20; lenList(2)-20;  lenList(3)-20;  lenList(4)-20];
    c=[c;-alphaList(1)-pi;     -alphaList(2)-pi;      -alphaList(3)-pi;      -alphaList(4)-pi;...
          alphaList(1)-pi; alphaList(2)-pi;  alphaList(3)-pi;  alphaList(4)-pi];
    c=[c;-abs(offsetList(1))+1e-3;-abs(offsetList(2))+1e-3;-abs(offsetList(3))+1e-3;...
         -abs(offsetList(4))+1e-3;...
          abs(offsetList(1))-20;abs(offsetList(2))-20;abs(offsetList(3))-20;...
          abs(offsetList(4))-20];
    c=[c;-thetaListC1(1)-2*pi;-thetaListC1(2)-2*pi;-thetaListC1(3)-2*pi;-thetaListC1(4)-2*pi;...
          thetaListC1(1)-2*pi; thetaListC1(2)-2*pi; thetaListC1(3)-2*pi; thetaListC1(4)-2*pi];
    c=[c;-thetaListC2(1)-2*pi;-thetaListC2(2)-2*pi;-thetaListC2(3)-2*pi;-thetaListC2(4)-2*pi;...
          thetaListC2(1)-2*pi; thetaListC2(2)-2*pi; thetaListC2(3)-2*pi; thetaListC2(4)-2*pi];
end

function dhPara=standardOCM(theta1C1,theta1C2)
    alpha12=pi/4; alpha23=5*pi/6;   % alpha12=pi/4; alpha23=2*pi/3;
    alpha34=alpha12;     alpha41=alpha23;
    l12=4; l23=sin(alpha23)/sin(alpha12)*l12;    l34=l12;   l41=l23;
    offset1=0.001;      offset2=0.001;     offset3=0.001;       offset4=0.001;

    if abs(theta1C1)<1e-5
        theta2C1=pi;
    else
        theta2C1=2*atan(sin(alpha23/2+alpha12/2)/sin(alpha23/2-alpha12/2)/tan(theta1C1/2));
    end
    theta3C1=-theta1C1; theta4C1=-theta2C1;

    if abs(theta1C2)<1e-5
        theta2C2=pi;
    else
        theta2C2=2*atan(sin(alpha23/2+alpha12/2)/sin(alpha23/2-alpha12/2)/tan(theta1C2/2));
    end
    theta3C2=-theta1C2; theta4C2=-theta2C2;

    dhPara=[l12,l23,l34,l41;...
            alpha12,alpha23,alpha34,alpha41;...
            offset1,offset2,offset3,offset4;...
            theta1C1,theta2C1,theta3C1,theta4C1;...
            theta1C2,theta2C2,theta3C2,theta4C2];
end

function MATX_FlagContact=ContactDetect4R(MATX_PointsOrigin)
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
        FLAG_ContactTetra(i)=GJK(CELL_StrucTetraStretch{INDX_FlagContact(i,1)},...
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

function [pointsDO,pointsDP,pointsSO,pointsSP]=form4RStruc4(lVars,alphaVars,offsetVars,...
                                                            slipEdgeVars,thetaVars)

    l12=lVars(1);    l23=lVars(2);    l34=lVars(3);   l41=lVars(4);
    alpha12=alphaVars(1);   alpha23=alphaVars(2);  alpha34=alphaVars(3);
    alpha41=alphaVars(4);
    offset1=offsetVars(1);  offset2=offsetVars(2);   offset3=offsetVars(3);
    offset4=offsetVars(4);
    slipEdgeC1=slipEdgeVars(1);   slipEdgeC2=slipEdgeVars(2);
    slipEdgeC3=slipEdgeVars(3);   slipEdgeC4=slipEdgeVars(4);
    theta1=thetaVars(1);  theta2=thetaVars(2);  theta3=thetaVars(3);
    theta4=thetaVars(4);

    % P1
    O1=[0;0;0];
    sys1X=[1;0;0]; sys1Z=[0;0;1];
    p1=O1+offset1*sys1Z;
    vec12=rotawithNormVec(sys1X,sys1Z,theta1);
    C1=1/2*(O1+p1)+slipEdgeC1*sys1Z;;

    % P2
    O2=p1+l12*vec12;
    sys2X=vec12;   sys2Z=rotawithNormVec(sys1Z,sys2X,alpha12);
    p2=O2+offset2*sys2Z;
    vec23=rotawithNormVec(sys2X,sys2Z,theta2);
    C2=1/2*(O2+p2)+slipEdgeC2*sys2Z;

    % P3
    O3=p2+l23*vec23;
    sys3X=vec23;   sys3Z=rotawithNormVec(sys2Z,sys3X,alpha23);
    p3=O3+offset3*sys3Z;
    vec34=rotawithNormVec(sys3X,sys3Z,theta3);
    C3=1/2*(O3+p3)+slipEdgeC3*sys3Z;

    % P4
    O4=p3+l34*vec34;
    sys4X=vec34;   sys4Z=rotawithNormVec(sys3Z,sys4X,alpha34);
    p4=O4+offset4*sys4Z;
    vec41=rotawithNormVec(sys4X,sys4Z,theta4);
    C4=1/2*(O4+p4)+slipEdgeC4*sys4Z;

%    lenHinge=mean(lVars);
    lenHinge=4; % Follow the standarxdOCM. Modified at 2023.04.02.
    O1H=C1-lenHinge/2*sys1Z;    p1H=C1+lenHinge/2*sys1Z;
    O2H=C2-lenHinge/2*sys2Z;    p2H=C2+lenHinge/2*sys2Z;
    O3H=C3-lenHinge/2*sys3Z;    p3H=C3+lenHinge/2*sys3Z;
    O4H=C4-lenHinge/2*sys4Z;    p4H=C4+lenHinge/2*sys4Z;

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

function eValStiff=calStiff4R(points,indexFixnodes,K)
    topoLink=[1,5;2,6;3,7;4,8;...
              1,2;1,6;5,2;5,6;...
              2,3;2,7;3,6;6,7;...
              3,4;3,8;4,7;7,8;...
              1,4;4,5;1,8;5,8];
    matStiff=calStiffMat(points,topoLink,indexFixnodes,K);
    eValStiff=min(eig(matStiff));
end

function Matrix_K=calStiffMat(points,topoBars,indexFixNodes,constElastic)
    Bars=topoBars;

    K=[];
    for i=1:size(topoBars,1)
%        K=[K;constElastic/norm(points(Bars(i,1),:)-points(Bars(i,2),:))];
        K=[K;constElastic]; % Modified/Delete at 2023.03.29
    end

    points_b=zeros(size(points,1),3);
    for i=1:length(indexFixNodes)
        points_b(indexFixNodes(i),:)=[1,1,1];
    end
    points_all=[points,points_b];
    Matrix_H=eqtr3m4Bars(points_all,Bars);
    Matrix_C=Matrix_H';

    Matrix_K=Matrix_H*diag(K)*Matrix_C; % stiffness matrix

end

function mat=transMatDH(theta,alpha,lenBar,offset,type)
    if type=='DES'
        mat=[cos(theta),-cos(alpha)*sin(theta), sin(alpha)*sin(theta),lenBar*cos(theta);...
             sin(theta), cos(alpha)*cos(theta),-sin(alpha)*cos(theta),lenBar*sin(theta);...
             0,sin(alpha),cos(alpha),offset;...
             0,0,0,1];
    elseif type=='ASC'
        mat=[ cos(theta),            sin(theta),           0,         -lenBar;
             -cos(alpha)*sin(theta), cos(alpha)*cos(theta),sin(alpha),-offset*sin(alpha);...
              sin(alpha)*sin(theta),-sin(alpha)*cos(theta),cos(alpha),-offset*cos(alpha);...
              0,0,0,1];
    end
end

function plot4RDH(points,figureTitle,color)
    offsetDist=norm(points(1,:)-points(2,:))/15;
    topoBars=[1,5;5,2;2,6;6,3;3,7;7,4;4,8;8,1];

    pointsDO=points(1:4,:); pointsDP=points(5:8,:);

    for i=1:size(pointsDO,1)
        plot3(pointsDO(i,1),pointsDO(i,2),pointsDO(i,3),'ko')
        hold on
        text(pointsDO(i,1),pointsDO(i,2),pointsDO(i,3)+offsetDist,string(i))
        hold on
    end

    for i=1:size(pointsDP,1)
        plot3(pointsDP(i,1),pointsDP(i,2),pointsDP(i,3),'k+')
        hold on
    end

    for i=2:2:size(topoBars,1)
        plot3([points(topoBars(i,1),1),points(topoBars(i,2),1)],...
              [points(topoBars(i,1),2),points(topoBars(i,2),2)],...
              [points(topoBars(i,1),3),points(topoBars(i,2),3)],color)
        hold on
    end
    for i=1:2:size(topoBars,1)
        plot3([points(topoBars(i,1),1),points(topoBars(i,2),1)],...
              [points(topoBars(i,1),2),points(topoBars(i,2),2)],...
              [points(topoBars(i,1),3),points(topoBars(i,2),3)],'k--')
        hold on
    end

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    title(figureTitle)
end

function plot4RStruc(points,figureTitle,color)
    topoTetra=[1,5,2,6;2,6,3,7;3,7,4,8;4,8,1,5];

    offsetDist=norm(points(1,:)-points(2,:))/15;

    tetramesh(topoTetra,points,'FaceColor',color,'FaceAlpha',0.2,'EdgeColor','k')
    hold on

    for i=1:size(points,1)/2
        plot3(points(i,1),points(i,2),points(i,3),'ko')
        hold on
        text(points(i,1),points(i,2),points(i,3)+offsetDist,'O'+string(i));
        hold on
    end
    for i=1+size(points,1)/2:size(points,1)
        text(points(i,1),points(i,2),points(i,3)+offsetDist,'P'+string(i-size(points,1)/2));
        hold on
    end

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    title(figureTitle)
end

function flag = GJK(shape1,shape2,iterations)

    %Point 1 and 2 selection (line segment)
    v = [0.8 0.5 1]; % Random intial vector
    [a,b] = pickLine(v,shape2,shape1);

    %Point 3 selection (triangle)
    [a,b,c,flag] = pickTriangle(a,b,shape2,shape1,iterations);

    %Point 4 selection (tetrahedron)
    if flag == 1 %Only bother if we could find a viable triangle.
        [a,b,c,d,flag] = pickTetrahedron(a,b,c,shape2,shape1,iterations);
    end

end

function [a,b] = pickLine(v,shape1,shape2)
    %Construct the first line of the simplex
    b = support(shape2,shape1,v);
    a = support(shape2,shape1,-v);
end

function [a,b,c,flag] = pickTriangle(a,b,shape1,shape2,IterationAllowed)
    flag = 0; %So far, we don't have a successful triangle.

    %First try:
    ab = b-a;
    ao = -a;
    v = cross(cross(ab,ao),ab); % v is perpendicular to ab pointing in the general direction of the origin.

    c = b;
    b = a;
    a = support(shape2,shape1,v);

    for i = 1:IterationAllowed %iterations to see if we can draw a good triangle.
        %Time to check if we got it:
        ab = b-a;
        ao = -a;
        ac = c-a;

        %Normal to face of triangle
        abc = cross(ab,ac);

        %Perpendicular to AB going away from triangle
        abp = cross(ab,abc);
        %Perpendicular to AC going away from triangle
        acp = cross(abc,ac);

        %First, make sure our triangle "contains" the origin in a 2d projection
        %sense.
        %Is origin above (outside) AB?
        if dot(abp,ao) > 0
            c = b; %Throw away the furthest point and grab a new one in the right direction
            b = a;
            v = abp; %cross(cross(ab,ao),ab);

            %Is origin above (outside) AC?
        elseif dot(acp, ao) > 0
            b = a;
            v = acp; %cross(cross(ac,ao),ac);

        else
            flag = 1;
            break; %We got a good one.
        end
        a = support(shape2,shape1,v);
    end
end

function [a,b,c,d,flag] = pickTetrahedron(a,b,c,shape1,shape2,IterationAllowed)
    %Now, if we're here, we have a successful 2D simplex, and we need to check
    %if the origin is inside a successful 3D simplex.
    %So, is the origin above or below the triangle?
    flag = 0;

    ab = b-a;
    ac = c-a;

    %Normal to face of triangle
    abc = cross(ab,ac);
    ao = -a;

    if dot(abc, ao) > 0 %Above
        d = c;
        c = b;
        b = a;

        v = abc;
        a = support(shape2,shape1,v); %Tetrahedron new point

    else %below
        d = b;
        b = a;
        v = -abc;
        a = support(shape2,shape1,v); %Tetrahedron new point
    end

    for i = 1:IterationAllowed %Allowing 10 tries to make a good tetrahedron.
        %Check the tetrahedron:
        ab = b-a;
        ao = -a;
        ac = c-a;
        ad = d-a;

        %We KNOW that the origin is not under the base of the tetrahedron based on
        %the way we picked a. So we need to check faces ABC, ABD, and ACD.

        %Normal to face of triangle
        abc = cross(ab,ac);

        if dot(abc, ao) > 0 %Above triangle ABC
            %No need to change anything, we'll just iterate again with this face as
            %default.
        else
            acd = cross(ac,ad);%Normal to face of triangle

            if dot(acd, ao) > 0 %Above triangle ACD
                %Make this the new base triangle.
                b = c;
                c = d;
                ab = ac;
                ac = ad;
                abc = acd;
            elseif dot(acd, ao) < 0
                adb = cross(ad,ab);%Normal to face of triangle

                if dot(adb, ao) > 0 %Above triangle ADB
                    %Make this the new base triangle.
                    c = b;
                    b = d;
                    ac = ab;
                    ab = ad;
                    abc = adb;
                else
                    flag = 1;
                    break; %It's inside the tetrahedron.
                end
            end
        end

        %try again:
        if dot(abc, ao) > 0 %Above
            d = c;
            c = b;
            b = a;
            v = abc;
            a = support(shape2,shape1,v); %Tetrahedron new point
        else %below
            d = b;
            b = a;
            v = -abc;
            a = support(shape2,shape1,v); %Tetrahedron new point
        end
    end

end

function point = getFarthestInDir(shape, v)
    %Find the furthest point in a given direction for a shape
    XData = shape.XData; % Making it more compatible with previous MATLAB releases.
    YData = shape.YData;
    ZData = shape.ZData;
    dotted = XData*v(1) + YData*v(2) + ZData*v(3);
    [maxInCol,rowIdxSet] = max(dotted);
    [maxInRow,colIdx] = max(maxInCol);
    rowIdx = rowIdxSet(colIdx);
    point = [XData(rowIdx,colIdx), YData(rowIdx,colIdx), ZData(rowIdx,colIdx)];
end

function point = support(shape1,shape2,v)
    %Support function to get the Minkowski difference.
    point1 = getFarthestInDir(shape1, v);
    point2 = getFarthestInDir(shape2, -v);
    point = point1 - point2;
end

function EQUIL=eqtr3m4Bars(NODE,ELEM)

    [nnode,var]=size(NODE);
    if (var ~= 6),
       error('ERROR in EQTR3: NODE matrix of incorrect size')
    end
    [nelem,~]=size(ELEM);

    icount=0;
    for i=1:nnode
        for j=1:3
            if (NODE(i,j+3) == 0),
               icount=icount+1;
               ROWNO(i,j)=icount;
            else
               ROWNO(i,j)=0;
            end
        end
    end

    ndof=icount;
    EQUIL=zeros(ndof,nelem);
    icol=0;
    for ielem=1:nelem
        icol=icol+1;
        C=[NODE(ELEM(ielem,1),1:3),NODE(ELEM(ielem,2),1:3)];
        len=sqrt((C(1)-C(4))^2+(C(2)-C(5))^2+(C(3)-C(6))^2);
        x=(C(1)-C(4))/len;
        y=(C(2)-C(5))/len;
        z=(C(3)-C(6))/len;
        EQ1=[x y z -x -y -z];
        ii=0;
        for i=1:2
            for j=1:3
                ii=ii+1;
                irow=ROWNO(ELEM(ielem,i),j);
                if (irow > 0)
                   EQUIL(irow,icol)=EQ1(ii);
                end
            end
        end
    end
end
