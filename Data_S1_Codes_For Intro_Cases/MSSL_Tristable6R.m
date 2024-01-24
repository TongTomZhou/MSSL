%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Theme:     This is a supplementary code of the paper: Multi-stable spatial linkages.
%%%%            Run it to generate a single-loop tri-stable 6R linkage.
%%%% Author:    Tong Zhou
%%%% Date:      2024.01.01
%%%% Result of the introductory case in Fig. 2:
%%%%            l12=1.708427154564838, l23=3.049554489692892, l34=0.940115800225138, l45=1.456822659113409
%%%%            l56=2.870978530732858, l61=4.211947674405008.
%%%%            alpha12=1.804085518504014, alpha23=1.995166185591464, alpha34=-0.967082660175891, alpha45=2.035929271737798,
%%%%            alpha56=2.113645252971641, alpha61=0.947551516569043.
%%%%            r1=1.015828321667094, r2=-3.388358671691819, r3=-4.074611415490486, r4=4.435799173784098,
%%%%            r5=5.960533532825726, r6=-0.392647294669189.
%%%%            d1=-3,      d2=3,        d3=-2,        d4=-0.5,    d5=4,   d6=-3.5.
%%%%            theta1S1=-0.523598775003225, theta2S1=2.035215449154171, theta3S1=-0.897580035529690, theta4S1=-2.733087481478449,
%%%%            theta5S1=-1.299145216602555, theta6S1=-1.134469559018412.
%%%%            theta1S2=1.203004090259966, theta2S2=1.553600631808526, theta3S2=-2.590320744210307, theta4S2=1.220253836917497,
%%%%            theta5S2=3.953604815742074, theta6S2=-0.234553079788255.
%%%%            theta1S3=2.093058411975421, theta2S3=2.087280217417274, theta3S3=0.226917871719025, theta4S3=5.298134876352753,
%%%%            theta5S3=2.586997478233854, theta6S3=-1.430268870561106.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
format short g

LISTr_ThetaTarget=[-30/180*pi,60/180*pi,120/180*pi];
INDX_FixedNodes=[1,6,7,12]; VALE_LimitDevia=1e-2; VALE_InteraMax=100000;

%% Initial Guess
MATX_ParaStrucInit=standardOCM(LISTr_ThetaTarget(1),LISTr_ThetaTarget(2),LISTr_ThetaTarget(3));

%% Generate 6R Tristable
A=[];b=[];  Aeq=[];beq=[];  lb=[];ub=[];
options=optimoptions('fmincon',... %  'ConstraintTolerance',1e-12,...
                     'MaxFunctionEvaluations',VALE_InteraMax,...
                     'MaxIterations',VALE_InteraMax);

LISTr_ParaStrucInit=reshape(MATX_ParaStrucInit,[],1);

obj6R=@(vars)objFunc6R(vars);
nonlcon6R=@(vars)conFunc6R(vars,LISTr_ThetaTarget,VALE_LimitDevia);

[LISTv_Vars6R,fval6R,exitflag6R]=...
    fmincon(obj6R,LISTr_ParaStrucInit,A,b,Aeq,beq,lb,ub,nonlcon6R,options);

MATX_ParaStruc=reshape(LISTv_Vars6R,6,6);

%% Post-Analysis
[PointsDOC1,PointsDPC1,PointsSOC1,PointsSPC1]=...
            form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                              zeros(1,6),MATX_ParaStruc(4,:));
[PointsDOC2,PointsDPC2,PointsSOC2,PointsSPC2]=...
            form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                              zeros(1,6),MATX_ParaStruc(5,:));
[PointsDOC3,PointsDPC3,PointsSOC3,PointsSPC3]=...
            form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                              zeros(1,6),MATX_ParaStruc(6,:));

MATX_PointsD(:,:,1)=[PointsDOC1;PointsDPC1];    MATX_PointsD(:,:,2)=[PointsDOC2;PointsDPC2];
MATX_PointsD(:,:,3)=[PointsDOC3;PointsDPC3];
MATX_PointsS(:,:,1)=[PointsSOC1;PointsSPC1];    MATX_PointsS(:,:,2)=[PointsSOC2;PointsSPC2];
MATX_PointsS(:,:,3)=[PointsSOC3;PointsSPC3];

VALE_StiffC1=calStiff6R(MATX_PointsS(:,:,1),INDX_FixedNodes,1);
VALE_StiffC2=calStiff6R(MATX_PointsS(:,:,2),INDX_FixedNodes,1);
VALE_StiffC3=calStiff6R(MATX_PointsS(:,:,3),INDX_FixedNodes,1);

VALE_DeviaTarget=abs(MATX_ParaStruc(4:6,1)-reshape(LISTr_ThetaTarget,[],1));

figure(1);clf;
subplot(3,3,1)
plot6RDH(MATX_PointsD(:,:,1),'','c');
title("Tristable 6R with DH notation at Stable state 1")
subplot(3,3,2)
plot6RDH(MATX_PointsD(:,:,2),'','m');
title("Tristable 6R with DH notation at Stable state 2")
subplot(3,3,3)
plot6RDH(MATX_PointsD(:,:,3),'','g');
title("Tristable 6R with DH notation at Stable state 3")
subplot(3,3,4)
plot6RStruc(MATX_PointsS(:,:,1),'','c');
title("Tristable 6R with form A at Stable state 1. Stiffness: "+string(VALE_StiffC1))
subplot(3,3,5)
plot6RStruc(MATX_PointsS(:,:,2),'','m');
title("Tristable 6R with form A at Stable state 2. Stiffness: "+string(VALE_StiffC2))
subplot(3,3,6)
plot6RStruc(MATX_PointsS(:,:,3),'','g');
title("Tristable 6R with form A at Stable state 3. Stiffness: "+string(VALE_StiffC3))


%% No-contact design
LISTr_SlipDist=[-3,3,-2,-0.5,4,-3.5]; VALE_LenHinge=[5,6,2,3,4,4.5];
MATX_PointsC_C1=1/2*(PointsDOC1+PointsDPC1);  MATX_PointsC_C2=1/2*(PointsDOC2+PointsDPC2);
MATX_PointsC_C3=1/2*(PointsDOC3+PointsDPC3);
MATX_VecHinge_C1=(PointsSPC1-PointsSOC1)/6; MATX_VecHinge_C2=(PointsSPC2-PointsSOC2)/6;
MATX_VecHinge_C3=(PointsSPC3-PointsSOC3)/6;

MATX_PointsC_C1=MATX_PointsC_C1+LISTr_SlipDist'.*MATX_VecHinge_C1;
MATX_PointsC_C2=MATX_PointsC_C2+LISTr_SlipDist'.*MATX_VecHinge_C2;
MATX_PointsC_C3=MATX_PointsC_C3+LISTr_SlipDist'.*MATX_VecHinge_C3;

MATX_PointsSP_C1=MATX_PointsC_C1+VALE_LenHinge'/2.*MATX_VecHinge_C1;
MATX_PointsSO_C1=MATX_PointsC_C1-VALE_LenHinge'/2.*MATX_VecHinge_C1;
MATX_PointsSP_C2=MATX_PointsC_C2+VALE_LenHinge'/2.*MATX_VecHinge_C2;
MATX_PointsSO_C2=MATX_PointsC_C2-VALE_LenHinge'/2.*MATX_VecHinge_C2;
MATX_PointsSP_C3=MATX_PointsC_C3+VALE_LenHinge'/2.*MATX_VecHinge_C3;
MATX_PointsSO_C3=MATX_PointsC_C3-VALE_LenHinge'/2.*MATX_VecHinge_C3;

MATX_PointsS(:,:,1)=[MATX_PointsSO_C1;MATX_PointsSP_C1];
MATX_PointsS(:,:,2)=[MATX_PointsSO_C2;MATX_PointsSP_C2];
MATX_PointsS(:,:,3)=[MATX_PointsSO_C3;MATX_PointsSP_C3];

VALE_StiffC1=calStiff6R(MATX_PointsS(:,:,1),INDX_FixedNodes,1);
VALE_StiffC2=calStiff6R(MATX_PointsS(:,:,2),INDX_FixedNodes,1);
VALE_StiffC3=calStiff6R(MATX_PointsS(:,:,3),INDX_FixedNodes,1);


figure(1)
subplot(3,3,7)
plot6RStruc(MATX_PointsS(:,:,1),'','c');
title("Tristable 6R with form B at Stable state 1. Stiffness: "+string(VALE_StiffC1))
subplot(3,3,8)
plot6RStruc(MATX_PointsS(:,:,2),'','m');
title("Tristable 6R with form B at Stable state 2. Stiffness: "+string(VALE_StiffC2))
subplot(3,3,9)
plot6RStruc(MATX_PointsS(:,:,3),'','g');
title("Tristable 6R with form B at Stable state 3. Stiffness: "+string(VALE_StiffC3))





%% Addtional Functions
function obj=objFunc6R(vars)
     FixedNodes=[1,6,7,12];
    MATX_ParaStruc=reshape(vars',6,6);

    [~,~,pointsSOC1,pointsSPC1]=form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                                                zeros(1,6),MATX_ParaStruc(4,:));
    [~,~,pointsSOC2,pointsSPC2]=form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                                                zeros(1,6),MATX_ParaStruc(5,:));
    [~,~,pointsSOC3,pointsSPC3]=form6RStruc2(MATX_ParaStruc(1,:),MATX_ParaStruc(2,:),MATX_ParaStruc(3,:),...
                                                zeros(1,6),MATX_ParaStruc(6,:));

    stiffC1Val=calStiff6R([pointsSOC1;pointsSPC1],FixedNodes,100);
    stiffC2Val=calStiff6R([pointsSOC2;pointsSPC2],FixedNodes,100);
    stiffC3Val=calStiff6R([pointsSOC3;pointsSPC3],FixedNodes,100);

    obj=-(real(stiffC1Val)+real(stiffC2Val)+real(stiffC3Val))/3;

end

function [c,ceq]=conFunc6R(vars,ThetaTarget,LimitDevia)
    c=[];   ceq=[];
    ParaStruc=reshape(vars',6,6);

    [cClosC1,ceqClosC1]=con6RClosure(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),ParaStruc(4,:));
    [cClosC2,ceqClosC2]=con6RClosure(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),ParaStruc(5,:));
    [cClosC3,ceqClosC3]=con6RClosure(ParaStruc(1,:),ParaStruc(2,:),ParaStruc(3,:),ParaStruc(6,:));

    [cShape,ceqShape]=con6RShape(ParaStruc,"a");
    [cTarget,ceqTarget]=con6RTarget(ParaStruc,ThetaTarget,LimitDevia);

    c=[c;cClosC1;cClosC2;cClosC3;cShape;cTarget];
    ceq=[ceq;ceqClosC1;ceqClosC2;ceqClosC3;ceqShape;ceqTarget];
end

function [c,ceq]=con6RTarget(ParaStruc,ThetaTarget,LimitDevia)
    c=[];   ceq=[];
    DeviaTarget1=abs(ParaStruc(4,1)-ThetaTarget(1));
    DeviaTarget2=abs(ParaStruc(5,1)-ThetaTarget(2));
    DeviaTarget3=abs(ParaStruc(6,1)-ThetaTarget(3));

    c=(DeviaTarget1+DeviaTarget2+DeviaTarget3)/3/(2*pi)-LimitDevia;
end

function [c,ceq]=con6RShape(MATX_ParaStruc,type)
    c=[];   ceq=[];

    lenList=MATX_ParaStruc(1,:);  alphaList=MATX_ParaStruc(2,:);    offsetList=MATX_ParaStruc(3,:);
    thetaListC1=MATX_ParaStruc(4,:);    thetaListC2=MATX_ParaStruc(5,:);    thetaListC3=MATX_ParaStruc(6,:);

    c=[c;-reshape(lenList,[],1)+1e-3;reshape(lenList,[],1)-20];
    c=[c;-reshape(alphaList,[],1)-pi;reshape(alphaList,[],1)-pi];
%         -abs(reshape(alphaList,[],1))+1e-3];
    c=[c;reshape(offsetList,[],1)-20];

%    c=[c;-reshape(thetaListC1,[],1)-2*pi;reshape(thetaListC1,[],1)-2*pi;...
%         -reshape(thetaListC2,[],1)-2*pi;reshape(thetaListC2,[],1)-2*pi;...
%         -reshape(thetaListC3,[],1)-2*pi;reshape(thetaListC3,[],1)-2*pi];

    if type=="a"
        c1=abs(reshape(lenList(1:3),[],1)-reshape(lenList(4:6),[],1));
        c2=abs(reshape(alphaList(1:3),[],1)-reshape(alphaList(4:6),[],1));
        c3=abs(reshape(offsetList(1:3),[],1)-reshape(offsetList(4:6),[],1));
        c=[c;-(sum(c1)+sum(c2)+sum(c3))+1e-3];
    elseif type=="b"
        c1=-abs(reshape(lenList(1:3),[],1)-reshape(lenList([6,5,4]),[],1));
        c2=-abs(reshape(alphaList(1:3),[],1)+reshape(alphaList(4:6),[],1)-2*pi*ones(3,1));
        c3=[-abs(offsetList(1));-abs(offsetList(4));...
            -abs(offsetList(2)+offsetList(6));-abs(offsetList(3)+offsetList(5))];
        c=[c;-(sum(c1)+sum(c2)+sum(c3))+1e-3];
    elseif type=="c"
        c1=-abs(lenList(1)^2+lenList(2)^2+lenList(5)^2-(lenList(2)^2+lenList(4)^2+lenList(5)^2));
        c2=-abs(reshape(alphaList,[],1)-[pi/2;3*pi/2;pi/2;3*pi/2;pi/2;3*pi/2]);
        c3=-abs(reshape(offsetList,[],1));
        c=[c;-(sum(c1)+sum(c2)+sum(c3))+1e-3];
    elseif type=="d"
        c1=-abs(reshape(lenList,[],1));
        c2=-abs(reshape(offsetList(1:3),[],1)+reshape(offsetList(4:6),[],1));
        c=[c;-(sum(c1)+sum(c2))+1e-3];
    elseif type=="e"
        c1=-abs(reshape(lenList,[],1));
        c2=[-abs(offsetList(1)+offsetList(4));...
            -abs(offsetList(2)+offsetList(1)*sin(alphaList(3))/sin(alphaList(1)+alphaList(2)));...
            -abs(offsetList(3)-offsetList(1)*sin(alphaList(1))/sin(alphaList(1)+alphaList(2)));...
            -abs(offsetList(5)-offsetList(1)*sin(alphaList(6))/sin(alphaList(4)+alphaList(6)));...
            -abs(offsetList(6)+offsetList(1)*sin(alphaList(4))/sin(alphaList(4)+alphaList(6)))];
        c=[c;-(sum(c1)+sum(c2))+1e-3];
    elseif type=="f"
        c1=-abs(reshape(lenList,[],1));
        c2=-abs(offsetList(1)*offsetList(3)*offsetList(5)+offsetList(2)*offsetList(4)*offsetList(6));
        c=[c;-(sum(c1)+sum(c2))+1e-3];
    end
end

function dhPara=standardOCM(theta1C1,theta1C2,theta1C3)
    MATX_ThetaTarget=ones(3,6);
    MATX_ThetaTarget=[theta1C1;theta1C2;theta1C3].*MATX_ThetaTarget;
    thetaListInit=reshape(MATX_ThetaTarget',[],1);

    alpha12=1/3*pi; alpha23=1/2*pi;   alpha34=2/3*pi;
    alpha45=alpha12;    alpha56=alpha23;    alpha61=alpha34;
    alphaListInit=[alpha12;alpha23;alpha34;alpha45;alpha56;alpha61];

    l12=4;   l23=6;   l34=8;   l45=l12;   l56=l23;   l61=l34;
    lenListInit=[l12;l23;l34;l45;l56;l61];

    offset1=1;      offset2=2;     offset3=3;
    offset4=offset1;    offset5=offset2;    offset6=offset3;
    offsetListInit=[offset1;offset2;offset3;offset4;offset5;offset6];

    varsInit=[lenListInit;alphaListInit;offsetListInit;thetaListInit];
    A=[];b=[];  Aeq=[];beq=[];  lb=[];ub=[];
    options=optimoptions('fmincon',... %  'ConstraintTolerance',1e-12,...
                         'MaxFunctionEvaluations',100000,...
                         'MaxIterations',100000);

    obj6R=@(vars)objFuncInit6R(vars,MATX_ThetaTarget);
    nonlcon6R=@(vars)conFuncInit6R(vars);

    [vars6R,fval6R,exitflag6R]=...
                    fmincon(obj6R,varsInit,A,b,Aeq,beq,lb,ub,nonlcon6R,options);
    dhPara=[vars6R(1:6),vars6R(7:12),vars6R(13:18),vars6R(19:24),vars6R(25:30),vars6R(31:36)]';
end

function obj=objFuncInit6R(vars,MATX_ThetaTarget)
    thetaList=[vars(19:24),vars(25:30),vars(31:36)]';
    obj=MATX_ThetaTarget(:,1)-thetaList(:,1);
    obj=abs(obj(1))+abs(obj(2))+abs(obj(3));
end

function [c,ceq]=conFuncInit6R(vars)
    c=[];   ceq=[];

    lList=vars(1:6);  alphaList=vars(7:12);   offsetList=vars(13:18);
    thetaList=[vars(19:24),vars(25:30),vars(31:36)]';

    [cClosC1,ceqClosC1]=con6RClosure(lList,alphaList,offsetList,thetaList(1,:));
    [cClosC2,ceqClosC2]=con6RClosure(lList,alphaList,offsetList,thetaList(2,:));
    [cClosC3,ceqClosC3]=con6RClosure(lList,alphaList,offsetList,thetaList(3,:));

    [cShape,ceqShape]=conFuncBricard(lList,alphaList,offsetList,"a");

    c=[c;cClosC1;cClosC2;cClosC3;cShape];
    ceq=[ceq;ceqClosC1;ceqClosC2;ceqClosC3;ceqShape];
end

function [c,ceq]=conFuncBricard(lenList,alphaList,offsetList,type)
    %   type: "a" the general line-symmetric case.     "b" the general plane-symmetric case.
    %         "c" the trihedral case.                  "d" The line-symmetric octahedral case.
    %         "e" The plane-symmetric octahedral case. "f" The doubly collapsible octahedral case.
    c=[];   ceq=[];
    if type=="a"
        ceq1=reshape(lenList(1:3),[],1)-reshape(lenList(4:6),[],1);
        ceq2=reshape(alphaList(1:3),[],1)-reshape(alphaList(4:6),[],1);
        ceq3=reshape(offsetList(1:3),[],1)-reshape(offsetList(4:6),[],1);
        ceq=[ceq1;ceq2;ceq3];
    elseif type=="b"
        ceq1=reshape(lenList(1:3),[],1)-reshape(lenList([6,5,4]),[],1);
        ceq2=reshape(alphaList(1:3),[],1)+reshape(alphaList(4:6),[],1)-2*pi*ones(3,1);
        ceq3=[offsetList(1);offsetList(4);offsetList(2)+offsetList(6);offsetList(3)+offsetList(5)];
        ceq=[ceq1;ceq2;ceq3];
    elseif type=="c"
        ceq1=lenList(1)^2+lenList(2)^2+lenList(5)^2-(lenList(2)^2+lenList(4)^2+lenList(5)^2);
        ceq2=reshape(alphaList,[],1)-[pi/2;3*pi/2;pi/2;3*pi/2;pi/2;3*pi/2];
        ceq3=reshape(offsetList,[],1);
        ceq=[ceq1;ceq2;ceq3];
    elseif type=="d"
        ceq1=reshape(lenList,[],1);
        ceq2=reshape(offsetList(1:3),[],1)+reshape(offsetList(4:6),[],1);
        ceq=[ceq1;ceq2];
    elseif type=="e"
        ceq1=reshape(lenList,[],1);
        ceq2=[offsetList(1)+offsetList(4);...
              offsetList(2)+offsetList(1)*sin(alphaList(3))/sin(alphaList(1)+alphaList(2));...
              offsetList(3)-offsetList(1)*sin(alphaList(1))/sin(alphaList(1)+alphaList(2));...
              offsetList(5)-offsetList(1)*sin(alphaList(6))/sin(alphaList(4)+alphaList(6));...
              offsetList(6)+offsetList(1)*sin(alphaList(4))/sin(alphaList(4)+alphaList(6))];
        ceq=[ceq1;ceq2];
    elseif type=="f"
        ceq1=reshape(lenList,[],1);
        ceq2=offsetList(1)*offsetList(3)*offsetList(5)+offsetList(2)*offsetList(4)*offsetList(6);
        ceq=[ceq1;ceq2];
    end

    % Basic Shape Constraints
    c=[c;-reshape(lenList,[],1)+1e-3;reshape(lenList,[],1)-20];
    c=[c;-reshape(alphaList,[],1)-pi-1e-3;reshape(alphaList,[],1)-pi+1e-3;...
         -abs(reshape(alphaList,[],1))+1e-3];
    c=[c;reshape(offsetList,[],1)-20];
end

function [c,ceq]=con6RClosure(lenList,alphaList,offsetList,thetaList)

    c=[];   ceq=[];

    % Closure Condition
    trans21=transMatDH(thetaList(1),alphaList(1),lenList(1),offsetList(1),'DES');
    trans32=transMatDH(thetaList(2),alphaList(2),lenList(2),offsetList(2),'DES');
    trans43=transMatDH(thetaList(3),alphaList(3),lenList(3),offsetList(3),'DES');

    trans45=transMatDH(thetaList(4),alphaList(4),lenList(4),offsetList(4),'ASC');
    trans56=transMatDH(thetaList(5),alphaList(5),lenList(5),offsetList(5),'ASC');
    trans61=transMatDH(thetaList(6),alphaList(6),lenList(6),offsetList(6),'ASC');

    conditionClosure6R=trans21*trans32*trans43-trans61*trans56*trans45;
    conditionClosure6R=reshape(conditionClosure6R,[],1);
    ceq=[ceq;conditionClosure6R];
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

function [pointsDO,pointsDP,pointsSO,pointsSP]=form6RStruc2(lVars,alphaVars,offsetVars,...
                                                            slipEdgeVars,thetaVars)

    l12=lVars(1);   l23=lVars(2);   l34=lVars(3);   l45=lVars(4);   l56=lVars(5);   l61=lVars(6);
    alpha12=alphaVars(1);   alpha23=alphaVars(2);   alpha34=alphaVars(3);
    alpha45=alphaVars(4);   alpha56=alphaVars(5);   alpha61=alphaVars(6);
    r1=offsetVars(1);   r2=offsetVars(2);   r3=offsetVars(3);   r4=offsetVars(4);
    r5=offsetVars(5);   r6=offsetVars(6);
    theta1=thetaVars(1);   theta2=thetaVars(2);   theta3=thetaVars(3);   theta4=thetaVars(4);
    theta5=thetaVars(5);   theta6=thetaVars(6);

    slipEdgeC1=slipEdgeVars(1);   slipEdgeC2=slipEdgeVars(2);   slipEdgeC3=slipEdgeVars(3);
    slipEdgeC4=slipEdgeVars(4);   slipEdgeC5=slipEdgeVars(5);   slipEdgeC6=slipEdgeVars(6);

    % P1
    O1=[0;0;0];
    sys1X=[1;0;0]; sys1Z=[0;0;1];
    p1=O1+r1*sys1Z;
    vec12=rotawithNormVec(sys1X,sys1Z,theta1);
    C1=1/2*(O1+p1);

    % P2
    O2=p1+l12*vec12;
    sys2X=vec12;   sys2Z=rotawithNormVec(sys1Z,sys2X,alpha12);
    p2=O2+r2*sys2Z;
    vec23=rotawithNormVec(sys2X,sys2Z,theta2);
    C2=1/2*(O2+p2);

    % P3
    O3=p2+l23*vec23;
    sys3X=vec23;   sys3Z=rotawithNormVec(sys2Z,sys3X,alpha23);
    p3=O3+r3*sys3Z;
    vec34=rotawithNormVec(sys3X,sys3Z,theta3);
    C3=1/2*(O3+p3);

    % P4
    O4=p3+l34*vec34;
    sys4X=vec34;   sys4Z=rotawithNormVec(sys3Z,sys4X,alpha34);
    p4=O4+r4*sys4Z;
    vec45=rotawithNormVec(sys4X,sys4Z,theta4);
    C4=1/2*(O4+p4);

    % P5
    O5=p4+l45*vec45;
    sys5X=vec45;   sys5Z=rotawithNormVec(sys4Z,sys5X,alpha45);
    p5=O5+r5*sys5Z;
    vec56=rotawithNormVec(sys5X,sys5Z,theta5);
    C5=1/2*(O5+p5);

    % P6
    O6=p5+l56*vec56;
    sys6X=vec56;   sys6Z=rotawithNormVec(sys5Z,sys6X,alpha56);
    p6=O6+r6*sys6Z;
    vec61=rotawithNormVec(sys6X,sys6Z,theta6);
    C6=1/2*(O6+p6);

%    lenHinge=mean(lVars);
    lenHinge=6; % Modified at 2023.04.02.
    O1H=C1-lenHinge/2*sys1Z;    p1H=C1+lenHinge/2*sys1Z;
    O2H=C2-lenHinge/2*sys2Z;    p2H=C2+lenHinge/2*sys2Z;
    O3H=C3-lenHinge/2*sys3Z;    p3H=C3+lenHinge/2*sys3Z;
    O4H=C4-lenHinge/2*sys4Z;    p4H=C4+lenHinge/2*sys4Z;
    O5H=C5-lenHinge/2*sys5Z;    p5H=C5+lenHinge/2*sys5Z;
    O6H=C6-lenHinge/2*sys6Z;    p6H=C6+lenHinge/2*sys6Z;

    pointsDO=[O1';O2';O3';O4';O5';O6'];
    pointsDP=[p1';p2';p3';p4';p5';p6'];
    pointsSO=[O1H';O2H';O3H';O4H';O5H';O6H'];
    pointsSP=[p1H';p2H';p3H';p4H';p5H';p6H'];
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

function eValStiff=calStiff6R(points,indexFixnodes,K)
    topoLink=[];topoBricardTetra=[1,7,2,8;2,8,3,9;3,9,4,10;4,10,5,11;5,11,6,12;6,12,1,7];
    for i=1:size(topoBricardTetra,1)
        topoLink=[topoLink;topoBricardTetra(i,1),topoBricardTetra(i,3);...
                           topoBricardTetra(i,1),topoBricardTetra(i,4);...
                           topoBricardTetra(i,2),topoBricardTetra(i,3);...
                           topoBricardTetra(i,2),topoBricardTetra(i,4)];
    end
    topoLink=[topoLink;topoBricardTetra(:,1:2)];

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

function plot6RDH(points,figureTitle,color)
    topoBars=[1,7;7,2;2,8;8,3;3,9;9,4;4,10;10,5;5,11;11,6;6,12;12,1];
    offsetDist=norm(points(1,:)-points(2,:))/10;

    for i=1:size(points,1)/2
        plot3(points(i,1),points(i,2),points(i,3),'ko')
        hold on
        text(points(i,1),points(i,2)+offsetDist,points(i,3)+offsetDist,string(i))
        hold on
    end
    for i=size(points,1)/2+1:size(points,1)
        plot3(points(i,1),points(i,2),points(i,3),'k+')
        hold on
    end

    for i=1:2:size(topoBars,1)
        plot3([points(topoBars(i,1),1),points(topoBars(i,2),1)],...
              [points(topoBars(i,1),2),points(topoBars(i,2),2)],...
              [points(topoBars(i,1),3),points(topoBars(i,2),3)],'k--')
        hold on
    end

    for i=2:2:size(topoBars,1)
        plot3([points(topoBars(i,1),1),points(topoBars(i,2),1)],...
              [points(topoBars(i,1),2),points(topoBars(i,2),2)],...
              [points(topoBars(i,1),3),points(topoBars(i,2),3)],color)
        hold on
    end

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    title(figureTitle)
end

function plot6RStruc(points,figureTitle,color)
    topoBarStrucTetra=[1,7,2,8;2,8,3,9;3,9,4,10;4,10,5,11;5,11,6,12;6,12,1,7];

    tetramesh(topoBarStrucTetra,points,'FaceColor',color,'FaceAlpha',0.2,'EdgeColor','k')
    hold on

    for i=1:size(points,1)/2
        plot3(points(i,1),points(i,2),points(i,3),'ko')
        hold on
        text(points(i,1),points(i,2),points(i,3)+0.05,'O'+string(i));
        hold on
    end

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    title(figureTitle)
end