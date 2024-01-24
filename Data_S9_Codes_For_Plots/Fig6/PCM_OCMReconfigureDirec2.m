clear; clc;%close all;
tool=TOOLS_SharedFunction();

%%%% General
theta1Init=0.567591958728399; theta2Init=2.954971782738739;
theta1End=2.617993877991494; theta2End=1.391583917266936;
LISTr_Theta1Range=theta1Init:(1/20*pi):theta1End;

MATX_Info=[];
for alpha=(-pi-0.3001):0.3:(pi+0.3001)
    %alpha=1.2441684642559037; beta=pi-alpha;
    constTheta1Theta2=tan(theta1Init/2)*tan(theta2Init/2);
    beta=2*acos(((-1+constTheta1Theta2)*cos(alpha/2))/(sqrt(1+constTheta1Theta2^2-2*constTheta1Theta2*cos(alpha))));
    constAlphaBeta=sin(beta/2+alpha/2)/sin(beta/2-alpha/2);

    TSOR_Points_Motion=Simu_BennettLinkage(alpha,beta,LISTr_Theta1Range);
    for i=1:size(TSOR_Points_Motion,3)
        theta2=2*atan(constAlphaBeta/tan(LISTr_Theta1Range(i)/2));

        [rho1_i,rho2_i]=Cal_ReconfigureDirec(i,TSOR_Points_Motion);
        MATX_Info=[MATX_Info;alpha,beta,LISTr_Theta1Range(i),theta2,rho1_i,rho2_i];
    end
end

%writematrix(MATX_Info,"MATX_PCM_OCMReconfigureDirec.csv");
papercolormap=tool.customcolormap(linspace(0,1,7), {'#AA1E29','#DC4638','#FDB26E','#E7F1D6','#7AB0D6','#3D4DA2','#3E51A1'});
[rho1_1,rho2_1,VecMotHinge1_1,VecMotHinge2_1]=Cal_ReconfigureDirec(10,TSOR_Points_Motion);
[rho1_2,rho2_2,VecMotHinge1_2,VecMotHinge2_2]=Cal_ReconfigureDirec(11,TSOR_Points_Motion);


MATX_PointsMot_OCM=TSOR_Points_Motion(:,:,10);
MATX_PointsMot_OCM_=TSOR_Points_Motion(:,:,11);
LISTv_Origin_OCM=transpose(1/4*sum(MATX_PointsMot_OCM([1,5,3,7],:)));
LISTv_Xdirc_OCM=transpose(1/2*sum(MATX_PointsMot_OCM([1,5],:))-LISTv_Origin_OCM');
LISTv_Ydirc_OCM=transpose(1/2*sum(MATX_PointsMot_OCM([2,6],:))-1/2*sum(MATX_PointsMot_OCM([4,8],:)));
LISTv_Zdirc_OCM=cross(LISTv_Xdirc_OCM,LISTv_Ydirc_OCM);
LISTv_Xdirc_OCM=LISTv_Xdirc_OCM/norm(LISTv_Xdirc_OCM);
LISTv_Ydirc_OCM=LISTv_Ydirc_OCM/norm(LISTv_Ydirc_OCM);
LISTv_Zdirc_OCM=LISTv_Zdirc_OCM/norm(LISTv_Zdirc_OCM);

LISTv_Origin_OCM_=transpose(1/4*sum(MATX_PointsMot_OCM_([1,5,3,7],:)));
LISTv_Xdirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([1,5],:))-LISTv_Origin_OCM_');
LISTv_Ydirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([2,6],:))-1/2*sum(MATX_PointsMot_OCM_([4,8],:)));
LISTv_Zdirc_OCM_=cross(LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_);
LISTv_Xdirc_OCM_=LISTv_Xdirc_OCM_/norm(LISTv_Xdirc_OCM_);
LISTv_Ydirc_OCM_=LISTv_Ydirc_OCM_/norm(LISTv_Ydirc_OCM_);
LISTv_Zdirc_OCM_=LISTv_Zdirc_OCM_/norm(LISTv_Zdirc_OCM_);
        

for ii=1:size(MATX_PointsMot_OCM,1)
    MATX_PointsMot_OCM(ii,:)=transpose(inv([LISTv_Xdirc_OCM,LISTv_Ydirc_OCM,LISTv_Zdirc_OCM])*(MATX_PointsMot_OCM(ii,:)'-...
                                           LISTv_Origin_OCM));
    MATX_PointsMot_OCM_(ii,:)=transpose(inv([LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_,LISTv_Zdirc_OCM_])*(MATX_PointsMot_OCM_(ii,:)'-...
                                           LISTv_Origin_OCM_));
end
        


figure(1);clf
tool.plot4RStruc(MATX_PointsMot_OCM,'','c');
tool.plot4RStruc(MATX_PointsMot_OCM_,'','m');
plot3([MATX_PointsMot_OCM(1,1),MATX_PointsMot_OCM(1,1)+VecMotHinge1_1(1)],...
      [MATX_PointsMot_OCM(1,2),MATX_PointsMot_OCM(1,2)+VecMotHinge1_1(2)],...
      [MATX_PointsMot_OCM(1,3),MATX_PointsMot_OCM(1,3)+VecMotHinge1_1(3)],'k--'); hold on;
plot3([MATX_PointsMot_OCM(2,1),MATX_PointsMot_OCM(2,1)+VecMotHinge2_1(1)],...
      [MATX_PointsMot_OCM(2,2),MATX_PointsMot_OCM(2,2)+VecMotHinge2_1(2)],...
      [MATX_PointsMot_OCM(2,3),MATX_PointsMot_OCM(2,3)+VecMotHinge2_1(3)],'k--'); hold on;



figure(3);clf
subplot(1,2,1)
plotRho1=scatter(MATX_Info(:,3),MATX_Info(:,1),200,MATX_Info(:,5),'filled');
%plotRho1.MarkerEdgeColor='k';
plotRho1.LineWidth=0.25;
colormap(papercolormap);
cb1=colorbar;
caxis([0, pi]);
%cb1.Ticks = [20,60,100,140,180]/180*pi;
%axis equal
xlim([30,150]/180*pi);   ylim([-pi-0.005,pi]);

figure(3);subplot(1,2,2)
plotRho2=scatter(MATX_Info(:,3),MATX_Info(:,1),200,MATX_Info(:,6),'filled');
%plotRho2.MarkerEdgeColor='k';
plotRho2.LineWidth=0.25;
colormap(papercolormap);
cb2=colorbar;
caxis([0, pi]);
%cb2.Ticks = [20,60,100,140,180]/180*pi;
%axis equal
xlim([30,150]/180*pi);   ylim([-pi-0.005,pi]);


%%%% Specific-Case 1
theta1_init_OCM=150/180*pi+(4.227727051568600-6.278128970831696);
theta1_end_OCM=150/180*pi;
LISTr_Theta1Range=theta1_init_OCM:(1/180*pi):theta1_end_OCM;

MATX_Info_=[];
alpha=1.2441684642559037;   beta=pi-alpha;

TSOR_Points_Motion=Simu_BennettLinkage(alpha,beta,LISTr_Theta1Range);
for ii=1:size(TSOR_Points_Motion,3)
    [rho1_i,rho2_i]=Cal_ReconfigureDirec(ii,TSOR_Points_Motion);
    MATX_Info_=[MATX_Info_;LISTr_Theta1Range(ii),rho1_i,rho2_i];
end

figure(2);clf;
plot(MATX_Info_(:,1),MATX_Info_(:,2),'c-'); hold on
plot(MATX_Info_(:,1),MATX_Info_(:,3),'m-'); hold on
plot(MATX_Info_(:,1),pi/2*ones(size(MATX_Info_(:,1))),'k--'); hold on
%xlim([30,150]/180*pi);
%ylim([15,105]/180*pi);
legend

%writematrix(MATX_Info_,"MATX_PCM_OCMReconfigureDirec_Spec.csv");


%figure(2)
%for i=1:size(TSOR_Points_Motion,3)
%    tool.plot4RStruc(TSOR_Points_Motion(:,:,i),'','m')
%    pause(0.5)
%    hold off;
%end

%calculate_angle_from_vector_to_vector(LISTv_VecHinge1_Right,LISTv_VecMotHinge1)

%figure(3)
%tool.plot4RStruc(MATX_Points_Right,'','m')
%plot3([LISTv_MidPotHinge1_Right(1),LISTv_MidPotHinge1_Right(1)+LISTv_VecMotHinge1(1)],...
%        [LISTv_MidPotHinge1_Right(2),LISTv_MidPotHinge1_Right(2)+LISTv_VecMotHinge1(2)],...
%        [LISTv_MidPotHinge1_Right(3),LISTv_MidPotHinge1_Right(3)+LISTv_VecMotHinge1(3)],'b-');hold on
%plot3([LISTv_MidPotHinge2_Right(1),LISTv_MidPotHinge2_Right(1)+LISTv_VecMotHinge2(1)],...
%        [LISTv_MidPotHinge2_Right(2),LISTv_MidPotHinge2_Right(2)+LISTv_VecMotHinge2(2)],...
%        [LISTv_MidPotHinge2_Right(3),LISTv_MidPotHinge2_Right(3)+LISTv_VecMotHinge2(3)],'b-');hold on


function TSOR_Points_Motion=Simu_BennettLinkage(alpha,beta,theta1Range)
    tool=TOOLS_SharedFunction();

    for i=1:length(theta1Range)
        theta1=theta1Range(i);
        theta2=2*atan(sin(beta/2+alpha/2)/sin(beta/2-alpha/2)/tan(theta1/2));

        MATX_ParaBenn=[4,4,4,4;...
                       alpha,beta,alpha,beta;...
                       zeros(1,4);...
                       theta1,theta2,-theta1,-theta2;...
                       1*ones(1,4);-1*ones(1,4)];

        [~,~,MATX_PointsSO_Benn,MATX_PointsSP_Benn]=...
                      tool.form4RStruc(MATX_ParaBenn(1,:),MATX_ParaBenn(2,:),MATX_ParaBenn(3,:),...
                                       MATX_ParaBenn(5,:),MATX_ParaBenn(6,:),MATX_ParaBenn(4,:));

        LISTv_SlipDist_Hinge1=3;
        LISTv_Hinge1=MATX_PointsSP_Benn(1,:)-MATX_PointsSO_Benn(1,:);
        LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
        MATX_PointsSO_Benn(1,:)=MATX_PointsSO_Benn(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
        MATX_PointsSP_Benn(1,:)=MATX_PointsSP_Benn(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;

        LISTv_SlipDist_Hinge3=-3;
        LISTv_Hinge3=MATX_PointsSP_Benn(3,:)-MATX_PointsSO_Benn(3,:);
        LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
        MATX_PointsSO_Benn(3,:)=MATX_PointsSO_Benn(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
        MATX_PointsSP_Benn(3,:)=MATX_PointsSP_Benn(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

        TSOR_Points_Motion(:,:,i)=[MATX_PointsSO_Benn;MATX_PointsSP_Benn];
    end
end

function [rho1_,rho2_,LISTv_VecMotHinge1,LISTv_VecMotHinge2]=... 
                    Cal_ReconfigureDirec(INDX_Points_Motion,TSOR_Points_Motion)
    VALE_ScaleRatio_OCM=10;  VALE_ScaleRatio_MSS=4;
    for i=1:size(TSOR_Points_Motion,3)
        MATX_PointsMot_OCM_=TSOR_Points_Motion(:,:,i)*VALE_ScaleRatio_OCM;

        LISTv_Origin_OCM_=transpose(1/4*sum(MATX_PointsMot_OCM_([1,5,3,7],:)));
        LISTv_Xdirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([1,5],:))-LISTv_Origin_OCM_');
        LISTv_Ydirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([2,6],:))-1/2*sum(MATX_PointsMot_OCM_([4,8],:)));
        LISTv_Zdirc_OCM_=cross(LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_);
        LISTv_Xdirc_OCM_=LISTv_Xdirc_OCM_/norm(LISTv_Xdirc_OCM_);
        LISTv_Ydirc_OCM_=LISTv_Ydirc_OCM_/norm(LISTv_Ydirc_OCM_);
        LISTv_Zdirc_OCM_=LISTv_Zdirc_OCM_/norm(LISTv_Zdirc_OCM_);

        for ii=1:size(MATX_PointsMot_OCM_,1)
            MATX_PointsMot_OCM(ii,:)=transpose(inv([LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_,LISTv_Zdirc_OCM_])*(MATX_PointsMot_OCM_(ii,:)'-...
                                                   LISTv_Origin_OCM_));
        end
        TSOR_Points_Motion(:,:,i)=MATX_PointsMot_OCM;
    end

    if INDX_Points_Motion<size(TSOR_Points_Motion,3)
        MATX_Points=TSOR_Points_Motion(:,:,INDX_Points_Motion);
        MATX_PointsR=TSOR_Points_Motion(:,:,INDX_Points_Motion+1);

        LISTv_VecHinge1=MATX_Points(1+4,:)'-MATX_Points(1,:)';
        LISTv_VecHinge1=LISTv_VecHinge1/norm(LISTv_VecHinge1);
        LISTv_VecHinge2=MATX_Points(2+4,:)'-MATX_Points(2,:)';
        LISTv_VecHinge2=LISTv_VecHinge2/norm(LISTv_VecHinge2);

        LISTv_MidPotHinge1=MATX_Points(1,:)';   LISTv_MidPotHinge1R=MATX_PointsR(1,:)';
        LISTv_MidPotHinge2=MATX_Points(2,:)';   LISTv_MidPotHinge2R=MATX_PointsR(2,:)';

        LISTv_VecMotHinge1=LISTv_MidPotHinge1R-LISTv_MidPotHinge1;
        LISTv_VecMotHinge1=LISTv_VecMotHinge1/norm(LISTv_VecMotHinge1);
        LISTv_VecMotHinge2=LISTv_MidPotHinge2R-LISTv_MidPotHinge2;
        LISTv_VecMotHinge2=LISTv_VecMotHinge2/norm(LISTv_VecMotHinge2);

        rho1=calculate_angle_from_vector_to_vector(LISTv_VecHinge1,LISTv_VecMotHinge1);
        rho2=calculate_angle_from_vector_to_vector(LISTv_VecHinge2,LISTv_VecMotHinge2);
        
        rho1_=calculate_angle_from_vector_to_vector([0;0;1],LISTv_VecMotHinge1);
        rho2_=calculate_angle_from_vector_to_vector([0;0;1],LISTv_VecMotHinge2);
        
        if rho1<=0
            rho1=-rho1;
        end
        if rho2<=0
            rho2=-rho2;
        end

    else
        rho1=0/0; rho2=0/0;
        rho1_=0/0; rho2_=0/0;
    end
end


function angle_rad = calculate_angle_from_vector_to_vector(start_vector, end_vector)
    dot_product = dot(start_vector, end_vector);
    magnitude_start = norm(start_vector);
    magnitude_end = norm(end_vector);

    cosine_angle = dot_product / (magnitude_start * magnitude_end);
    angle_rad = acos(cosine_angle);

    cross_product = cross(start_vector, end_vector);
    if dot(cross_product, [0; 0; 1]) < 0
        angle_rad = -angle_rad;
    end

%    angle_deg = rad2deg(angle_rad);
end


