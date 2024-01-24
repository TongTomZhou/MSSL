clear;clc;%close all
tool=TOOLS_SharedFunction;

%% Input the MSS's motion data
MATX_PointsMot_MSS=readmatrix('MATX_Points_Motion_Bistable4R.csv');
j=1;
for i=1:8:size(MATX_PointsMot_MSS,1)
    TSOR_PointsMot_MSS(:,:,j)=MATX_PointsMot_MSS(i:i+7,:);
    j=j+1;
end
TSOR_PointsMot_MSS(:,:,103)=[]; TSOR_PointsMot_MSS(:,:,102)=[];


for i=1:size(TSOR_PointsMot_MSS,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_PointsMot_MSS(:,:,i));
    theta1=MATX_DHStruc(4,1); theta2=MATX_DHStruc(4,2);theta3=MATX_DHStruc(4,3); theta4=MATX_DHStruc(4,4);
    if theta1<0
        theta1=theta1+2*pi;
    elseif theta1>pi
        theta1=theta1-2*pi;
    end
    if theta2<0
        theta2=theta2+2*pi;
    elseif theta2>pi
        theta2=theta2-2*pi;
    end
    if theta3<-0
        theta3=theta3+2*pi;
    elseif theta3>pi
        theta3=theta3-2*pi;
    end
    if theta4<0
        theta4=theta4+2*pi;
    elseif theta4>pi
        theta4=theta4-2*pi;
    end

    LISTr_Theta1(i)=theta1; LISTr_Theta2(i)=theta2; LISTr_Theta3(i)=theta3; LISTr_Theta4(i)=theta4;
end

figure(1);clf
plot(LISTr_Theta1);hold on
plot(LISTr_Theta2);hold on
plot(LISTr_Theta3);hold on
plot(LISTr_Theta4);hold on
legend()


%% Design the OCM parameter to be coupled
INDX_init_MSS=size(TSOR_PointsMot_MSS,3);   INDX_end_MSS=1;
theta1_init_MSS=LISTr_Theta3(INDX_init_MSS);
theta2_init_MSS=LISTr_Theta4(INDX_init_MSS);
theta1_end_MSS=LISTr_Theta3(INDX_end_MSS);
theta2_end_MSS=LISTr_Theta4(INDX_end_MSS);

alpha_OCM=1.2441684642559037; beta_OCM=pi-alpha_OCM;
theta1_init_OCM=150/180*pi-(theta1_end_MSS-theta1_init_MSS);
theta1_end_OCM=150/180*pi;

theta2_init_OCM=2*atan(1/cos(alpha_OCM)/tan(theta1_init_OCM/2));
theta2_end_OCM=2*atan(1/cos(alpha_OCM)/tan(theta1_end_OCM/2));

LISTr_Theta1_OCM=linspace(theta1_init_OCM,theta1_end_OCM,size(TSOR_PointsMot_MSS,3));
LISTr_Theta2_OCM=linspace(theta2_init_OCM,theta2_end_OCM,size(TSOR_PointsMot_MSS,3));

TSOR_PointsMot_OCM=Simu_BennettLinkage_(alpha_OCM,beta_OCM,LISTr_Theta1_OCM);

% Check the design
theta1_init_OCM-theta1_end_OCM;
theta1_init_MSS-theta1_end_MSS;

theta2_init_OCM-theta2_end_OCM;
theta2_init_MSS-theta2_end_MSS;

%figure(2);clf
for INDX_MSS=INDX_init_MSS:-1:INDX_end_MSS
    INDX_OCM=size(TSOR_PointsMot_MSS,3)-INDX_MSS+1;

    VALE_ScaleRatio_OCM=10;  VALE_ScaleRatio_MSS=4;
    MATX_PointsMot_OCM_=Simu_BennettLinkage(alpha_OCM,beta_OCM,LISTr_Theta1_OCM(INDX_OCM))*VALE_ScaleRatio_OCM;
    MATX_PointsMot_MSS=TSOR_PointsMot_MSS(:,:,INDX_MSS)*VALE_ScaleRatio_MSS;

    % Transform the OCM
    LISTv_Origin_OCM_=transpose(1/4*sum(MATX_PointsMot_OCM_([1,5,3,7],:)));
    LISTv_Xdirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([1,5],:))-LISTv_Origin_OCM_');
    LISTv_Ydirc_OCM_=transpose(1/2*sum(MATX_PointsMot_OCM_([2,6],:))-1/2*sum(MATX_PointsMot_OCM_([4,8],:)));
    LISTv_Zdirc_OCM_=cross(LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_);
    LISTv_Xdirc_OCM_=LISTv_Xdirc_OCM_/norm(LISTv_Xdirc_OCM_);
    LISTv_Ydirc_OCM_=LISTv_Ydirc_OCM_/norm(LISTv_Ydirc_OCM_);
    LISTv_Zdirc_OCM_=LISTv_Zdirc_OCM_/norm(LISTv_Zdirc_OCM_);

    for i=1:size(MATX_PointsMot_MSS,1)
        MATX_PointsMot_OCM(i,:)=transpose(inv([LISTv_Xdirc_OCM_,LISTv_Ydirc_OCM_,LISTv_Zdirc_OCM_])*(MATX_PointsMot_OCM_(i,:)'-...
                                               LISTv_Origin_OCM_));
    end
    TSOR_PointsMot_OCM_(:,:,INDX_OCM)=MATX_PointsMot_OCM;
    
    % Get the target orientation of OCM
    % for MSS1
    LISTv_Origin_OCM=MATX_PointsMot_OCM(1,:)';
    LISTv_Zdirc_OCM=MATX_PointsMot_OCM(1+4,:)'-MATX_PointsMot_OCM(1,:)';

    LISTv_Xdirc1_OCM=(MATX_PointsMot_OCM(2,:)'-MATX_PointsMot_OCM(1,:)'+...
                      (-(MATX_PointsMot_OCM(5,:)-MATX_PointsMot_OCM(1,:))*(MATX_PointsMot_OCM(2,:)'-MATX_PointsMot_OCM(1,:)')/...
                        (norm(MATX_PointsMot_OCM(5,:)-MATX_PointsMot_OCM(1,:))^2))*...
                        (MATX_PointsMot_OCM(5,:)'-MATX_PointsMot_OCM(1,:)'));
    LISTv_Xdirc1_OCM=LISTv_Xdirc1_OCM/norm(LISTv_Xdirc1_OCM);
    LISTv_Xdirc2_OCM=(MATX_PointsMot_OCM(8,:)'-MATX_PointsMot_OCM(1,:)'+...
                      (-(MATX_PointsMot_OCM(5,:)-MATX_PointsMot_OCM(1,:))*(MATX_PointsMot_OCM(8,:)'-MATX_PointsMot_OCM(1,:)')/...
                        (norm(MATX_PointsMot_OCM(5,:)-MATX_PointsMot_OCM(1,:))^2))*...
                        (MATX_PointsMot_OCM(5,:)'-MATX_PointsMot_OCM(1,:)'));
    LISTv_Xdirc2_OCM=LISTv_Xdirc2_OCM/norm(LISTv_Xdirc2_OCM);

    LISTv_Xdirc_OCM=LISTv_Xdirc1_OCM+LISTv_Xdirc2_OCM;
    if INDX_MSS>90
        LISTv_Xdirc_OCM=-LISTv_Xdirc_OCM/norm(LISTv_Xdirc_OCM);
    else
        LISTv_Xdirc_OCM=LISTv_Xdirc_OCM/norm(LISTv_Xdirc_OCM);
    end

    LISTv_Ydirc_OCM=cross(LISTv_Zdirc_OCM,LISTv_Xdirc_OCM);
    LISTv_Ydirc_OCM=LISTv_Ydirc_OCM/norm(LISTv_Ydirc_OCM);

    % for MSS2
    LISTv_Origin_OCM2=MATX_PointsMot_OCM(2,:)';
    LISTv_Zdirc_OCM2=MATX_PointsMot_OCM(2+4,:)'-MATX_PointsMot_OCM(2,:)';

    LISTv_Xdirc1_OCM2=(MATX_PointsMot_OCM(5,:)'-MATX_PointsMot_OCM(2,:)'+...
                      (-(MATX_PointsMot_OCM(6,:)-MATX_PointsMot_OCM(2,:))*(MATX_PointsMot_OCM(5,:)'-MATX_PointsMot_OCM(2,:)')/...
                        (norm(MATX_PointsMot_OCM(6,:)-MATX_PointsMot_OCM(2,:))^2))*...
                        (MATX_PointsMot_OCM(6,:)'-MATX_PointsMot_OCM(2,:)'));
    LISTv_Xdirc1_OCM2=LISTv_Xdirc1_OCM2/norm(LISTv_Xdirc1_OCM2);
    LISTv_Xdirc2_OCM2=(MATX_PointsMot_OCM(3,:)'-MATX_PointsMot_OCM(2,:)'+...
                      (-(MATX_PointsMot_OCM(6,:)-MATX_PointsMot_OCM(2,:))*(MATX_PointsMot_OCM(3,:)'-MATX_PointsMot_OCM(2,:)')/...
                        (norm(MATX_PointsMot_OCM(6,:)-MATX_PointsMot_OCM(2,:))^2))*...
                        (MATX_PointsMot_OCM(6,:)'-MATX_PointsMot_OCM(2,:)'));
    LISTv_Xdirc2_OCM2=LISTv_Xdirc2_OCM2/norm(LISTv_Xdirc2_OCM2);

    LISTv_Xdirc_OCM2=LISTv_Xdirc1_OCM2+LISTv_Xdirc2_OCM2;
    if INDX_MSS>7
        LISTv_Xdirc_OCM2=-LISTv_Xdirc_OCM2/norm(LISTv_Xdirc_OCM2);
    else
        LISTv_Xdirc_OCM2=LISTv_Xdirc_OCM2/norm(LISTv_Xdirc_OCM2);
    end

    LISTv_Ydirc_OCM2=cross(LISTv_Zdirc_OCM2,LISTv_Xdirc_OCM2);
    LISTv_Ydirc_OCM2=LISTv_Ydirc_OCM2/norm(LISTv_Ydirc_OCM2);

    % Transform the MSS1
    LISTv_Origin_MSS=MATX_PointsMot_MSS(3,:)';
    LISTv_Zdirc_MSS=MATX_PointsMot_MSS(3+4,:)'-MATX_PointsMot_MSS(3,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_PointsMot_MSS(6,:)'-MATX_PointsMot_MSS(3,:)'+...
                      (-(MATX_PointsMot_MSS(7,:)-MATX_PointsMot_MSS(3,:))*(MATX_PointsMot_MSS(6,:)'-MATX_PointsMot_MSS(3,:)')/...
                        (norm(MATX_PointsMot_MSS(7,:)-MATX_PointsMot_MSS(3,:))^2))*...
                        (MATX_PointsMot_MSS(7,:)'-MATX_PointsMot_MSS(3,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_PointsMot_MSS(4,:)'-MATX_PointsMot_MSS(3,:)'+...
                      (-(MATX_PointsMot_MSS(7,:)-MATX_PointsMot_MSS(3,:))*(MATX_PointsMot_MSS(4,:)'-MATX_PointsMot_MSS(3,:)')/...
                        (norm(MATX_PointsMot_MSS(7,:)-MATX_PointsMot_MSS(3,:))^2))*...
                        (MATX_PointsMot_MSS(7,:)'-MATX_PointsMot_MSS(3,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    % Transform the MSS2
    LISTv_Origin_MSS2=MATX_PointsMot_MSS(8,:)';
    LISTv_Zdirc_MSS2=MATX_PointsMot_MSS(4+4,:)'-MATX_PointsMot_MSS(4,:)';
    LISTv_Zdirc_MSS2=LISTv_Zdirc_MSS2/norm(LISTv_Zdirc_MSS2);

    LISTv_Xdirc1_MSS2=(MATX_PointsMot_MSS(1,:)'-MATX_PointsMot_MSS(8,:)'+...
                      (-(MATX_PointsMot_MSS(4,:)-MATX_PointsMot_MSS(8,:))*(MATX_PointsMot_MSS(1,:)'-MATX_PointsMot_MSS(8,:)')/...
                        (norm(MATX_PointsMot_MSS(4,:)-MATX_PointsMot_MSS(8,:))^2))*...
                        (MATX_PointsMot_MSS(4,:)'-MATX_PointsMot_MSS(8,:)'));
    LISTv_Xdirc1_MSS2=LISTv_Xdirc1_MSS2/norm(LISTv_Xdirc1_MSS2);
    LISTv_Xdirc2_MSS2=(MATX_PointsMot_MSS(7,:)'-MATX_PointsMot_MSS(8,:)'+...
                      (-(MATX_PointsMot_MSS(4,:)-MATX_PointsMot_MSS(8,:))*(MATX_PointsMot_MSS(7,:)'-MATX_PointsMot_MSS(8,:)')/...
                        (norm(MATX_PointsMot_MSS(4,:)-MATX_PointsMot_MSS(8,:))^2))*...
                        (MATX_PointsMot_MSS(4,:)'-MATX_PointsMot_MSS(8,:)'));
    LISTv_Xdirc2_MSS2=LISTv_Xdirc2_MSS2/norm(LISTv_Xdirc2_MSS2);

    LISTv_Xdirc_MSS2=LISTv_Xdirc1_MSS2+LISTv_Xdirc2_MSS2;
    LISTv_Xdirc_MSS2=LISTv_Xdirc_MSS2/norm(LISTv_Xdirc_MSS2);
    LISTv_Ydirc_MSS2=cross(LISTv_Zdirc_MSS2,LISTv_Xdirc_MSS2);
    LISTv_Ydirc_MSS2=LISTv_Ydirc_MSS2/norm(LISTv_Ydirc_MSS2);

    % Align the MSSs to OCM
    LISTv_OriginTarg_MSS=LISTv_Origin_OCM;
    LISTv_XdircTarg_MSS=LISTv_Xdirc_OCM;   LISTv_YdircTarg_MSS=LISTv_Ydirc_OCM;
    LISTv_ZdircTarg_MSS=LISTv_Zdirc_OCM;
    LISTv_XdircTarg_MSS=LISTv_XdircTarg_MSS/norm(LISTv_XdircTarg_MSS);
    LISTv_YdircTarg_MSS=LISTv_YdircTarg_MSS/norm(LISTv_YdircTarg_MSS);
    LISTv_ZdircTarg_MSS=LISTv_ZdircTarg_MSS/norm(LISTv_ZdircTarg_MSS);

    LISTv_OriginTarg_MSS2=LISTv_Origin_OCM2;
    LISTv_XdircTarg_MSS2=LISTv_Xdirc_OCM2;   LISTv_YdircTarg_MSS2=LISTv_Ydirc_OCM2;
    LISTv_ZdircTarg_MSS2=LISTv_Zdirc_OCM2;
    LISTv_XdircTarg_MSS2=LISTv_XdircTarg_MSS2/norm(LISTv_XdircTarg_MSS2);
    LISTv_YdircTarg_MSS2=LISTv_YdircTarg_MSS2/norm(LISTv_YdircTarg_MSS2);
    LISTv_ZdircTarg_MSS2=LISTv_ZdircTarg_MSS2/norm(LISTv_ZdircTarg_MSS2);

    for i=1:size(MATX_PointsMot_MSS,1)
        MATX_PointsMot_MSS_(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_PointsMot_MSS(i,:)'-...
                                               LISTv_Origin_MSS));
        MATX_PointsMot_MSS2_(i,:)=transpose(inv([LISTv_Xdirc_MSS2,LISTv_Ydirc_MSS2,LISTv_Zdirc_MSS2])*(MATX_PointsMot_MSS(i,:)'-...
                                               LISTv_Origin_MSS2));

        MATX_PointsMot_MSS_(i,:)=transpose([LISTv_XdircTarg_MSS,...
                                            LISTv_YdircTarg_MSS,...
                                            LISTv_ZdircTarg_MSS]*MATX_PointsMot_MSS_(i,:)'+LISTv_OriginTarg_MSS);
        MATX_PointsMot_MSS2_(i,:)=transpose([LISTv_XdircTarg_MSS2,...
                                            LISTv_YdircTarg_MSS2,...
                                            LISTv_ZdircTarg_MSS2]*MATX_PointsMot_MSS2_(i,:)'+LISTv_OriginTarg_MSS2);
    end

%    tool.plot4RStruc(MATX_PointsMot_OCM,'','c')
%    tool.plot4RStruc(MATX_PointsMot_MSS_,'','m')
%    tool.plot4RStruc(MATX_PointsMot_MSS2_,'','g')
%    pause(0.1)
%    hold off

    TSOR_PointsMot_OCM(:,:,INDX_OCM)=MATX_PointsMot_OCM;
    TSOR_PointsMot_MSS(:,:,INDX_MSS)=MATX_PointsMot_MSS_;
    TSOR_PointsMot_MSS2(:,:,INDX_MSS)=MATX_PointsMot_MSS2_;
end

figure(3);clf;
%subplot(1,2,1);
tool.plot4RStruc(MATX_PointsMot_OCM,'','c')
%plot3([1,0],[0,0],[0,0],'r'); hold on;
%plot3([0,0],[1,0],[0,0],'g'); hold on;
%plot3([0,0],[0,0],[1,0],'b'); hold on;
%subplot(1,2,2);
tool.plot4RStruc(MATX_PointsMot_MSS_,'','m')
tool.plot4RStruc(MATX_PointsMot_MSS2_,'','g')
%plot3([LISTv_Origin_MSS2(1),LISTv_Origin_MSS2(1)+LISTv_Xdirc_MSS2(1)*30],...
%      [LISTv_Origin_MSS2(2),LISTv_Origin_MSS2(2)+LISTv_Xdirc_MSS2(2)*30],...
%      [LISTv_Origin_MSS2(3),LISTv_Origin_MSS2(3)+LISTv_Xdirc_MSS2(3)*30],'r'); hold on;
%plot3([LISTv_Origin_MSS2(1),LISTv_Origin_MSS2(1)+LISTv_Ydirc_MSS2(1)*30],...
%      [LISTv_Origin_MSS2(2),LISTv_Origin_MSS2(2)+LISTv_Ydirc_MSS2(2)*30],...
%      [LISTv_Origin_MSS2(3),LISTv_Origin_MSS2(3)+LISTv_Ydirc_MSS2(3)*30],'g'); hold on;
%plot3([LISTv_Origin_MSS2(1),LISTv_Origin_MSS2(1)+LISTv_Zdirc_MSS2(1)*30],...
%      [LISTv_Origin_MSS2(2),LISTv_Origin_MSS2(2)+LISTv_Zdirc_MSS2(2)*30],...
%      [LISTv_Origin_MSS2(3),LISTv_Origin_MSS2(3)+LISTv_Zdirc_MSS2(3)*30],'b'); hold on;


for INDX_OCM=1:(INDX_init_MSS-1)
    INDX_MSS=size(TSOR_PointsMot_MSS,3)-INDX_OCM+1;

    LISTv_VecHinge1_OCM=TSOR_PointsMot_OCM_(1+4,:,INDX_OCM)'-TSOR_PointsMot_OCM_(1,:,INDX_OCM)';
    LISTv_VecHinge1_OCM=LISTv_VecHinge1_OCM/norm(LISTv_VecHinge1_OCM);
    LISTv_VecHinge2_OCM=TSOR_PointsMot_OCM_(2+4,:,INDX_OCM)'-TSOR_PointsMot_OCM_(2,:,INDX_OCM)';
    LISTv_VecHinge2_OCM=LISTv_VecHinge2_OCM/norm(LISTv_VecHinge2_OCM);

    LISTv_RefHinge1_OCM=TSOR_PointsMot_OCM_(1,:,INDX_OCM)';
    LISTv_RefHinge1R_OCM=TSOR_PointsMot_OCM_(1,:,INDX_OCM+1)';
    LISTv_RefHinge2_OCM=TSOR_PointsMot_OCM_(2,:,INDX_OCM)';
    LISTv_RefHinge2R_OCM=TSOR_PointsMot_OCM_(2,:,INDX_OCM+1)';

    LISTv_VecMotHinge1_OCM=LISTv_RefHinge1R_OCM-LISTv_RefHinge1_OCM;
    LISTv_VecMotHinge1_OCM=LISTv_VecMotHinge1_OCM/norm(LISTv_VecMotHinge1_OCM);
    LISTv_VecMotHinge2_OCM=LISTv_RefHinge2R_OCM-LISTv_RefHinge2_OCM;
    LISTv_VecMotHinge2_OCM=LISTv_VecMotHinge2_OCM/norm(LISTv_VecMotHinge2_OCM);

    rho1=calculate_angle_from_vector_to_vector(LISTv_VecHinge1_OCM,LISTv_VecMotHinge1_OCM);
    rho2=calculate_angle_from_vector_to_vector(LISTv_VecHinge2_OCM,LISTv_VecMotHinge2_OCM);

    rho1_=calculate_angle_from_vector_to_vector([0;0;1],LISTv_VecMotHinge1_OCM);
    rho2_=calculate_angle_from_vector_to_vector([0;0;1],LISTv_VecMotHinge2_OCM);

    if rho1<=0
            rho1=-rho1;
        end
    if rho2<=0
            rho2=-rho2;
    end
    LISTr_rho1(INDX_OCM)=rho1;  LISTr_rho2(INDX_OCM)=rho2;

    if rho1_<=0
            rho1_=-rho1_;
        end
    if rho2_<=0
            rho2_=-rho2_;
    end
    LISTr_rho1_(INDX_OCM)=rho1_;  LISTr_rho2_(INDX_OCM)=rho2_;

    p1AVG=1/2*(TSOR_PointsMot_MSS(1,:,INDX_MSS)'+TSOR_PointsMot_MSS(1+4,:,INDX_MSS)')/100;
    p2AVG=1/2*(TSOR_PointsMot_MSS(2,:,INDX_MSS)'+TSOR_PointsMot_MSS(2+4,:,INDX_MSS)')/100;
    p3AVG=1/2*(TSOR_PointsMot_MSS(3,:,INDX_MSS)'+TSOR_PointsMot_MSS(3+4,:,INDX_MSS)')/100;
    p4AVG=1/2*(TSOR_PointsMot_MSS(4,:,INDX_MSS)'+TSOR_PointsMot_MSS(4+4,:,INDX_MSS)')/100;

    p1AVG2=1/2*(TSOR_PointsMot_MSS2(1,:,INDX_MSS)'+TSOR_PointsMot_MSS2(1+4,:,INDX_MSS)')/100;
    p2AVG2=1/2*(TSOR_PointsMot_MSS2(2,:,INDX_MSS)'+TSOR_PointsMot_MSS2(2+4,:,INDX_MSS)')/100;
    p3AVG2=1/2*(TSOR_PointsMot_MSS2(3,:,INDX_MSS)'+TSOR_PointsMot_MSS2(3+4,:,INDX_MSS)')/100;
    p4AVG2=1/2*(TSOR_PointsMot_MSS2(4,:,INDX_MSS)'+TSOR_PointsMot_MSS2(4+4,:,INDX_MSS)')/100;

    areaMSS(INDX_OCM)=calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],LISTv_VecMotHinge1_OCM);
    areaMSS2(INDX_OCM)=calProjectionArea([p1AVG2';p2AVG2';p3AVG2';p4AVG2'],LISTv_VecMotHinge2_OCM);

    areaMSS_(INDX_OCM)=sign(rho1_-pi/2)*calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],[0;0;-1]);
    areaMSS2_(INDX_OCM)=sign(rho2_-pi/2)*calProjectionArea([p1AVG2';p2AVG2';p3AVG2';p4AVG2'],[0;0;-1]);

    LISTv_VecMotHinge1_OCM=LISTv_RefHinge1R_OCM-LISTv_RefHinge1_OCM;
    LISTv_VecMotHinge2_OCM=LISTv_RefHinge2R_OCM-LISTv_RefHinge2_OCM;
    volumnMSS(INDX_OCM)=-calProjectionArea([p1AVG';p2AVG';p3AVG';p4AVG'],[0;0;-1])*LISTv_VecMotHinge1_OCM(3);
    volumnMSS_(INDX_OCM)=-calProjectionArea([p1AVG2';p2AVG2';p3AVG2';p4AVG2'],[0;0;-1])*LISTv_VecMotHinge2_OCM(3);

end

for i=1:length(volumnMSS)
    volumnMSS_Sum(i)=sum(volumnMSS(1:i));
    volumnMSS_Sum_(i)=sum(volumnMSS_(1:i));
end

figure(4);clf;
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),LISTr_rho1); hold on
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),LISTr_rho2); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),LISTr_rho1_); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),LISTr_rho2_); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),pi/2*ones(size(LISTr_rho1)));hold on
legend()
figure(5);clf;
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),areaMSS,'c--'); hold on
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),areaMSS2,'m--'); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),areaMSS_,'c-'); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),areaMSS2_,'m-'); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),areaMSS_+areaMSS2_,'k-'); hold on
legend()

figure(6);clf;
tool.plot4RStruc(TSOR_PointsMot_OCM(:,:,INDX_OCM),'','c')
tool.plot4RStruc(TSOR_PointsMot_MSS(:,:,INDX_MSS),'','m')
tool.plot4RStruc(TSOR_PointsMot_MSS2(:,:,INDX_MSS),'','g')
plot3([LISTv_RefHinge1_OCM(1),LISTv_RefHinge1_OCM(1)+LISTv_VecMotHinge1_OCM(1)*30],...
      [LISTv_RefHinge1_OCM(2),LISTv_RefHinge1_OCM(2)+LISTv_VecMotHinge1_OCM(2)*30],...
      [LISTv_RefHinge1_OCM(3),LISTv_RefHinge1_OCM(3)+LISTv_VecMotHinge1_OCM(3)*30],'k--'); hold on;
plot3([LISTv_RefHinge2_OCM(1),LISTv_RefHinge2_OCM(1)+LISTv_VecMotHinge2_OCM(1)*30],...
      [LISTv_RefHinge2_OCM(2),LISTv_RefHinge2_OCM(2)+LISTv_VecMotHinge2_OCM(2)*30],...
      [LISTv_RefHinge2_OCM(3),LISTv_RefHinge2_OCM(3)+LISTv_VecMotHinge2_OCM(3)*30],'k--'); hold on;

figure(7);clf;
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),volumnMSS); hold on
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),volumnMSS_); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),volumnMSS_Sum); hold on
plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),volumnMSS_Sum_); hold on
legend()
%plot(LISTr_Theta1_OCM(1:(INDX_init_MSS-1)),pi/2*ones(size(LISTr_rho1)));hold on


















function MATX_Points_Motion=Simu_BennettLinkage(alpha,beta,theta1)
    tool=TOOLS_SharedFunction();

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

        MATX_Points_Motion=[MATX_PointsSO_Benn;MATX_PointsSP_Benn];
end

function TSOR_Points_Motion=Simu_BennettLinkage_(alpha,beta,theta1Range)
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