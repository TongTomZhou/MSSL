clear;clc;close all
tool=TOOLS_SharedFunction;

MATX_Points_Motion=readmatrix('MATX_Points_Motion_Bistable4R.csv');
j=1;
for i=1:8:size(MATX_Points_Motion,1)
    TSOR_Points_Motion(:,:,j)=MATX_Points_Motion(i:i+7,:);
    j=j+1;
end
TSOR_Points_Motion(:,:,103)=[]; TSOR_Points_Motion(:,:,102)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for i=1:size(TSOR_Points_Motion,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_Points_Motion(:,:,i));
    theta1=MATX_DHStruc(4,1);
    if theta1<-0
        theta1=theta1+2*pi;
    elseif theta1>pi
        theta1=theta1-2*pi;
    end

    LISTr_theta1(i)=theta1;
end
deltaTheta1=LISTr_theta1(end)-LISTr_theta1(1);

plot(LISTr_theta1)


%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
VALE_a=4; VALE_b=4;
LISTr_theta1OCM_C1=[15,75,120,150,165]/180*pi;
for i=1:length(LISTr_theta1OCM_C1)
    theta1OCM_C1= LISTr_theta1OCM_C1(i);
    theta1OCM_C2= theta1OCM_C1-deltaTheta1;

    plot([LISTr_theta1(1),LISTr_theta1(end)],[theta1OCM_C1,theta1OCM_C2]); hold on;
end
legend('15','75','120','150','165')
xlim([-0.1,3.5]); ylim([-3.1,3.1]);


%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
VALE_ScaleRatio_OCM=10;  VALE_ScaleRatio_MSS=4;
i_alpha=1; i_beta=1;

iii=3; MATX_Rcmax=[];
theta1OCMStart=LISTr_theta1OCM_C1(iii);
for alphaOCM=-pi:2*pi/9:pi
    for betaOCM=-pi:2*pi/9:pi
        rc=[];
        for ii=1:length(LISTr_theta1)
            theta1=LISTr_theta1(ii);
            theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

            MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

            % Transform the MSS
            LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
            LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
            LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

            LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                              (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                                (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                                (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
            LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
            LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                              (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                                (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                                (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
            LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

            LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
            LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
            LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
            LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

            for i=1:size(MATX_Points_MSS,1)
                MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                                 LISTv_Origin_MSS));
            end

            theta2OCM=2*atan(sin(betaOCM/2+alphaOCM/2)/sin(betaOCM/2-alphaOCM/2)/tan(theta1OCM/2));
            offset=0;
            MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                          alphaOCM,betaOCM,alphaOCM,betaOCM;...
                          offset*ones(1,4);...
                          theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                          1*ones(1,4);-1*ones(1,4)];
            [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
                      tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                                       MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
            % Slip the OCM
            LISTv_SlipDist_Hinge1=3;
            LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
            LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
            MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
            MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
            LISTv_SlipDist_Hinge3=-3;
            LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
            LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
            MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
            MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

            MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
            MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

            %% rotate the MSS to align the OCM
            LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
            LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
            LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
            LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

            LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                              (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                                (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                                (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
            LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
            LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                              (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                                (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                                (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
            LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
            LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
            LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

            LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                              (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                                (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                                (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
            LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
            LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                              (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                                (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                                (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
            LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
            LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
            LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

            LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
            LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
            LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
            LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

            LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

            for i=1:size(MATX_Points_MSS,1)
                MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                                   MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
                MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                                   MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
            end
            rc(ii)=norm(MATX_Points_Wing1(7,:)-MATX_Points_Wing2(7,:));
            
%             TSOR_Points_Wing1(:,:,ii)=MATX_Points_Wing1;    
%             TSOR_Points_Wing2(:,:,ii)=MATX_Points_Wing2;
%             TSOR_Points_OCM(:,:,ii)=[MATX_PointsSO_OCM;MATX_PointsSP_OCM];
    
%            tool.plot4RStruc([MATX_PointsSO_OCM;MATX_PointsSP_OCM],'','c');
%            tool.plot4RStruc(MATX_Points_Wing1,'','m');
%            tool.plot4RStruc(MATX_Points_Wing2,'','g');
%            pause(0.05);hold off;
        end
        rcmax(i_alpha,i_beta)=max(rc);
        MATX_Rcmax=[MATX_Rcmax;alphaOCM,betaOCM,max(rc)];

        rcTensor(i_alpha,i_beta,:)=rc;
        i_beta=i_beta+1;
    end
i_alpha=i_alpha+1;
end

%writematrix(MATX_Rcmax,"MATX_GripperAnalysis_Rcmax_165.csv");


figure(4)
i_alpha=1; i_beta=1;
for alphaOCM=-pi:2*pi/9:pi
    for betaOCM=-pi:2*pi/9:pi
        scatter(alphaOCM,betaOCM,[],rcmax(i_alpha,i_beta),'filled'); hold on
        i_beta=i_beta+1;
    end
    i_alpha=i_alpha+1;
end
%caxis([0,300])
colorbar
axis equal


%% Special case
figure(5);clf;
alphaSpec=1.415471991295120;    betaSpec=1.726120662294673;
rc_Spec=[];
theta1OCMStart=120/180*pi;
for ii=1:length(LISTr_theta1)
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
    rc_Spec(ii)=norm(MATX_Points_Wing1(7,:)-MATX_Points_Wing2(7,:));
    LISTr_theta1OCM(ii)=theta1OCM;
    
%     TSOR_Points_Wing1(:,:,ii)=MATX_Points_Wing1;    
%     TSOR_Points_Wing2(:,:,ii)=MATX_Points_Wing2;
%     TSOR_Points_OCM(:,:,ii)=[MATX_PointsSO_OCM;MATX_PointsSP_OCM];
            
%     tool.plot4RStruc([MATX_PointsSO_OCM;MATX_PointsSP_OCM],'','c');
%     tool.plot4RStruc(MATX_Points_Wing1,'','m');
%     tool.plot4RStruc(MATX_Points_Wing2,'','g');
%     pause(0.1);hold off;

end
plot(LISTr_theta1,rc_Spec)

INDX_RcMax=find(rc_Spec==max(rc_Spec));
for ii=INDX_RcMax
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
end
MATX_GripperPoints_ThisMoment=[MATX_PointsSO_OCM;MATX_PointsSP_OCM;MATX_Points_Wing1;MATX_Points_Wing2];
writematrix(MATX_GripperPoints_ThisMoment,"MATX_MATX_GripperPoints_ThisMoment_120.csv");

figure(6);
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(1:8,:),'','c');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(9:16,:),'','m');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(17:24,:),'','g');


%% Special case 2
figure(7);clf;
alphaSpec=1.415471991295120;    betaSpec=1.726120662294673;
rc_Spec=[];
theta1OCMStart=15/180*pi;
for ii=1:length(LISTr_theta1)
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
    rc_Spec(ii)=norm(MATX_Points_Wing1(7,:)-MATX_Points_Wing2(7,:));
    LISTr_theta1OCM(ii)=theta1OCM;
    
    TSOR_Points_Wing1(:,:,ii)=MATX_Points_Wing1;    
    TSOR_Points_Wing2(:,:,ii)=MATX_Points_Wing2;
    TSOR_Points_OCM(:,:,ii)=[MATX_PointsSO_OCM;MATX_PointsSP_OCM];

%    tool.plot4RStruc([MATX_PointsSO_OCM;MATX_PointsSP_OCM],'','c');
%    tool.plot4RStruc(MATX_Points_Wing1,'','m');
%    tool.plot4RStruc(MATX_Points_Wing2,'','g');
%    pause(0.1);hold off;

end
plot(LISTr_theta1,rc_Spec)

INDX_RcMax=find(rc_Spec==max(rc_Spec));
for ii=INDX_RcMax
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
end
MATX_GripperPoints_ThisMoment=[MATX_PointsSO_OCM;MATX_PointsSP_OCM;MATX_Points_Wing1;MATX_Points_Wing2];
writematrix(MATX_GripperPoints_ThisMoment,"MATX_MATX_GripperPoints_ThisMoment_15.csv");

figure(8);
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(1:8,:),'','c');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(9:16,:),'','m');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(17:24,:),'','g');


%% Special case 3
figure(9);clf;
alphaSpec=1.415471991295120;    betaSpec=1.726120662294673;
rc_Spec=[];
theta1OCMStart=165/180*pi;
for ii=1:length(LISTr_theta1)
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
    rc_Spec(ii)=norm(MATX_Points_Wing1(7,:)-MATX_Points_Wing2(7,:));
    LISTr_theta1OCM(ii)=theta1OCM;
    
%     TSOR_Points_Wing1(:,:,ii)=MATX_Points_Wing1;    
%     TSOR_Points_Wing2(:,:,ii)=MATX_Points_Wing2;
%     TSOR_Points_OCM(:,:,ii)=[MATX_PointsSO_OCM;MATX_PointsSP_OCM];

   tool.plot4RStruc([MATX_PointsSO_OCM;MATX_PointsSP_OCM],'','c');
   tool.plot4RStruc(MATX_Points_Wing1,'','m');
   tool.plot4RStruc(MATX_Points_Wing2,'','g');
   pause(0.1);hold off;

end
plot(LISTr_theta1,rc_Spec)

INDX_RcMax=find(rc_Spec==max(rc_Spec));
for ii=INDX_RcMax
    theta1=LISTr_theta1(ii);
    theta1OCM=theta1OCMStart-(theta1-LISTr_theta1(1));

    MATX_Points_MSS=VALE_ScaleRatio_MSS*TSOR_Points_Motion(:,:,ii);

    % Transform the MSS
    LISTv_Origin_MSS=MATX_Points_MSS(5,:)';
    LISTv_Zdirc_MSS=MATX_Points_MSS(5,:)'-MATX_Points_MSS(1,:)';
    LISTv_Zdirc_MSS=LISTv_Zdirc_MSS/norm(LISTv_Zdirc_MSS);

    LISTv_Xdirc1_MSS=(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(6,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc1_MSS=LISTv_Xdirc1_MSS/norm(LISTv_Xdirc1_MSS);
    LISTv_Xdirc2_MSS=(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)'+...
                      (-(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))*(MATX_Points_MSS(4,:)'-MATX_Points_MSS(5,:)')/...
                        (norm(MATX_Points_MSS(1,:)-MATX_Points_MSS(5,:))^2))*...
                        (MATX_Points_MSS(1,:)'-MATX_Points_MSS(5,:)'));
    LISTv_Xdirc2_MSS=LISTv_Xdirc2_MSS/norm(LISTv_Xdirc2_MSS);

    LISTv_Xdirc_MSS=LISTv_Xdirc1_MSS+LISTv_Xdirc2_MSS;
    LISTv_Xdirc_MSS=LISTv_Xdirc_MSS/norm(LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=cross(LISTv_Zdirc_MSS,LISTv_Xdirc_MSS);
    LISTv_Ydirc_MSS=LISTv_Ydirc_MSS/norm(LISTv_Ydirc_MSS);

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_MSS(i,:)=transpose(inv([LISTv_Xdirc_MSS,LISTv_Ydirc_MSS,LISTv_Zdirc_MSS])*(MATX_Points_MSS(i,:)'-...
                                         LISTv_Origin_MSS));
    end

    theta2OCM=2*atan(sin(betaSpec/2+alphaSpec/2)/sin(betaSpec/2-alphaSpec/2)/tan(theta1OCM/2));
    offset=0;
    MATX_ParaOCM=[VALE_a,VALE_b,VALE_a,VALE_b;...
                  alphaSpec,betaSpec,alphaSpec,betaSpec;...
                  offset*ones(1,4);...
                  theta1OCM,theta2OCM,-theta1OCM,-theta2OCM;...
                  1*ones(1,4);-1*ones(1,4)];
    [~,~,MATX_PointsSO_OCM,MATX_PointsSP_OCM]=...
              tool.form4RStruc(MATX_ParaOCM(1,:),MATX_ParaOCM(2,:),MATX_ParaOCM(3,:),...
                               MATX_ParaOCM(5,:),MATX_ParaOCM(6,:),MATX_ParaOCM(4,:));
    % Slip the OCM
    LISTv_SlipDist_Hinge1=3;
    LISTv_Hinge1=MATX_PointsSP_OCM(1,:)-MATX_PointsSO_OCM(1,:);
    LISTv_Hinge1=LISTv_Hinge1/norm(LISTv_Hinge1);
    MATX_PointsSO_OCM(1,:)=MATX_PointsSO_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    MATX_PointsSP_OCM(1,:)=MATX_PointsSP_OCM(1,:)+LISTv_SlipDist_Hinge1*LISTv_Hinge1;
    LISTv_SlipDist_Hinge3=-3;
    LISTv_Hinge3=MATX_PointsSP_OCM(3,:)-MATX_PointsSO_OCM(3,:);
    LISTv_Hinge3=LISTv_Hinge3/norm(LISTv_Hinge3);
    MATX_PointsSO_OCM(3,:)=MATX_PointsSO_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;
    MATX_PointsSP_OCM(3,:)=MATX_PointsSP_OCM(3,:)+LISTv_SlipDist_Hinge3*LISTv_Hinge3;

    MATX_PointsSO_OCM=VALE_ScaleRatio_OCM*MATX_PointsSO_OCM;
    MATX_PointsSP_OCM=VALE_ScaleRatio_OCM*MATX_PointsSP_OCM;

    %% rotate the MSS to align the OCM
    LISTv_Zdirc_Wing1=MATX_PointsSP_OCM(1,:)'-MATX_PointsSO_OCM(1,:)';
    LISTv_Zdirc_Wing1=LISTv_Zdirc_Wing1/norm(LISTv_Zdirc_Wing1);
    LISTv_Zdirc_Wing2=MATX_PointsSO_OCM(3,:)'-MATX_PointsSP_OCM(3,:)';
    LISTv_Zdirc_Wing2=LISTv_Zdirc_Wing2/norm(LISTv_Zdirc_Wing2);

    LISTv_Xdirc1_Wing1=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc1_Wing1=LISTv_Xdirc1_Wing1/norm(LISTv_Xdirc1_Wing1);
    LISTv_Xdirc2_Wing1=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSP_OCM(1,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSP_OCM(1,:))*(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')/...
                        (norm(MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)')^2))*...
                        (MATX_PointsSO_OCM(1,:)'-MATX_PointsSP_OCM(1,:)'));
    LISTv_Xdirc2_Wing1=LISTv_Xdirc2_Wing1/norm(LISTv_Xdirc2_Wing1);
    LISTv_Xdirc_Wing1=LISTv_Xdirc1_Wing1+LISTv_Xdirc2_Wing1;
    LISTv_Xdirc_Wing1=LISTv_Xdirc_Wing1/norm(LISTv_Xdirc_Wing1);

    LISTv_Xdirc1_Wing2=(MATX_PointsSP_OCM(2,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSP_OCM(2,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc1_Wing2=LISTv_Xdirc1_Wing2/norm(LISTv_Xdirc1_Wing2);
    LISTv_Xdirc2_Wing2=(MATX_PointsSO_OCM(4,:)'-MATX_PointsSO_OCM(3,:)'+...
                      (-(MATX_PointsSO_OCM(4,:)-MATX_PointsSO_OCM(3,:))*(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')/...
                        (norm(MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)')^2))*...
                        (MATX_PointsSP_OCM(3,:)'-MATX_PointsSO_OCM(3,:)'));
    LISTv_Xdirc2_Wing2=LISTv_Xdirc2_Wing2/norm(LISTv_Xdirc2_Wing2);
    LISTv_Xdirc_Wing2=LISTv_Xdirc1_Wing2+LISTv_Xdirc2_Wing2;
    LISTv_Xdirc_Wing2=LISTv_Xdirc_Wing2/norm(LISTv_Xdirc_Wing2);

    LISTv_Ydirc_Wing1=cross(LISTv_Zdirc_Wing1,LISTv_Xdirc_Wing1);
    LISTv_Ydirc_Wing1=LISTv_Ydirc_Wing1/norm(LISTv_Ydirc_Wing1);
    LISTv_Ydirc_Wing2=cross(LISTv_Zdirc_Wing2,LISTv_Xdirc_Wing2);
    LISTv_Ydirc_Wing2=LISTv_Ydirc_Wing2/norm(LISTv_Ydirc_Wing2);

    LISTv_Origin_Wing1=MATX_PointsSP_OCM(1,:)'; LISTv_Origin_Wing2=MATX_PointsSO_OCM(3,:)';

    for i=1:size(MATX_Points_MSS,1)
        MATX_Points_Wing1(i,:)=transpose([LISTv_Xdirc_Wing1,LISTv_Ydirc_Wing1,LISTv_Zdirc_Wing1]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing1);
        MATX_Points_Wing2(i,:)=transpose([LISTv_Xdirc_Wing2,LISTv_Ydirc_Wing2,LISTv_Zdirc_Wing2]*...
                                           MATX_Points_MSS(i,:)'+LISTv_Origin_Wing2);
    end
end
MATX_GripperPoints_ThisMoment=[MATX_PointsSO_OCM;MATX_PointsSP_OCM;MATX_Points_Wing1;MATX_Points_Wing2];
writematrix(MATX_GripperPoints_ThisMoment,"MATX_MATX_GripperPoints_ThisMoment_165.csv");

figure(10);
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(1:8,:),'','c');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(9:16,:),'','m');
tool.plot4RStruc(MATX_GripperPoints_ThisMoment(17:24,:),'','g');