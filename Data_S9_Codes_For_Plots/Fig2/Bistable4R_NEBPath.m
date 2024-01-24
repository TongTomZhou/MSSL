clear;clc;%close all;
tool=TOOLS_SharedFunction;

load('RUST_SpecialCase_4R_Miao.mat');

MATX_Points_Start=MATX_Points_SpecialCase(1:8,:);
MATX_Points_End=MATX_Points_SpecialCase(9:16,:);

INDX_FixedNodes=[1,4,5,8];
topoBarsStruc=[1,5;2,6;3,7;4,8;1,2;1,6;5,2;5,6;...
               2,3;2,7;3,6;6,7;3,4;3,8;4,7;7,8;...
               1,4;4,5;1,8;5,8];

LISTv_Lengths_Ori=Lengths_Bars(MATX_Points_Start,topoBarsStruc);
[TSOR_Points_Motion,~]=ANS_Simulation_NEB(MATX_Points_Start,MATX_Points_End,INDX_FixedNodes,topoBarsStruc);

for i=1:size(TSOR_Points_Motion,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_Points_Motion(:,:,i));
    theta1=MATX_DHStruc(4,1); theta2=MATX_DHStruc(4,2);theta3=MATX_DHStruc(4,3); theta4=MATX_DHStruc(4,4);
    LISTr_Theta1(i)=theta1; LISTr_Theta2(i)=theta2; LISTr_Theta3(i)=theta3; LISTr_Theta4(i)=theta4;

    % Energy analysis
    LISTv_Lengths_i=Lengths_Bars(TSOR_Points_Motion(:,:,i),topoBarsStruc);
    LISTv_DeltaLengths_i=0;
    for ii=1:length(LISTv_Lengths_i)
        LISTv_DeltaLengths_i=LISTv_DeltaLengths_i+(LISTv_Lengths_i(ii)-LISTv_Lengths_Ori(ii))^2;
    end
    LISTv_DeltaLengths(i)=LISTv_DeltaLengths_i;
end

figure(1)
plot([LISTr_Theta1',LISTr_Theta2',LISTr_Theta3',LISTr_Theta4']);
figure(2)
plot(LISTv_DeltaLengths/max(LISTv_DeltaLengths),'r--');
figure(3)
plotEnergy=plot(LISTr_Theta3,LISTv_DeltaLengths/max(LISTv_DeltaLengths),'-'); hold on
plotEnergy.LineWidth=3; plotEnergy.Color='#4B6EB2';
plotIncomp.LineWidth=3; plotIncomp.Color='#EE6F42';
xlim([-2.5,0.5]);   ylim([0,1.2]);
legend('Energy','Incompability')
grid on

figure(4)
for i=1:size(TSOR_Points_Motion,3)
    tool.plot4RStruc(TSOR_Points_Motion(:,:,i),1,'m')
    pause(0.1)
    hold off
end

for i=1:size(TSOR_Points_Motion,3)
    MATX_Points_Motion(8*i-7:8*i,:)=TSOR_Points_Motion(:,:,i);
end
writematrix(MATX_Points_Motion,'MATX_Points_Motion_Bistable4R.csv');



function [TSOR_Points_Motion,fPhyProj]=ANS_Simulation_NEB(MATX_Points_Start,MATX_Points_End,INDX_FixedNodes,topoBarsStruc)
    tool=TOOLS_SharedFunction;

    numNodes=size(MATX_Points_Start,1);   numDirec=3*numNodes;    numInvBars=size(topoBarsStruc,1);
    configC1=reshape(MATX_Points_Start',[],1);   configC2=reshape(MATX_Points_End',[],1);
    K_nod=10^(0.1)*ones(numInvBars,1);   K_neb=10^(4.18); %0.1,1.5
    step_m=0.01;
    limitIter=5000; limitError=1e-3;
    lengthOri=Lengths_Bars(MATX_Points_Start,topoBarsStruc);

    position_zeros=zeros(1,numDirec);
    for i=1:length(INDX_FixedNodes)
        position_zeros(1,3*INDX_FixedNodes(i)-2:3*INDX_FixedNodes(i))=[1,1,1];
    end
    nodesB=reshape(position_zeros,3,[])';

    % Arrange the initial location of balls
    numConfigBalls=200; diffConfig=configC2-configC1; locListConfigBall=[];
    for i=1:numConfigBalls
        locConfigBalli=configC1+i/(numConfigBalls+1)*diffConfig;
        locListConfigBall=[locListConfigBall;locConfigBalli];
    end

    % Iteration for the minimum location of configuration balls
    eTolList=[]; numCount=1;
    while 1
        % Collect the Energy
        ePhy=0; ePhyList=[];
        for i=1:numConfigBalls
            locConfigBalli=locListConfigBall((i-1)*numDirec+1:i*numDirec);
            nodesi=reshape(locConfigBalli,3,[])';
            lengthi=Lengths_Bars(nodesi,topoBarsStruc);

            ePhyi=1/2*K_nod'*(lengthi-lengthOri);
            ePhy=ePhy+ePhyi; ePhyList=[ePhyList;ePhyi];
        end

        eNeb=0; eNebList=[];
        locListConfigBall_aC1aC2=[configC1;locListConfigBall;configC2];
        for i=1:numConfigBalls+1
            locConfigBalli=locListConfigBall_aC1aC2((i-1)*numDirec+1:i*numDirec);
            locConfigBalli_1=locListConfigBall_aC1aC2((i)*numDirec+1:(i+1)*numDirec);
            eNeb=eNeb+1/2*K_neb*norm(locConfigBalli_1-locConfigBalli)^2;
        end

        eTol=norm(ePhy)+norm(eNeb); eTolList=[eTolList;eTol];

        % Calculate the Excessive Force
        fPhy=[];
        for i=1:numConfigBalls
            locConfigBalli=locListConfigBall((i-1)*numDirec+1:i*numDirec);
            nodesi=reshape(locConfigBalli,3,[])';
            nodesAll=[nodesi,nodesB];

            Matrix_H=tool.eqtr3m4Bars(nodesAll,topoBarsStruc);
            lengthi=Lengths_Bars(nodesi,topoBarsStruc);
            ti=K_nod.*(lengthi-lengthOri); fi=Matrix_H*ti;

            for ii=1:numDirec
                if position_zeros(ii)==1
                    if ii==1
                        fi=[0;fi];
                    else
                        fi=[fi(1:ii-1);0;fi(ii:end)];
                    end
                end
            end
            fPhy=[fPhy;fi];
        end
        fNeb=[];
        for i=2:numConfigBalls+1
            locConfigBalli_m1=locListConfigBall_aC1aC2((i-2)*numDirec+1:(i-1)*numDirec);
            locConfigBalli=locListConfigBall_aC1aC2((i-1)*numDirec+1:(i)*numDirec);
            locConfigBalli_p1=locListConfigBall_aC1aC2((i)*numDirec+1:(i+1)*numDirec);
            fi=-K_neb*(locConfigBalli_m1+locConfigBalli_p1-2*locConfigBalli);

            for ii=1:numDirec
                if position_zeros(ii)==1
                    fi(ii)=0;
                end
            end
            fNeb=[fNeb;fi];
        end
        fTol=fPhy+fNeb;

        % Project the Excessive Force
        taoList=[]; %proj direction
        fPhyProj=[];    fNebProj=[];    ePhyList=[0;ePhyList;0];
        for i=2:numConfigBalls+1
            locConfigBalli_m1=locListConfigBall_aC1aC2((i-2)*numDirec+1:(i-1)*numDirec);
            locConfigBalli=locListConfigBall_aC1aC2((i-1)*numDirec+1:(i)*numDirec);
            locConfigBalli_p1=locListConfigBall_aC1aC2((i)*numDirec+1:(i+1)*numDirec);
            ePhyi_m1=ePhyList(i-1); ePhyi_1=ePhyList(i); ePhyi_p1=ePhyList(i+1);

            if ePhyi_p1>ePhyi_1 && ePhyi_1>ePhyi_m1
                tao_i=locConfigBalli_p1-locConfigBalli;
            elseif ePhyi_p1<ePhyi_1 && ePhyi_1<ePhyi_m1
                tao_i=locConfigBalli-locConfigBalli_m1;
            else
                del_E_max=max(norm(ePhyi_p1-ePhyi_1),norm(ePhyi_m1-ePhyi_1));
                del_E_min=min(norm(ePhyi_p1-ePhyi_1),norm(ePhyi_m1-ePhyi_1));
                del_x_left=locConfigBalli_p1-locConfigBalli;
                del_x_right=locConfigBalli-locConfigBalli_m1;
                if ePhyi_p1>ePhyi_m1
                    tao_i=del_x_left*del_E_max+del_x_right*del_E_min;
                else
                    tao_i=del_x_left*del_E_min+del_x_right*del_E_max;
                end
            end
            tao_i=tao_i/norm(tao_i);
            taoList=[taoList;tao_i];

            fPhyi=fPhy((i-2)*numDirec+1:(i-1)*numDirec);
            fPhyProj_i=fPhyi-(tao_i'*fPhyi)*tao_i;
            for ii=numDirec:-1:1 %remove the bounded nodes
                if position_zeros(ii)==1
                    fPhyProj_i(ii)=[];
                end
            end
            fPhyProj=[fPhyProj;fPhyProj_i];

            fNebi=fNeb((i-2)*numDirec+1:(i-1)*numDirec);
            fNebProj_i=(tao_i'*fNebi)*tao_i;
            fNebProj_i'*fNebi;%test>0
            for ii=numDirec:-1:1%remove the bounded nodes
                if position_zeros(ii)==1
                    fNebProj_i(ii)=[];
                end
            end
            fNebProj=[fNebProj;fNebProj_i];
        end
        [dim_reduced,~]=size(fNebProj_i);
        fTolProj=fPhyProj+fNebProj; % if not alike, adjust the K

        fprintf("At NEB Iterat. %d, Track if equal:\n \t Physical Force: %2.4f\n  \t NEB Force: %2.4f\n",...
                numCount,norm(fPhyProj),norm(fNebProj));
        numCount=numCount+1;

        if abs(norm(fPhyProj)-norm(fNebProj))<limitError
            fprintf("Converged!  \n");
            break
        elseif numCount>limitIter
            fprintf("Reach to the Iterative LimitL: %d. Stop and Restart.\n",limitIter);
            break
        else
            d_c=-fTolProj;    d_c=d_c/norm(d_c);
            d_c_a=[];
            for i=1:numConfigBalls
                d_c_a_i=d_c((i-1)*dim_reduced+1:i*dim_reduced);
                for ii=1:numDirec
                    if position_zeros(ii)==1
                        if ii==1
                            d_c_a_i=[0;d_c_a_i];
                        else
                            d_c_a_i=[d_c_a_i(1:ii-1);0;d_c_a_i(ii:end)];
                        end
                    end
                end
                d_c_a=[d_c_a;d_c_a_i];
            end
            locListConfigBall=locListConfigBall+d_c_a*step_m;
        end
    end

    for i=1:numConfigBalls
        locConfigBalli=locListConfigBall((i-1)*numDirec+1:i*numDirec);
        nodesi=reshape(locConfigBalli,3,[])';
        TSOR_Points_Motion_temp(:,:,i)=nodesi;
    end
    TSOR_Points_Motion(:,:,1)=MATX_Points_Start;
    TSOR_Points_Motion(:,:,2:numConfigBalls+1)=TSOR_Points_Motion_temp;
    TSOR_Points_Motion(:,:,numConfigBalls+2)=MATX_Points_End;
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