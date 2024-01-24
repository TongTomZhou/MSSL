clear;clc;close all;
tool=TOOLS_SharedFunction();

papercolormap=tool.customcolormap([0 0.25 0.5 0.75 1], {'#DC4638','#FDB26E','#E7F1D6','#7AB0D6','#3D4DA2'});
papercolormapRed =tool.customcolormap(linspace(0,1,4),{'#AA1E29','#DC4638','#FDB26E','#E7F1D6'});
papercolormapBlue=tool.customcolormap(linspace(0,1,4),{'#3E51A1','#3D4DA2','#7AB0D6','#E7F1D6'});

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

for i=1:size(TSOR_Points_Motion,3)
    MATX_DHStruc=tool.Nodes2DH_4R(TSOR_Points_Motion(:,:,i));
    theta3=-MATX_DHStruc(4,3);
    LISTr_theta3(i)=theta3;
end

for i=1:size(TSOR_Points_Motion,3)
    p1AVG=1/2*(TSOR_Points_Motion(1,:,i)'+TSOR_Points_Motion(1+4,:,i)')/100;
    p2AVG=1/2*(TSOR_Points_Motion(2,:,i)'+TSOR_Points_Motion(2+4,:,i)')/100;
    p3AVG=1/2*(TSOR_Points_Motion(3,:,i)'+TSOR_Points_Motion(3+4,:,i)')/100;
    p4AVG=1/2*(TSOR_Points_Motion(4,:,i)'+TSOR_Points_Motion(4+4,:,i)')/100;

    NodesProj=calProjectionPoints([p1AVG';p2AVG';p3AVG';p4AVG'],...
                                   TSOR_Points_Motion(3,:,i)'-TSOR_Points_Motion(3+4,:,i)');
    TSOR_PotProj_Motion(:,:,i)=NodesProj;

%    plot3([NodesProj(1,1),NodesProj(2,1)],[NodesProj(1,2),NodesProj(2,2)],[NodesProj(1,3),NodesProj(2,3)],'m-');
%    hold on
%    plot3([NodesProj(2,1),NodesProj(3,1)],[NodesProj(2,2),NodesProj(3,2)],[NodesProj(2,3),NodesProj(3,3)],'m-');
%    hold on
%    plot3([NodesProj(3,1),NodesProj(4,1)],[NodesProj(3,2),NodesProj(4,2)],[NodesProj(3,3),NodesProj(4,3)],'m-');
%    hold on
%    plot3([NodesProj(1,1),NodesProj(4,1)],[NodesProj(1,2),NodesProj(4,2)],[NodesProj(1,3),NodesProj(4,3)],'m-');
%    hold on
%
%    xlabel('X');    ylabel('Y');    zlabel('Z')
%    xlim([-0.3,0.1]); ylim([-0.3,0.15]); zlim([-0.05,0.2]);
%
%    pause(0.1)
%    hold off
end

figure(1); subplot(1,2,1);
i=202;
NodesProj=TSOR_PotProj_Motion(:,:,i);
%Nodes=TSOR_Points_Motion(:,:,i);
vecOrigin=NodesProj(3,:)';
vecNormZ=cross(NodesProj(3,:)'-NodesProj(1,:)',NodesProj(2,:)'-NodesProj(1,:)')+NodesProj(3,:)';
vecNormX=NodesProj(4,:)';
for i=1:size(NodesProj,1)
    NodesProj_(i,:)=transpose(tool.mapGlobal(NodesProj(i,:)',[vecOrigin';vecNormZ';vecNormX'],'-1'));
end
plot([NodesProj_(1,2),NodesProj_(2,2)],[NodesProj_(1,3),NodesProj_(2,3)],'k-');
hold on
plot([NodesProj_(2,2),NodesProj_(3,2)],[NodesProj_(2,3),NodesProj_(3,3)],'k-');
hold on
plot([NodesProj_(3,2),NodesProj_(4,2)],[NodesProj_(3,3),NodesProj_(4,3)],'k-');
hold on
plot([NodesProj_(1,2),NodesProj_(4,2)],[NodesProj_(1,3),NodesProj_(4,3)],'k-');
hold on

i=1;
NodesProj=TSOR_PotProj_Motion(:,:,i);
%Nodes=TSOR_Points_Motion(:,:,i);
vecOrigin=NodesProj(3,:)';
vecNormZ=cross(NodesProj(3,:)'-NodesProj(1,:)',NodesProj(2,:)'-NodesProj(1,:)')+NodesProj(3,:)';
vecNormX=NodesProj(4,:)';
for i=1:size(NodesProj,1)
    NodesProj_(i,:)=transpose(tool.mapGlobal(NodesProj(i,:)',[vecOrigin';vecNormZ';vecNormX'],'-1'));
end
plot([NodesProj_(1,2),NodesProj_(2,2)],[NodesProj_(1,3),NodesProj_(2,3)],'k-');
hold on
plot([NodesProj_(2,2),NodesProj_(3,2)],[NodesProj_(2,3),NodesProj_(3,3)],'k-');
hold on
plot([NodesProj_(3,2),NodesProj_(4,2)],[NodesProj_(3,3),NodesProj_(4,3)],'k-');
hold on
plot([NodesProj_(1,2),NodesProj_(4,2)],[NodesProj_(1,3),NodesProj_(4,3)],'k-');
hold on

subplot(1,2,2)


for i=1:201
    NodesProj=TSOR_PotProj_Motion(:,:,i);   NodesProj1=TSOR_PotProj_Motion(:,:,i+1);

    vecOriginj=NodesProj(3,:)';
    vecNormZj=cross(NodesProj(3,:)'-NodesProj(1,:)',NodesProj(2,:)'-NodesProj(1,:)')+NodesProj(3,:)';
    vecNormXj=NodesProj(4,:)';
    vecOriginj1=NodesProj1(3,:)';
    vecNormZj1=cross(NodesProj1(3,:)'-NodesProj1(1,:)',NodesProj1(2,:)'-NodesProj1(1,:)')+NodesProj1(3,:)';
    vecNormXj1=NodesProj1(4,:)';

    for ii=1:size(NodesProj,1)
        NodesProj_(ii,:)=transpose(tool.mapGlobal(NodesProj(ii,:)',[vecOriginj';vecNormZj';vecNormXj'],'-1'));
        NodesProj1_(ii,:)=transpose(tool.mapGlobal(NodesProj1(ii,:)',[vecOriginj1';vecNormZj1';vecNormXj1'],'-1'));
    end
    subplot(1,2,1);
    plot([NodesProj_(1,2),NodesProj1_(1,2)],[NodesProj_(1,3),NodesProj1_(1,3)],'c-');
    hold on
    plot([NodesProj_(2,2),NodesProj1_(2,2)],[NodesProj_(2,3),NodesProj1_(2,3)],'m-');
    hold on
    
    ovalityList(i)=calculateOvality(NodesProj(:,2:3));
    
end
axis equal
xlabel('X');    ylabel('Y');    zlabel('Z')
xlim([-0.15,0.2]); ylim([-0.25,0.05]);

    
subplot(1,2,2);
plot(LISTr_theta3(1:201),ovalityList); hold on
xlabel('X');    ylabel('Y');
xlim([-0.15,0.2]); ylim([-0.25,0.05]);

%gamma3List=[-pi*1/6,-pi*1/4,-pi*1/3,-pi*1/2,-pi*2/3,-pi*3/4,-pi*5/6,-pi];
gamma3List=-pi:pi/10:0;
for iii=1:length(gamma3List)
    gamma3=gamma3List(iii);
    area_union=[]; area_intersection=[];
    for ii=1:size(TSOR_PotProj_Motion,3)
        NodesProj=TSOR_PotProj_Motion(:,:,ii);
        vecOrigin=NodesProj(3,:)';
        vecNormZ=cross(NodesProj(3,:)'-NodesProj(1,:)',NodesProj(2,:)'-NodesProj(1,:)')+NodesProj(3,:)';
        vecNormX=NodesProj(4,:)';
        for i=1:size(NodesProj,1)
            NodesProj1(i,:)=transpose(tool.mapGlobal(NodesProj(i,:)',[vecOrigin';vecNormZ';vecNormX'],'-1'));
        end
        NodesProj_=NodesProj1(:,2:3);
        NodesProj__=tool.RotaPoints2(NodesProj_,[0;0;0],[0;0;1],gamma3);

        [x_union, y_union]=polybool('union', NodesProj_(:,1), NodesProj_(:,2), NodesProj__(:,1), NodesProj__(:,2));
        area_union_ = polyarea(x_union, y_union);
        if isnan(area_union_)
            [x_union, y_union]=polybool('union', NodesProj_(:,1), NodesProj_(:,2), NodesProj_(:,1), NodesProj_(:,2));
            area_union(ii)=2*polyarea(x_union, y_union);
        else
            area_union(ii)=area_union_;
        end
        [x_intersection, y_intersection]=polybool('intersection', NodesProj_(:,1), NodesProj_(:,2), NodesProj__(:,1), NodesProj__(:,2));
        area_intersection(ii) = polyarea(x_intersection, y_intersection);
    end
    MATX_Area_Union(:,iii)=area_union'; MATX_Area_Intersection(:,iii)=area_intersection';
end

[Gamma3List,LISTr_Theta3]=meshgrid(gamma3List, LISTr_theta3);
margin=0.1;

figure(2)
contourf(Gamma3List, LISTr_Theta3, MATX_Area_Union);
f2c=colorbar;
colormap(papercolormapBlue)
xticks([-pi,-5*pi/6,-2*pi/3,-pi/2,-pi/3,-pi/6,0]);
xticklabels({'-\pi','-5/6pi','-2/3pi','-1/2pi','-1/3pi','-1/6pi',0})
yticks([0,pi/6,pi/3,pi/2,2*pi/3]);
yticklabels({'0','1/6pi','1/3pi','1/2pi','2/3pi'})
xlim([-pi-margin,0+margin]);
ylim([0-margin,2*pi/3+margin/2]);
caxis([0,0.04]);
set(f2c,'YTick',0:0.01:0.04);
axis equal

figure(3)
contourf(Gamma3List, LISTr_Theta3, MATX_Area_Intersection);
f3c=colorbar;
colormap(papercolormapRed)
xticks([-pi,-5*pi/6,-2*pi/3,-pi/2,-pi/3,-pi/6,0]);
xticklabels({'-\pi','-5/6pi','-2/3pi','-1/2pi','-1/3pi','-1/6pi',0})
yticks([0,pi/6,pi/3,pi/2,2*pi/3]);
yticklabels({'0','1/6pi','1/3pi','1/2pi','2/3pi'})
xlim([-pi-margin,0+margin]);
ylim([0-margin,2*pi/3+margin/2]);
caxis([0,0.02]);
set(f3c,'YTick',0:0.005:0.02);
axis equal

figure(4)
contourf(Gamma3List, LISTr_Theta3, MATX_Area_Intersection./MATX_Area_Union);
f4c=colorbar;
colormap(papercolormap)
xticks([-pi,-5*pi/6,-2*pi/3,-pi/2,-pi/3,-pi/6,0]);
xticklabels({'-\pi','-5/6pi','-2/3pi','-1/2pi','-1/3pi','-1/6pi',0})
yticks([0,pi/6,pi/3,pi/2,2*pi/3]);
yticklabels({'0','1/6pi','1/3pi','1/2pi','2/3pi'})
xlim([-pi-margin,0+margin]);
ylim([0-margin,2*pi/3+margin/2]);
caxis([0,1]);
set(f4c,'YTick',0:0.25:1);
axis equal

figure(5)
gamma3List=-[3/4*pi-1*pi/10,3/4*pi-pi/20,3/4*pi,3/4*pi+pi/20,3/4*pi+pi/10;];
for iii=1:length(gamma3List)
    gamma3_Spec=gamma3List(iii);
    area_union=[]; area_intersection=[];
    for ii=1:size(TSOR_PotProj_Motion,3)
        NodesProj=TSOR_PotProj_Motion(:,:,ii);
        vecOrigin=NodesProj(3,:)';
        vecNormZ=cross(NodesProj(3,:)'-NodesProj(1,:)',NodesProj(2,:)'-NodesProj(1,:)')+NodesProj(3,:)';
        vecNormX=NodesProj(4,:)';
        for i=1:size(NodesProj,1)
            NodesProj1(i,:)=transpose(tool.mapGlobal(NodesProj(i,:)',[vecOrigin';vecNormZ';vecNormX'],'-1'));
        end
        NodesProj_=NodesProj1(:,2:3);
        NodesProj__=tool.RotaPoints2(NodesProj_,[0;0;0],[0;0;1],gamma3_Spec);

        [x_union, y_union]=polybool('union', NodesProj_(:,1), NodesProj_(:,2), NodesProj__(:,1), NodesProj__(:,2));
        area_union_ = polyarea(x_union, y_union);
        if isnan(area_union_)
            [x_union, y_union]=polybool('union', NodesProj_(:,1), NodesProj_(:,2), NodesProj_(:,1), NodesProj_(:,2));
            area_union(ii)=2*polyarea(x_union, y_union);
        else
            area_union(ii)=area_union_;
        end
        [x_intersection, y_intersection]=polybool('intersection', NodesProj_(:,1), NodesProj_(:,2), NodesProj__(:,1), NodesProj__(:,2));
        area_intersection(ii) = polyarea(x_intersection, y_intersection);
    end
    MATX_Area_Union_Spec(:,iii)=area_union'; MATX_Area_Intersection_Spec(:,iii)=area_intersection';
end

plot(LISTr_theta3,MATX_Area_Intersection_Spec(:,1)./MATX_Area_Union_Spec(:,1));hold on;
%plot(LISTr_theta3,MATX_Area_Intersection_Spec(:,2)./MATX_Area_Union_Spec(:,1));hold on;
%plot(LISTr_theta3,MATX_Area_Intersection_Spec(:,3)./MATX_Area_Union_Spec(:,1));hold on;
%plot(LISTr_theta3,MATX_Area_Intersection_Spec(:,4)./MATX_Area_Union_Spec(:,1));hold on;
%plot(LISTr_theta3,MATX_Area_Intersection_Spec(:,5)./MATX_Area_Union_Spec(:,1));hold on;
xticks([0,pi/6,pi/3,pi/2,2*pi/3]);
xticklabels({'0','1/6pi','1/3pi','1/2pi','2/3pi'})
ylim([0-margin/10,0.2+margin/10]);
xlim([0-margin,2*pi/3+margin/2]);


%plot(LISTr_Theta3,MATX_Area_Union(:,1),'r-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,2),'g-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,3),'b-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,4),'c-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,5),'m-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,6),'k-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,7),'y-'); hold on
%plot(LISTr_Theta3,MATX_Area_Union(:,8),'y--'); hold on
%
%plot(LISTr_Theta3,MATX_Area_Intersection(:,1),'r--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,2),'g--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,3),'b--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,4),'c--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,5),'m--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,6),'k--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,7),'y--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,8),'y--'); hold on
%
%figure(3)
%plot(LISTr_Theta3,MATX_Area_Intersection(:,1)./MATX_Area_Union(:,1),'r--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,2)./MATX_Area_Union(:,2),'g--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,3)./MATX_Area_Union(:,3),'b--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,4)./MATX_Area_Union(:,4),'c--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,5)./MATX_Area_Union(:,5),'m--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,6)./MATX_Area_Union(:,6),'k--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,7)./MATX_Area_Union(:,7),'y--'); hold on
%plot(LISTr_Theta3,MATX_Area_Intersection(:,8)./MATX_Area_Union(:,8),'y--'); hold on



function NodesProj=calProjectionPoints(NodesList,NormVec)
    n=NormVec/norm(NormVec);
    A=NodesList(1,:)';
    B=NodesList(2,:)'; C=NodesList(3,:)'; D=NodesList(4,:)';

    % Step 2: Project the vertices onto the plane
    proj_A = A - dot(A, n) * n;
    proj_B = B - dot(B, n) * n;
    proj_C = C - dot(C, n) * n;
    proj_D = D - dot(D, n) * n;

    NodesProj=[proj_A';proj_B';proj_C';proj_D'];
end

function ovality = calculateOvality(points)
    % Calculate the centroid
    Cx = mean(points(:, 1));
    Cy = mean(points(:, 2));

    % Calculate distances from each point to the centroid
    distances = sqrt((points(:, 1) - Cx).^2 + (points(:, 2) - Cy).^2);

    % Find the maximum distance
    Dmax = max(distances);

    % Calculate the average distance
    avgDistance = mean(distances);

    % Calculate the ovality
    ovality = Dmax / avgDistance;
end
