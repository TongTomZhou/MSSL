% PAPER NAME: Multi-stable spatial linkages
% THEME:      Data S1: Stiffness Analysis for 4R MSSL
% AUTHOR:     Tong Zhou <zhoutong_whu@whu.edu.cn>
% UPDATA:     2023-11-06

function eValStiff=calStiff4R(points,indexFixnodes,K)  % Main function
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
        K=[K;constElastic];
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