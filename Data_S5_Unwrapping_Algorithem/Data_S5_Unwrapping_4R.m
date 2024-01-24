



function [MATX_Points_Flat,MATX_Hinges]=FlattenTetra(MATX_Points_Tetra)

    tool=TOOLS_SharedFunction;
    MATX_Topo_Surface=[1,2,3;1,2,4;1,3,4;2,3,4];
    MATX_Topo_Hinges=[1,2;3,4];
    INDX_BottomNodes=[1,2,3];   INDX_RotaNodes=4;

    MATX_Topo_RotaPlate=[INDX_BottomNodes(1),INDX_BottomNodes(2),INDX_RotaNodes,INDX_BottomNodes(3);...
                         INDX_BottomNodes(2),INDX_BottomNodes(3),INDX_RotaNodes,INDX_BottomNodes(1);...
                         INDX_BottomNodes(3),INDX_BottomNodes(1),INDX_RotaNodes,INDX_BottomNodes(2)];

    for i=1:size(MATX_Topo_RotaPlate,1)
        MATX_Points_Rotai=MATX_Points_Tetra(MATX_Topo_RotaPlate(i,:),:);
        VALE_RotaAnglei=tool.CalProjANG(MATX_Points_Tetra(MATX_Topo_RotaPlate(i,3),:)',...
                                        MATX_Points_Tetra(MATX_Topo_RotaPlate(i,4),:)',...
                                        MATX_Points_Tetra(MATX_Topo_RotaPlate(i,1),:)',...
                                        MATX_Points_Tetra(MATX_Topo_RotaPlate(i,2),:)');
        VALE_RotaAnglei=pi-VALE_RotaAnglei;
        MATX_Points_Rota(i,:)=transpose(tool.RotaPoints(MATX_Points_Tetra(INDX_RotaNodes,:)',...
                                                        MATX_Points_Tetra(MATX_Topo_RotaPlate(i,1),:)',...
                                                        MATX_Points_Tetra(MATX_Topo_RotaPlate(i,2),:)',
                                                        VALE_RotaAnglei));
    end

    MATX_Points_Flat=[MATX_Points_Tetra(INDX_BottomNodes,:);MATX_Points_Rota];
    MATX_Hinges=[1,2;3];
end

