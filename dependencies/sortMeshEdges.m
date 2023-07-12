function [sortedEdgeElements,Ef] = sortMeshEdges(Ef)
    unsortedEdgeElements=Ef;
    sortedEdgeElements={};
    nLines=0;
    while size(unsortedEdgeElements,1)>0
        sortedLine=[];
        firstEdge=unsortedEdgeElements(1,:);
        sortedLine=[sortedLine;firstEdge];
        unsortedEdgeElements(1,:)=[];
        closureFlag=0;
        while ~closureFlag
            lastEdge=sortedLine(end,:);
            connectedEdgesB=find(unsortedEdgeElements(:,2)==lastEdge(2));
            unsortedEdgeElements(connectedEdgesB,:)=fliplr(unsortedEdgeElements(connectedEdgesB,:));
            connectedEdgesA=find(unsortedEdgeElements(:,1)==lastEdge(2));
            if ~isempty(connectedEdgesA)
                connectedEdge=connectedEdgesA(1);
                sortedLine=[sortedLine;unsortedEdgeElements(connectedEdge,:)];
                unsortedEdgeElements(connectedEdge,:)=[];
            else
%                 disp('One contact line sorted, moving to next line')
                closureFlag=1;
                nLines=nLines+1;
                sortedEdgeElements{nLines,1}=sortedLine;
            end
        end
    end
    if nLines>1
        Ef=cell2mat(sortedEdgeElements);
    else
        Ef=sortedEdgeElements{1};
    end
end