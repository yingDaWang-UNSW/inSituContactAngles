function [thetaMacroSun,thetaMacroKearney,directTheta,V,Ff,Fs,Effs,voxelInds,KCells,numUpsamples] = computeSmoothContactAnglesOnMeshv4(domain,minTriangles,numUpsamples,smoothFlag,smoothingOptVec,sunFlag,kearneyFlag,directFlag,curvaturesFlag,outlierFlag,fluids,solids,plotFlag)
    [thetaMacroSun,thetaMacroKearney,directTheta,V,Ff,Fs,Effs,KCells] = deal([],[],[],[],[],[],{},{});
    warning off
    if ~isstruct(domain)
%         CC2=bwconncomp(domain~=1);
%         chi=1-CC2.NumObjects+1;
        chi=imEuler3d(domain==1);
        [X,Y,Z]=meshgrid(1:size(domain,1),1:size(domain,2),1:size(domain,3));
        domain1=double(domain==1);
        
        [F,V] = MarchingCubes(X,Y,Z,permute(domain1,[2,1,3]),0.5);
        
%         [F,V] = isosurface(X,Y,Z,permute(domain1,[2,1,3]),0.5);
        
        nTri=0;
        while nTri<minTriangles
            nTri=size(F,1);
            if nTri<minTriangles
                [V,F] = subdivideMesh(V, F, 2);
%                 disp(['Upsampling Blob'])
                numUpsamples=numUpsamples+1;
            end
        end
        % extract the meshes
        % extract the fluid-solid interface
        domain1=double(domain>0);
        
        [Fs,Vs] = MarchingCubes(X,Y,Z,permute(domain1,[2,1,3]),0.5);view(3)
%         [Fs,Vs] = isosurface(X,Y,Z,permute(domain1,[2,1,3]),0.5);view(3)
        
        for nn=1:numUpsamples
            [Vs,Fs] = subdivideMesh(Vs, Fs, 2);
        end
        % find common nodes
        [~,nodef,~]=intersect(V,Vs,'rows');
        % find elements on the blob that only have 1 or 2 nodes intersecting the
        % solid surface
        Ff=F;
        Fs=F;
        % only keep elements that do not have all 3 corners in the intersection
        % region
        for elem=1:size(F,1)
            tri=F(elem,:);
            flags=[sum(tri(1)==nodef), sum(tri(2)==nodef), sum(tri(3)==nodef)];
            if sum(flags)==3
                Ff(elem,:)=[0 0 0];
            else
                Fs(elem,:)=[0 0 0];
            end
        end
        Ff( ~any(Ff,2), : ) = [];  %rows
        Fs( ~any(Fs,2), : ) = [];  %rows
        F=[Fs;Ff];
        voxelInds={};
        Effs={};
        if numel(Fs)<10 || numel(Ff)<10
            disp('Cluster not in contact with surface/fluid')
            return
        end
    else
        V=domain.V;
        Fs=domain.Fs;
        Ff=domain.Ff;
        F=domain.F;
        voxelInds=domain.voxelInds;
        Effs = domain.Effs;
        chi=1;
    end
    % find and sort the contact line nodes
    % find the edge nodes on mesh F
    Tf=triangulation(Ff,V);
    Ef = freeBoundary(Tf); % find the edge elements that are unbounded
    contactPerimeter=ceil(numel(Ef)./(2.^numUpsamples));
    removedLine=[];

    if isempty(Effs)
        [Effs,Ef] = sortMeshEdges(Ef);
        nContactLines=size(Effs,1);
        for n=1:nContactLines
            edgeNodeInds=unique(Effs{n}(:),'stable');
            X=V(edgeNodeInds,:);
            if norm(diff(X)) == 0 || size(X,1)<=3*2^(numUpsamples)
    %             disp('Contact line is duplicate points or is less than 1 voxel, removing line')
                removedLine=[removedLine;n];
            end
        end
    end
    if isempty(voxelInds)
%         disp('Mapping contact line nodes to voxels')
        for n=1:nContactLines
            edgeNodeInds=unique(Effs{n}(:),'stable');
            X=V(edgeNodeInds,:);            
            domainInds=ceil(X);
            domainInds=sub2ind(size(domain),domainInds(:,1),domainInds(:,2),domainInds(:,3));
            voxelInds{n,1}=domainInds;
        end
    end
    Effs(removedLine)=[];
    voxelInds(removedLine)=[];
    nContactLines=numel(Effs);

    if plotFlag
        figure(1);clf(1);        
        drawMesh(V,Fs);title('Original Mesh')
        hold on
        drawMesh(V,Ff,'b');view(3);axis equal tight;alpha(1)
        for n=1:nContactLines
            for e=1:size(Effs{n},1)
                edgeNodeInds=unique(Effs{n}(:),'stable');
                X=V(edgeNodeInds,:);
                plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',15);
                plot3(V(Effs{n}(e,:),1),V(Effs{n}(e,:),2),V(Effs{n}(e,:),3),'-g','LineWidth',3); 
            end
        end
    end
    if smoothFlag
        KNodes=zeros(size(V,1),3);
        if smoothingOptVec(1)>1
            nSurfIters=smoothingOptVec(1);
            itersPerLine=smoothingOptVec(2);
            shrinkIters=smoothingOptVec(3);
            expansionFactor=0.00;
            shrinkVolLim=0;
            expandVolLim=0;
        else
            nSurfIters=1e10;
            itersPerLine=smoothingOptVec(2);
            shrinkIters=1e10;
            shrinkVolLim=smoothingOptVec(3);
            expandVolLim=smoothingOptVec(1);
            expansionFactor=0.01;
            smoothingOptVec(4);
        end 
        lineShrinksPerVolumePreservation=10;
        surfShrinksPerVolumePreservation=10;
        volumePreservationFactorLine=1+expansionFactor;
        volumePreservationFactorSurface=1+expansionFactor;
        
        VOld=V;
        nv = size(V, 1);
        adjLine=meshAdjacencyMatrix(Ef);
        if size(adjLine, 1) < nv
            adjLine(nv, nv) = 0;
        end
        adjLine = adjLine + speye(nv);
        w = spdiags(full(sum(adjLine, 2).^(-1)), 0, nv, nv);
        adjLine = w * adjLine;

        adjSurf=meshAdjacencyMatrix(F);
        if size(adjSurf, 1) < nv
            adjSurf(nv, nv) = 0;
        end
        adjSurf = adjSurf + speye(nv);
        w = spdiags(full(sum(adjSurf, 2).^(-1)), 0, nv, nv);
        adjSurf = w * adjSurf;
        % 
        volVec=meshVolume(V,F);
        deformationIndex=[];
        for m=1:nSurfIters
            % smooth the contact line
            if mod(m,itersPerLine)==1 || itersPerLine == 1 %|| m<=5 
                Vl = adjLine*V;
                if m>shrinkIters
                    if lineShrinksPerVolumePreservation>1
                        for s=1:lineShrinksPerVolumePreservation
                            Vl = adjLine*Vl;
                        end
                    end
                    for n=1:nContactLines
                        edgeNodeInds=unique(Effs{n}(:),'stable');
                        X=Vl(edgeNodeInds,:);
                        X=[X(end,:);X;X(1,:)];
                        [LL,R,K] = curvature(X);
                        K=-K(2:end-1,:);
                        KNodes(edgeNodeInds,:)=K;
                    end
                    KNodes=KNodes./vecnorm(KNodes,2,2);
                    KNodes(isnan(KNodes))=0;
                    % if a few smoothing cycles have passed, activate volume
                    % preservation - calculate distance shift and move back in
                    % normal direction - scale distance by curvature
                    deltaDist=((Vl(:,1)-V(:,1)).^2+(Vl(:,2)-V(:,2)).^2+(Vl(:,3)-V(:,3)).^2).^0.5;
                    Vl = Vl + volumePreservationFactorLine.*deltaDist.*KNodes;
                end
            end
            % smooth the whole surface
            V2 = adjSurf*V;
            if m>shrinkIters
                if surfShrinksPerVolumePreservation>1
                    for s=1:surfShrinksPerVolumePreservation
                        V2 = adjSurf*V2;
                    end
                end
                KNodes = vertexNormal(V2,F);
                KNodes = KNodes./vecnorm(KNodes,2,2);
                KNodes(isnan(KNodes))=0;
                % if a few smoothing cycles have passed, activate volume
                % preservation - calculate distance shift and move back in
                % normal direction
                deltaDist=((V2(:,1)-V(:,1)).^2+(V2(:,2)-V(:,2)).^2+(V2(:,3)-V(:,3)).^2).^0.5;
                V2 = V2 + volumePreservationFactorSurface.*deltaDist.*KNodes;
%                 % confirm the movement was in the correct direction
%                 deltaDist2=((V2n(:,1)-V(:,1)).^2+(V2n(:,2)-V(:,2)).^2+(V2n(:,3)-V(:,3)).^2).^0.5;
%                 realDirections=(deltaDist-deltaDist2);
%                 realDirections=realDirections./abs(realDirections);
%                 realDirections(isnan(realDirections))=0;
%                 V2 = V2 + realDirections.*volumePreservationFactorSurface.*deltaDist.*KNodes;
            end
            V=V2;
            % overwrite the surface with the contact line
            if itersPerLine>0
                for n=1:nContactLines
                    edgeNodeInds=unique(Effs{n}(:),'stable');
                    V(edgeNodeInds,:)=Vl(edgeNodeInds,:);
                end
            end
        end
        % 
        Vol = meshVolume(V,F);
        volVec=[volVec;Vol];
        deformation=var(sum(sqrt((V-VOld).^2),2));
        deformationIndex=[deformationIndex;deformation];
    end
    faceAreaf= 1/2*vecnorm(cross(V(Ff(:,2),:)-V(Ff(:,1),:),V(Ff(:,3),:)-V(Ff(:,1),:))')';
    faceAreas= 1/2*vecnorm(cross(V(Fs(:,2),:)-V(Fs(:,1),:),V(Fs(:,3),:)-V(Fs(:,1),:))')';
    faceArea=[faceAreas;faceAreaf];
    if plotFlag
        figure(2);clf(2);title('Smoothed Mesh and Contact Line')
        drawMesh(V,Fs)
        hold on
        drawMesh(V,Ff,'b');view(3);axis equal tight;alpha(1)

        for n=1:nContactLines
            for e=1:size(Effs{n},1)
                edgeNodeInds=unique(Effs{n}(:),'stable');
                X=V(edgeNodeInds,:);
                plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',15);
                plot3(V(Effs{n}(e,:),1),V(Effs{n}(e,:),2),V(Effs{n}(e,:),3),'-g','LineWidth',3); 
            end
        end
    end
    if kearneyFlag || directFlag
               % get solid surface normals at contact edges
        FVs=triangulation(Fs,V);
        Ps = incenter(FVs);
        Ns = faceNormal(FVs);  
        Ns=-1.*Ns;
        if plotFlag
            figure(6);clf(6);
            trisurf(FVs, ...
                 'FaceColor','cyan','FaceAlpha',0.8);
            axis equal
            hold on  
            quiver3(Ps(:,1),Ps(:,2),Ps(:,3), ...
                 Ns(:,1),Ns(:,2),Ns(:,3),5,'color','r');
            title('Normals at Solid-Fluid Interface')
            axis equal tight
        end
        % interpolate normals from solid elements to contact line nodes
        if plotFlag
            figure(7);clf(7);
            trisurf(FVs);alpha(0.5);hold on;axis equal tight
            title('Solid-Fluid Interface Normals at Contact Line')
        end
        meanMeshFaceArea=mean(faceAreas);
        minTriangleArea=0;
        NsLineNodesVec=[];
        in=[];
        for n=1:nContactLines
            edgeNodeInds=unique(Effs{n}(:),'stable');
            X=V(edgeNodeInds,:);
            NsLineNodes=zeros(size(X));
            for node=1:numel(edgeNodeInds)
                attachedElems=vertexAttachments(FVs,edgeNodeInds(node));
                NsElems=Ns(attachedElems{:},:);
                PsElems=Ps(attachedElems{:},:);
                AsElems=faceAreas(attachedElems{:});
                AsElems(AsElems<minTriangleArea)=0;
                AsElems=AsElems./meanMeshFaceArea;
                NsLineNodes(node,:)=NsElems'*AsElems;
                if plotFlag
                    trisurf(FVs(attachedElems{:},:),V(:,1),V(:,2),V(:,3),'FaceColor','r');
                    quiver3(PsElems(:,1),PsElems(:,2),PsElems(:,3), NsElems(:,1).*AsElems,NsElems(:,2).*AsElems,NsElems(:,3).*AsElems,5,'color','g');
                end
            end
            NsLineNodes=NsLineNodes./vecnorm(NsLineNodes,2,2);
            for dim=1:3
                NsLineNodes(:,dim)=filloutliers(NsLineNodes(:,dim),'linear','quartiles');
            end
            for dim=1:3
                NsLineNodes(:,dim)=fillmissing(NsLineNodes(:,dim),'movmean',10);
            end               

            NsLineNodesVec=[NsLineNodesVec;NsLineNodes];
            if plotFlag
                quiver3(X(:,1),X(:,2),X(:,3),NsLineNodes(:,1),NsLineNodes(:,2),NsLineNodes(:,3),'color','b');
            end         
            % find which surface normal vectors point into the cluster and
            % which point out of the cluster - travel a very small distance
            % along the vector and check if inside or outside F,V  
            X2 = X + (NsLineNodes).*1e-5;
%                 scatter3(X2(:,1),X2(:,2),X2(:,3),'r');
            try
                in = [in;intriangulation(V,F,X)];
            catch
                in=[];
            end
        end
        NsLineNodesVec(isnan(NsLineNodesVec))=0; 
    end
    if sunFlag || kearneyFlag || curvaturesFlag
        thetaMacroSun=[];
        thetaMacroKearney=[];
        % measure curvatures on mesh
        FVs=struct();
        FVs.faces=Fs;
        FVs.vertices=V;
        FVf=struct();
        FVf.faces=Ff;
        FVf.vertices=V;
        [Cmeanf,Cgaussianf,~,~,Lambda1f,Lambda2f]=patchcurvature(FVf,0);
        [Cmeans,Cgaussians,~,~,Lambda1s,Lambda2s]=patchcurvature(FVs,0);
        if plotFlag
            figure(3);clf(3);
            trisurf(FVf.faces,V(:,1),V(:,2),V(:,3),log(abs(Cgaussianf)),'edgecolor','none');axis equal;axis tight;view(3);colormap jet
            hold on
            trisurf(FVs.faces,V(:,1),V(:,2),V(:,3),log(abs(Cgaussians)),'edgecolor','none');axis equal;axis tight;view(3);colormap jet
            colorbar
            title('AbsLog10 Curvatures')
            hold off
        end
        
        % interpolate curvatures from nodes to elements and list them out
        Kxyzf=zeros(size(FVf.faces,1),3);
        Kgevf=zeros(size(FVf.faces,1),1);
        Kmevf=zeros(size(FVf.faces,1),1);
        K1evf=zeros(size(FVf.faces,1),1);
        K2evf=zeros(size(FVf.faces,1),1);
        for elem=1:size(FVf.faces,1)
            nodes=FVf.vertices(FVf.faces(elem,:),:);
            centroid=mean(nodes,1);
            Kge=mean(Cgaussianf(FVf.faces(elem,:)));
            Kme=mean(Cmeanf(FVf.faces(elem,:)));
            K1e=mean(Lambda1f(FVf.faces(elem,:)));
            K2e=mean(Lambda2f(FVf.faces(elem,:)));
            if Kme<0 && Kge>0
                Kge=-1*Kge;
            end
            Kxyzf(elem,:)=centroid;
            Kgevf(elem)=Kge;
            Kmevf(elem)=Kme;
            K1evf(elem)=K1e;
            K2evf(elem)=K2e;
        end
        Kxyzs=zeros(size(FVs.faces,1),3);
        Kgevs=zeros(size(FVs.faces,1),1);
        Kmevs=zeros(size(FVs.faces,1),1);
        K1evs=zeros(size(FVs.faces,1),1);
        K2evs=zeros(size(FVs.faces,1),1);
        for elem=1:size(FVs.faces,1)
            nodes=FVs.vertices(FVs.faces(elem,:),:);
            centroid=mean(nodes,1);
            Kge=mean(Cgaussians(FVs.faces(elem,:)));
            Kme=mean(Cmeans(FVs.faces(elem,:)));
            K1e=mean(Lambda1s(FVs.faces(elem,:)));
            K2e=mean(Lambda2s(FVs.faces(elem,:)));
            if Kme<0 && Kge>0
                Kge=-1*Kge;
            end
            Kxyzs(elem,:)=centroid;
            Kgevs(elem)=Kge;
            Kmevs(elem)=Kme;
            K1evs(elem)=K1e;
            K2evs(elem)=K2e;
        end
        Kgev=[Kgevs;Kgevf];
        if fluids && solids
            Kxyz=[Kxyzs;Kxyzf];
            KInds=[ones(size(Kgevs));ones(size(Kgevf)).*2];
            Kmev=[Kmevs;Kmevf];
            K1ev=[K1evs;K1evf];
            K2ev=[K2evs;K2evf];
        elseif fluids
            Kxyz=[Kxyzf];
            KInds=[ones(size(Kgevf)).*2];
            Kmev=[Kmevf];
            K1ev=[K1evf];
            K2ev=[K2evf];
        elseif solids
            Kxyz=[Kxyzs];
            KInds=[ones(size(Kgevs))];
            Kmev=[Kmevs];
            K1ev=[K1evs];
            K2ev=[K2evs];
        end

        % remove outliers by zeroing/nan them
        if outlierFlag
            interQuartileRange=quantile(Kgev,0.75)-quantile(Kgev,0.25);
            minLim=quantile(Kgev,0.25)-1.5*interQuartileRange;
            maxLim=quantile(Kgev,0.75)+1.5*interQuartileRange;
            Kgev(Kgev>maxLim)=0;
            Kgev(Kgev<minLim)=0;

            interQuartileRange=quantile(Kmev,0.75)-quantile(Kmev,0.25);
            minLim=quantile(Kmev,0.25)-1.5*interQuartileRange;
            maxLim=quantile(Kmev,0.75)+1.5*interQuartileRange;
            Kmev(Kmev>maxLim)=nan;
            Kmev(Kmev<minLim)=nan;

            interQuartileRange=quantile(K1ev,0.75)-quantile(K1ev,0.25);
            minLim=quantile(K1ev,0.25)-1.5*interQuartileRange;
            maxLim=quantile(K1ev,0.75)+1.5*interQuartileRange;
            K1ev(K1ev>maxLim)=nan;
            K1ev(K1ev<minLim)=nan;

            interQuartileRange=quantile(K2ev,0.75)-quantile(K2ev,0.25);
            minLim=quantile(K2ev,0.25)-1.5*interQuartileRange;
            maxLim=quantile(K2ev,0.75)+1.5*interQuartileRange;
            K2ev(K2ev>maxLim)=nan;
            K2ev(K2ev<minLim)=nan;
        end
        KCells{1,1}=Kgev;
        KCells{1,2}=Kmev;
        KCells{1,3}=K1ev;
        KCells{1,4}=K2ev;
        KCells{1,5}=Kxyz;
        KCells{1,6}=KInds;
        
        if sunFlag
            Kgev(isnan(Kgev))=0;
            Kgev(Kmev<0 & Kgev>0)=-1.*Kgev(Kmev<0 & Kgev>0);
%             Kgev(Kmev>0 & Kgev<0)=-1.*Kgev(Kmev>0 & Kgev<0);

            thetaMacroSun=abs((4*pi*chi-sum(Kgev.*faceArea))/(4*nContactLines)*180/pi);
            thetaMacroSun=repmat(thetaMacroSun,[contactPerimeter,1]);
        end        
        
        if kearneyFlag
            Kgevf(isnan(Kgevf))=0;
            numFluidSurfaces = 0;
            % compute number of fluid-fluid interfaces
            [~, ~, rmElems]=removeUnconnectedTri(Ff,V);
            numFluidSurfaces = numFluidSurfaces + 1;
%             drawMesh(V,Fs(rmElems,:));view(3);drawnow
            Ft=Ff(~rmElems,:);
            faceAreat= 1/2*vecnorm(cross(V(Ft(:,2),:)-V(Ft(:,1),:),V(Ft(:,3),:)-V(Ft(:,1),:))')';
            while sum(faceAreat)./sum(faceAreaf)>1e-5 && size(Ft,1)>2
                [~, ~, rmElems]=removeUnconnectedTri(Ft,V);
%                 drawMesh(V,Ft(rmElems,:));view(3);drawnow
                Ft=Ft(~rmElems,:);
                numFluidSurfaces = numFluidSurfaces + 1;
                faceAreat= 1/2*vecnorm(cross(V(Ft(:,2),:)-V(Ft(:,1),:),V(Ft(:,3),:)-V(Ft(:,1),:))')';
            end            
            % compute line curvatures at contact line nodes
            Kvec=[];
            if plotFlag
                figure(4);clf(4);
                figure(5);clf(5);
            end
            for n=1:nContactLines
                edgeNodeInds=unique(Effs{n}(:),'stable');
                % edgeNodeInds=edgeNodeInds(1:end-1);
                X=V(edgeNodeInds,:);
                X=[X(end,:);X;X(1,:)];
                [LL,R,K] = curvature(X);
                X=X(2:end-1,:);
                LL=LL(2:end-1,:);
                R=R(2:end-1,:);
                K=K(2:end-1,:);
                R=filloutliers(R,'linear','quartiles');
                K=filloutliers(K,'linear','quartiles');
                Kvec=[Kvec;K];
                if plotFlag
                    figure(4);
                    plot(LL,1./R)
                    hold on
                    xlabel L
                    ylabel R
                    title('Curvature vs. cumulative curve length')
                    figure(5);
                    h = plot3(X(:,1),X(:,2),X(:,3),'-o'); 
                    grid on; 
                    axis equal
                    set(h,'marker','.');
                    xlabel x
                    ylabel y
                    zlabel z
                    title('Contact Line curvature vectors')
                    hold on
                    quiver3(X(:,1),X(:,2),X(:,3),K(:,1),K(:,2),K(:,3));
                    axis equal tight
                end
            end

            
           
            in(isnan(NsLineNodesVec))=0;
%             meanIn=~~round(mean(in));
            in=~~in;
            thetaMacroKearney=zeros(numel(in),1);
            phi=(acos(dot(NsLineNodesVec,Kvec,2))-pi/2)*180/pi;
            thetaMacroKearney=180-abs(phi+(acos(1-(4*pi-sum(Kgevf.*faceAreaf))/(2*pi*nContactLines)))*180/pi);
%             if sum(Kgevf.*faceAreaf)>0
%                 if sum(phi==real(phi))/numel(phi)>0.5 && quantile(thetaMacroKearney,0.5)>90
%                     thetaMacroKearney=180-thetaMacroKearney;
%                 end
%                 if sum(phi==real(phi))/numel(phi)<0.5 && quantile(thetaMacroKearney,0.5)<90
%                     thetaMacroKearney=180-thetaMacroKearney;
%                 end
%             else
%                 if sum(phi==real(phi))/numel(phi)>0.5 && quantile(thetaMacroKearney,0.5)<90
%                     thetaMacroKearney=180-thetaMacroKearney;
%                 end
%                 if sum(phi==real(phi))/numel(phi)<0.5 && quantile(thetaMacroKearney,0.5)>90
%                     thetaMacroKearney=180-thetaMacroKearney;
%                 end
%             end
%             if mean(in)<0.5 && quantile(thetaMacroKearney,0.5)<90
%                 thetaMacroKearney=180-thetaMacroKearney;
%             end
%             if mean(in)>0.5 && quantile(thetaMacroKearney,0.5)>90
%                 thetaMacroKearney=180-thetaMacroKearney;
%             end
%             thetaMacroKearney(in)=180-abs(phi(in)+(acos(1-(4*pi*numFluidSurfaces-sum(Kgevf.*faceAreaf))/(2*pi*nContactLines)))*180/pi);
%             thetaMacroKearney(~in)=abs(phi(~in)+(acos(1-(4*pi*numFluidSurfaces-sum(Kgevf.*faceAreaf))/(2*pi*nContactLines)))*180/pi);
%             if sunFlag
%                 if ~meanIn
%                     thetaMacroSun=180-thetaMacroSun;
%                 end
%             end
        end
    end
    
    if directFlag
        
        FVf=triangulation(Ff,V);
        Pf = incenter(FVf);
        Nf = faceNormal(FVf);  
        Nf=-1.*Nf;
        if plotFlag
            figure(6);clf(6);
            trisurf(FVf, ...
                 'FaceColor','cyan','FaceAlpha',0.8);
            axis equal
            hold on  
            quiver3(Pf(:,1),Pf(:,2),Pf(:,3), ...
                 Nf(:,1),Nf(:,2),Nf(:,3),5,'color','b');
            title('Normals at Fluid-Fluid Interface')
            axis equal tight
        end

        NfLineNodesVec=[];
        in=[];
        for n=1:nContactLines
            edgeNodeInds=unique(Effs{n}(:),'stable');
            X=V(edgeNodeInds,:);
            NfLineNodes=zeros(size(X));
            for node=1:numel(edgeNodeInds)
                attachedElems=vertexAttachments(FVf,edgeNodeInds(node));
                NfElems=Nf(attachedElems{:},:);
                PfElems=Pf(attachedElems{:},:);
                AfElems=faceAreaf(attachedElems{:});
                AfElems(AfElems<minTriangleArea)=0;
                AfElems=AfElems./meanMeshFaceArea;
                NfLineNodes(node,:)=NfElems'*AfElems;
                if plotFlag
                    trisurf(FVf(attachedElems{:},:),V(:,1),V(:,2),V(:,3),'FaceColor','r');
                    quiver3(PfElems(:,1),PfElems(:,2),PfElems(:,3), NfElems(:,1).*AfElems,NfElems(:,2).*AfElems,NfElems(:,3).*AfElems,5,'color','g');
                end
            end
            NfLineNodes=NfLineNodes./vecnorm(NfLineNodes,2,2);
            for dim=1:3
                NfLineNodes(:,dim)=filloutliers(NfLineNodes(:,dim),'linear','quartiles');
            end
            for dim=1:3
                NfLineNodes(:,dim)=fillmissing(NfLineNodes(:,dim),'movmean',10);
            end               

            NfLineNodesVec=[NfLineNodesVec;NfLineNodes];
            if plotFlag
                quiver3(X(:,1),X(:,2),X(:,3),NfLineNodes(:,1),NfLineNodes(:,2),NfLineNodes(:,3),'color','b');
            end         
            % find which surface normal vectors point into the cluster and
            % which point out of the cluster - travel a very small distance
            % along the vector and check if inside or outside F,V  
            X2 = X + (NfLineNodes).*1e-5;
%                 scatter3(X2(:,1),X2(:,2),X2(:,3),'r');
            try
                in = [in;intriangulation(V,F,X)];
            catch
                in=[];
            end
        end
        NfLineNodesVec(isnan(NfLineNodesVec))=0;

        directTheta=acos(dot(NfLineNodesVec,NsLineNodesVec,2))*180/pi;
        directTheta=180-real(directTheta);

%         % calculate tangents and projected tangents at contact line nodes
%         % for each line node, calculate centred secant line, and plane with that as
%         % normal. intersect plane and meshes, and determine angles up to x lines
%         % away
%         if plotFlag
%             figure(8);clf(8);
%             drawMesh(V,F);axis equal;axis tight;view(3);alpha(0.1)
%             title('Direct Measurement Selected Planes and Points')
%         end
%         directTheta=[];
%         e = meshEdges(F);
%         edges = [ V(e(:,1), :) V(e(:,2), :) ];
%         for n=1:nContactLines
%             edgeNodeInds=unique(Effs{n}(:),'stable');
%             % edgeNodeInds=edgeNodeInds(1:end-1);
%             X=V(edgeNodeInds,:);
%             X=[X(end,:);X;X(1,:)];
%             for node=2:size(edgeNodeInds)+1
%                 % find the tangent line, and normal plane
%                 NN=X(node+1,:)-X(node-1,:);
%                 PLANE = createPlane(X(node,:), NN);
%                 if any(isnan(PLANE))
%                     directTheta=[directTheta;nan];
%                     continue
%                 end
%         %         subplot(1,2,1)
%                 try
%                     polys = intersectPlaneMeshYDW(PLANE, V, F, e, edges);
%             %         polyf = intersectPlaneMesh(PLANE, V, Ff);
%                     if plotFlag
%     %                     if mod(node,10)==1
%                             drawPoint3d(X(node,:));hold on
%                             drawPlane3d(PLANE);alpha(0)
%                             drawPolygon3d(polys, 'LineWidth', 2);
%     %                     end
%                     end
%                     % sweep the normal plane to find the angle  - find when line passes at
%                     % least once, and up to max distance bwteen 2 points
%     %                         plot3(polys{1}(:,1),polys{1}(:,2),polys{1}(:,3),'-')
%             %         polys=cell2mat(polys');
%                     for p=1:numel(polys)
%                         polysTemp=polys{p};
%                         distsq=(X(node,1)-polysTemp(:,1)).^2+(X(node,2)-polysTemp(:,2)).^2+(X(node,3)-polysTemp(:,3)).^2;
%                         [minDist,i]=min(distsq);
%                         if minDist<1
%                             polys=polys{p};
%                             break
%                         end
%                     end
%             %         [p,i]=ismember(X(node,:),polys,'rows');
%                     angle=nan;
%                     oldAngle=0;
%                     inc=3;
%                     converged=0;
%                     while isnan(angle) || ~converged
%                         bef=i-inc;
%                         aft=i+inc;
%                         if bef<1
%                             bef=size(polys,1)+bef;
%                         elseif aft>size(polys,1)
%                             aft=aft-size(polys,1);
%                         end
%                 %         % use distance logic to avoid index overflow
%                 %         distsq=(polys(i,1)-polys(:,1)).^2+(polys(i,2)-polys(:,2)).^2+(polys(i,3)-polys(:,3)).^2;
%                 %         [~,srtinds]=sort(distsq,'ascend');
%                 %         bef=srtinds(2);
%                 %         aft=srtinds(3);
%                         P0=polys(bef,:);
%                         P1=polys(i,:);
%                         P2=polys(aft,:);
%                         v1=P0-P1;
%                         v2=P2-P1;
%                         v1=v1./norm(v1);
%                         v2=v2./norm(v2);
%                         res=dot(v1,v2);
%                         angle=acos(res)*180/pi;
%                         delta=(angle-90)*(oldAngle-90);
%                         inc=inc+1;
%                         if delta>0
%                             converged=1;
%             %                 angle=oldAngle;
%                         end
%                         oldAngle=angle;
%                     end
%                     if plotFlag
%                         c=jet(180);
%                         c=c(ceil(angle),:);
%                         plot3(P1(1),P1(2),P1(3),'.','MarkerSize',25,'Color',c)
%                         colormap jet
%                         c=colorbar;
%                         c.TickLabels={'0', '18', '36', '54', '72', '90', '108', '126', '144', '162', '180'};
%                     end
%                     directTheta=[directTheta;real(angle)];
%                 catch
%                     directTheta=[directTheta;nan];
% 
%                 end
%             end
%         end
        
    end
    
end