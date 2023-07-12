function thetaVec=computeMeshContactAnglesOnImage(fullDomain,geometricFlag,topologicalFlag,extendedTopologicalFlag,clusterInd,solidInd,minsize,maxsize,varargin)
    %%
    [Nx,Ny,Nz]=size(fullDomain);
 
    % get the domain and isolate the clusters of interest
    % get the oil blobs
    oil=fullDomain==clusterInd;
    
    CC=bwconncomp(oil);
    
    [~,I] = sort(cellfun(@length,CC.PixelIdxList));
    CC.PixelIdxList = CC.PixelIdxList(I);
    
    minSize=minsize^3;
    maxSize=maxsize^3;
    %
    thetaMeanVecSun=[];
    thetaMeanVecKearney=[];
    thetaMeanVecDirect=[];
    thetaMeanVecActual=[];
    thetaVec=cell(CC.NumObjects,1);
    tic;
    figure;

    for n=1:CC.NumObjects-1
        
        tempDomain=zeros(Nx,Ny,Nz);
        thetaVec{n}=struct();
        if numel(CC.PixelIdxList{n})>minSize && numel(CC.PixelIdxList{n})<maxSize
            % if blob within size range, isolate the blob with the solid
            [i,j,k]=ind2sub([Nx,Ny,Nz],CC.PixelIdxList{n});
            dL=1;
            DX=[max(min(i)-dL,1):min(max(i)+dL,Nx)];
            DY=[max(min(j)-dL,1):min(max(j)+dL,Ny)];
            DZ=[max(min(k)-dL,1):min(max(k)+dL,Nz)];
            locOrigin=[DX(1)-1,DY(1)-1,DZ(1)-1];
            tempDomain(CC.PixelIdxList{n})=1;
            tempDomainBB=tempDomain(DX,DY,DZ);
            domain=fullDomain(DX,DY,DZ);
            domain(domain~=solidInd)=2;
            domain(domain==solidInd)=0;
            domain(tempDomainBB==1)=1;
            % remove falsely generated isolated oil drops, keep the coalesed
            % ones
            CCw=bwconncomp(domain==2,6);
            for nn=1:CCw.NumObjects
                if numel(CCw.PixelIdxList{nn})<5
                    domain(CCw.PixelIdxList{nn})=1;
                end
            end

            CCw=bwconncomp(domain==1,6);
            for nn=1:CCw.NumObjects
                if numel(CCw.PixelIdxList{nn})<5
                    domain(CCw.PixelIdxList{nn})=2;
                end
            end
            thetaActual=0;
            if numel(varargin)>0
                bbField=varargin{1}(DX,DY,DZ);
                
                dClusters=bwdist(domain==1);
                dOtherPhase=bwdist(domain==2);
                dInternalSolid=bwdist(~(domain==0));
                % the contact line voxels within the solid are 
                solidContactLineVoxels=dInternalSolid<1.8 & dInternalSolid>0 & dClusters<2.8 & dClusters>0 & dOtherPhase<1.8 & dOtherPhase>0;
                thetaActual=bbField(solidContactLineVoxels);
            end
            
            if sum(domain(:)==0)>4 && sum(domain(:)==2)>4
%                 disp('Generating Mesh');
                % generate the mesh for this blob
                minTriangles=1000;
                smoothFlag=0;
                smoothingOptVec=[100,3,1000];
                sunFlag=0;
                kearneyFlag=0;
                directFlag=0;
                curvaturesFlag=0;
                outlierFlag=0;
                fluids=0;
                solids=0;
                plotFlag=0;
                [thetaMacroSun,thetaMacroKearney,directTheta,V,Ff,Fs,Effs,voxelInds,KCells,numUpsamples] = ...
                    computeSmoothContactAnglesOnMeshv4(...
                    domain,minTriangles,0,smoothFlag,smoothingOptVec,sunFlag,kearneyFlag,directFlag,curvaturesFlag,outlierFlag,fluids,solids,plotFlag);
                if numel(Fs)<10 || numel(Ff)<10 || numel(Effs)<1
%                     disp('Cluster not in sufficient contact with surface/fluid')
                    continue
                end
                thetaVec2=[];
                F=[Fs;Ff];

                vol=meshVolume(V,F);
                volFactor=1;
%                 disp(['volFactor: ', num2str(volFactor),', Sun: ', num2str(180-mean(thetaMacroSun)), ', Kearney: ', num2str(quantile(thetaMacroKearney,0.5)), ' Actual: ', num2str(mean(thetaActual))])

                i=1;
                % smooth the mesh
%                 disp('Smoothing Mesh - tuning smoothing to match kearney and actual');
                while volFactor>0.9% && abs(mean(thetaActual)-quantile(thetaMacroKearney,0.5))>10
                    F=[Fs;Ff];
                    domain=struct();
                    domain.V=V;
                    domain.Fs=Fs;
                    domain.Ff=Ff;
                    domain.F=F;
                    domain.voxelInds=voxelInds;
                    domain.Effs=Effs;

                    minTriangles=1;
                    smoothFlag=1;
                %     if i<L/2
%                     if volFactor>0.95
%                         smoothingOptVec=[2,1,2];
%                     else
                        smoothingOptVec=[2,2,2];
%                     end
                %     else
                %         smoothingOptVec=[2,2,0];
                %     end
                
                    sunFlag=0;
                    kearneyFlag=0;
                    
                    directFlag=0;
                    curvaturesFlag=0;
                    outlierFlag=0;
                    fluids=0;
                    solids=0;
                    plotFlag=0;
                    [thetaMacroSun,thetaMacroKearney,directTheta,V,Ff,Fs,Effs,voxelInds,KCells,numUpsamples] = ...
                        computeSmoothContactAnglesOnMeshv4(...
                        domain,minTriangles,numUpsamples,smoothFlag,smoothingOptVec,sunFlag,kearneyFlag,directFlag,curvaturesFlag,outlierFlag,fluids,solids,plotFlag);
                    volFactor=meshVolume(V,F)./vol;
%                 	disp(['volFactor: ', num2str(volFactor),', Sun: ', num2str(180-mean(thetaMacroSun)), ', Kearney: ', num2str(quantile(thetaMacroKearney,0.5)), ' Actual: ', num2str(mean(thetaActual))])

    %                 thetaVec2=[thetaVec2;i,volFactor,mean(thetaMacroSun),mean(thetaMacroKearney),mean(directTheta)];
                    i=i+1;
                end
%                 % calculate the final contact angle
%                 disp('Calculating Angles');

                F=[Fs;Ff];
                domain=struct();
                domain.V=V;
                domain.Fs=Fs;
                domain.Ff=Ff;
                domain.F=F;
                domain.voxelInds=voxelInds;
                domain.Effs=Effs;

                minTriangles=1;
                smoothFlag=0;
            %     if i<L/2
                    smoothingOptVec=[4,4,2];
            %     else
            %         smoothingOptVec=[2,2,0];
            %     end
                sunFlag=topologicalFlag;
                kearneyFlag=extendedTopologicalFlag;
                directFlag=geometricFlag;
                curvaturesFlag=1;
                outlierFlag=1;
                fluids=1;
                solids=1;
                plotFlag=0;
                [thetaMacroSun,thetaMacroKearney,directTheta,V,Ff,Fs,Effs,voxelInds,KCells,numUpsamples] = ...
                    computeSmoothContactAnglesOnMeshv4(...
                    domain,minTriangles,numUpsamples,smoothFlag,smoothingOptVec,sunFlag,kearneyFlag,directFlag,curvaturesFlag,outlierFlag,fluids,solids,plotFlag);
                % confirm that the number of Effs equals the direct and
                % keraney calculations
                numEdgeNodes=sum(cellfun(@numel,voxelInds));
                if geometricFlag
                    assert(numEdgeNodes==numel(directTheta))
                end
                if extendedTopologicalFlag
                    assert(numEdgeNodes==numel(thetaMacroKearney))
                end
                % no clear way to fix this stupid bullshit for now
%                 if  mean(thetaMacroKearney)<90 
%                     thetaMacroSun=180-thetaMacroSun;
%                     thetaMacroKearney=180-thetaMacroKearney;
% 
%                     %                 elseif mean(thetaMacroSun)>90 && mean(thetaMacroKearney)<90 
% %                     thetaMacroSun=180-thetaMacroSun;
%                 end
                if ~isnan(mean(directTheta)) && ~isempty([mean(directTheta),mean(thetaMacroSun),mean(thetaMacroKearney)])
                    
%                     if abs(quantile(thetaMacroKearney,0.5)-mean(thetaActual))>abs(180-quantile(thetaMacroKearney,0.5)-mean(thetaActual))
%                         thetaMacroKearney=180-thetaMacroKearney;
%                         disp(['Acute-Obtuse Error'])

%                     end
                    thetaVec{n}.thetaMacroSun=180-thetaMacroSun;
                    thetaVec{n}.thetaMacroKearney=thetaMacroKearney;
                    thetaVec{n}.directTheta=directTheta;
                    thetaVec{n}.V=V;
                    thetaVec{n}.Ff=Ff;
                    thetaVec{n}.Fs=Fs;
                    thetaVec{n}.Effs=Effs;
                    thetaVec{n}.KCells=KCells;
                    thetaVec{n}.numUpsamples=numUpsamples;
                    thetaVec{n}.domainInds=CC.PixelIdxList{n};
                    thetaVec{n}.voxelInds=voxelInds;
% 
                    thetaMeanVecSun=[thetaMeanVecSun;180-thetaMacroSun];
                    thetaMeanVecKearney=[thetaMeanVecKearney;thetaMacroKearney];
                    thetaMeanVecDirect=[thetaMeanVecDirect;directTheta];
                    disp([])

                    disp(['Time: ',num2str(toc),'s, Computing Blob ',num2str(n),' of ', num2str(CC.NumObjects), '. EquivRadius: ', num2str((numel(CC.PixelIdxList{n}).^0.33)),', Direct: ', num2str(mean(directTheta)), ', Kearney: ', num2str(quantile(thetaMacroKearney,0.5)), ', Sun: ', num2str(mean(thetaMacroSun)), ', Actual: ', num2str(mean(thetaActual))])                    
                    clf;
                    view(2);xlim([0 180]);hold on
                    histogram(thetaMeanVecSun,'binwidth',7,'normalization','pdf','FaceColor','none','edgecolor','r');hold on
                    histogram(thetaMeanVecKearney,'binwidth',6,'normalization','pdf','FaceColor','none','edgecolor','g');hold on
                    histogram(thetaMeanVecDirect,'binwidth',5,'normalization','pdf','FaceColor','none','edgecolor','b');hold off
                    legend('Sun', 'Kearney', 'Direct')
                    drawnow;
                    
                    
%                     thetaMeanVecSun=[thetaMeanVecSun;180-thetaMacroSun(1)];
%                     thetaMeanVecActual=[thetaMeanVecActual;mean(thetaActual)];
%                     thetaMeanVecKearney=[thetaMeanVecKearney;quantile(thetaMacroKearney,0.5)];
%                     thetaMeanVecDirect=[thetaMeanVecDirect;mean(directTheta)];
%                     figure(2);
%                     clf;
% %                     radiusVsTheta=[radiusVsTheta;numel(cluster.domainInds).^(1/3),mean(cluster.thetaMacroSun),mean(cluster.thetaMacroKearney),mean(cluster.directTheta)];
%                     subplot 131
%                     plot(thetaMeanVecActual,[thetaMeanVecDirect],'o')
%                     xlim([0 180])
%                     ylim([0 180])
%                     hold on
%                     plot([0 180],[0 180],'-k')
%                     
%                     subplot 132
%                     plot(thetaMeanVecActual,[thetaMeanVecSun],'o')
%                     xlim([0 180])
%                     ylim([0 180])
%                     hold on
%                     plot([0 180],[0 180],'-k')
%                     
%                     subplot 133
%                     plot(thetaMeanVecActual,[thetaMeanVecKearney],'o')
%                     xlim([0 180])
%                     ylim([0 180])
%                     hold on
%                     plot([0 180],[0 180],'-k')
%                     drawnow
                end
                
            end
        end
    end


end
