addpath(genpath(pwd))

resamplingFactor=10;

%% load the contact angle data
thetaData0=load('thetaData0.mat');
thetaData0=thetaData0.thetaData0;
%% load in the geoetry as well

solidInd=1;

rockImg0=readTiffStackYDW('segmentedRockTestTheta30.tif');
solid=uint8(rockImg0==solidInd);
%% begin remapping - map clusters and contact loops

domainSun=solid;
domainKearney=solid;
domainDirect=solid;
domainSunLines=solid;
domainKearneyLines=solid;
domainDirectLines=solid;
thetaVec=[];
for i=numel(thetaData0)-1:-1:1
    cluster=thetaData0{i};
    if ~isempty(fieldnames(cluster))        
        [x,y,z]=ind2sub(size(solid),cluster.domainInds);

        domainSun(cluster.domainInds)=uint8(mean(cluster.thetaMacroSun));
        domainKearney(cluster.domainInds)=uint8(quantile(cluster.thetaMacroKearney,0.5));
        
        x=max([min(x)-3,1]):min([max(x)+3,size(solid,1)]);
        y=max([min(y)-3,1]):min([max(y)+3,size(solid,2)]);
        z=max([min(z)-3,1]):min([max(z)+3,size(solid,3)]);
        
        bbDomain=domainSun(x,y,z);
        bbSolid=solid(x,y,z);
        domainSunClusters=bbDomain>1;
        dClusters=bwdist(domainSunClusters);
        dOtherPhase=bwdist(bbDomain==0);
        dInternalSolid=bwdist(~(bbDomain==1));
        solidContactLineVoxels=dInternalSolid<1.8 & dInternalSolid>0 & dClusters<2.8 & dClusters>0 & dOtherPhase<1.8 & dOtherPhase>0;
        if ~isempty(cluster.directTheta)
            domainDirect(cluster.domainInds)=uint8(mean(cluster.directTheta));
            bbSolid(solidContactLineVoxels)=uint8(mean(cluster.directTheta));
            domainDirectLines(x,y,z)=bbSolid;
            
            
            bbSolid(solidContactLineVoxels)=uint8(mean(cluster.thetaMacroSun));
            domainSunLines(x,y,z)=bbSolid;

            if abs(quantile(cluster.thetaMacroKearney,0.5)-mean(cluster.directTheta))>abs((180-quantile(cluster.thetaMacroKearney,0.5))-mean(cluster.directTheta))
                cluster.thetaMacroKearney=180-cluster.thetaMacroKearney;
            end
            
            bbSolid(solidContactLineVoxels)=uint8(quantile(cluster.thetaMacroKearney,0.5));
            domainKearneyLines(x,y,z)=bbSolid;
            disp([num2str(i),' of ', num2str(numel(thetaData0)),', Direct: ', num2str(mean(cluster.directTheta)), ', Kearney: ', num2str(quantile(cluster.thetaMacroKearney,0.5)), ' Sun: ', num2str(mean(cluster.thetaMacroSun))])
        end

    end
end
writeTiffStackYDW('./domainSun.tif',domainSun);
writeTiffStackYDW('./domainSunLines.tif',domainSunLines);
writeTiffStackYDW('./domainKearney.tif',domainKearney);
writeTiffStackYDW('./domainKearneyLines.tif',domainKearneyLines);
writeTiffStackYDW('./domainDirect.tif',domainDirect);
writeTiffStackYDW('./domainDirectLines.tif',domainDirectLines);
%% use the loop data and interpolate over the solid


domainSunLinesAll=readTiffStackYDW('./domainSunLines.tif');
domainKearneyLinesAll=readTiffStackYDW('./domainKearneyLines.tif');
domainDirectLinesAll=readTiffStackYDW('./domainDirectLines.tif');

[domainSunInterp] = interpolateContactLoopstoSolid(domainSunLinesAll,resamplingFactor);
[domainKearneyInterp] = interpolateContactLoopstoSolid(domainKearneyLinesAll,resamplingFactor);
[domainDirectInterp] = interpolateContactLoopstoSolid(domainDirectLinesAll,resamplingFactor);

writeTiffStackYDW('./domainKearneyInterp.tif',domainKearneyInterp);
writeTiffStackYDW('./domainSunInterp.tif',domainSunInterp);
writeTiffStackYDW('./domainDirectInterp.tif',domainDirectInterp);
