
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
% ds=bwdist(solid);
thetaVec=[];
for i=numel(thetaData0)-1:-1:1
    cluster=thetaData0{i};
    if ~isempty(fieldnames(cluster))        
        [x,y,z]=ind2sub(size(solid),cluster.domainInds);

        domainSun(cluster.domainInds)=uint8(mean(cluster.thetaMacroSun));
        domainKearney(cluster.domainInds)=uint8(quantile(cluster.thetaMacroKearney,0.5));
        
        % identify the contact lines voxels and map value
        x=max([min(x)-3,1]):min([max(x)+3,size(solid,1)]);
        y=max([min(y)-3,1]):min([max(y)+3,size(solid,2)]);
        z=max([min(z)-3,1]):min([max(z)+3,size(solid,3)]);
        
        bbDomain=domainSun(x,y,z);
        bbSolid=solid(x,y,z);
        domainSunClusters=bbDomain>1;
        dClusters=bwdist(domainSunClusters);
        dOtherPhase=bwdist(bbDomain==0);
        dInternalSolid=bwdist(~(bbDomain==1));
        % the contact line voxels within the solid are 
        solidContactLineVoxels=dInternalSolid<1.8 & dInternalSolid>0 & dClusters<2.8 & dClusters>0 & dOtherPhase<1.8 & dOtherPhase>0;
        if ~isempty(cluster.directTheta)
            domainDirect(cluster.domainInds)=uint8(mean(cluster.directTheta));
            bbSolid(solidContactLineVoxels)=uint8(mean(cluster.directTheta));
            domainDirectLines(x,y,z)=bbSolid;
            
            
            bbSolid(solidContactLineVoxels)=uint8(mean(cluster.thetaMacroSun));
            domainSunLines(x,y,z)=bbSolid;
%             
            if abs(quantile(cluster.thetaMacroKearney,0.5)-mean(cluster.directTheta))>abs((180-quantile(cluster.thetaMacroKearney,0.5))-mean(cluster.directTheta))
                cluster.thetaMacroKearney=180-cluster.thetaMacroKearney;
            end
            
            bbSolid(solidContactLineVoxels)=uint8(quantile(cluster.thetaMacroKearney,0.5));
            domainKearneyLines(x,y,z)=bbSolid;
%             thetaVec=[thetaVec;mean(thetaActual),mean(cluster.directTheta),quantile(cluster.thetaMacroKearney,0.5),mean(cluster.thetaMacroSun)];

%         else
%             bbSolid(solidContactLineVoxels)=uint8(mean(cluster.thetaMacroSun));
%             domainSunLines(x,y,z)=bbSolid;
%             bbSolid(solidContactLineVoxels)=uint8(quantile(cluster.thetaMacroKearney,0.5));
%             domainKearneyLines(x,y,z)=bbSolid;
%             thetaVec=[thetaVec;mean(thetaActual),mean(cluster.directTheta),quantile(cluster.thetaMacroKearney,0.5),mean(cluster.thetaMacroSun)];
            disp([num2str(i),' of ', num2str(numel(thetaData0)),', Direct: ', num2str(mean(cluster.directTheta)), ', Kearney: ', num2str(quantile(cluster.thetaMacroKearney,0.5)), ' Sun: ', num2str(mean(cluster.thetaMacroSun))])
%             sum(solidContactLineVoxels(:))
        end
%         if ~isempty(cluster.directTheta)
%             disp([num2str(i),' of ', num2str(numel(thetaData0)),' Actual: ', num2str(mean(thetaActual)), ', Direct: ', num2str(mean(cluster.directTheta)), ', Kearney: ', num2str(quantile(cluster.thetaMacroKearney,0.5)), ' Sun: ', num2str(mean(cluster.thetaMacroSun))])
%             thetaVec=[thetaVec;mean(thetaActual),mean(cluster.directTheta),quantile(cluster.thetaMacroKearney,0.5),mean(cluster.thetaMacroSun)];
%         end
%         radiusVsTheta=[radiusVsTheta;numel(cluster.domainInds).^(1/3),mean(cluster.thetaMacroSun),mean(cluster.thetaMacroKearney),mean(cluster.directTheta)];
%         figure(1);
%         plot(thetaVec(:,1),thetaVec(:,2:4),'o')
%         xlim([0 180])
%         ylim([0 180])
% %         hold on
%         drawnow
%         step=100;
%         scatter3(x(1:step:end),y(1:step:end),z(1:step:end))
%         drawnow
%         hold on
    end
end
writeTiffStackYDW('./domainSun.tif',domainSun);
writeTiffStackYDW('./domainSunLines.tif',domainSunLines);
writeTiffStackYDW('./domainKearney.tif',domainKearney);
writeTiffStackYDW('./domainKearneyLines.tif',domainKearneyLines);
writeTiffStackYDW('./domainDirect.tif',domainDirect);
writeTiffStackYDW('./domainDirectLines.tif',domainDirectLines);
%% use the loop data and interpolate over the solid

resamplingFactor=10;

domainSunLinesAll=readTiffStackYDW('./domainSunLines.tif');
domainKearneyLinesAll=readTiffStackYDW('./domainKearneyLines.tif');
domainDirectLinesAll=readTiffStackYDW('./domainDirectLines.tif');

data=imresize3(domainKearneyLinesAll,1/resamplingFactor,'nearest');
[nzX, nzY, nzZ] = ind2sub(size(data), find(data > 1));

% get x, y, z indices of zero elements
[zX, zY, zZ] = ind2sub(size(data), find(data <= 1));

% build vectors of x, y, z coordinates and corresponding values
nonzero_points = [nzX, nzY, nzZ];
nonzero_values = data(data >1);

% create scattered interpolant with inverse distance weighting
F = scatteredInterpolant(nonzero_points, double(nonzero_values), 'natural');

% now interpolate the zero values
zero_points = [zX, zY, zZ];
interpolated_values=zeros(size(zero_points,1),1);
for i=1:size(zero_points,1)
    if mod(i,1000)==0
        disp(['Interpolating Points: ', num2str(i),' of ',num2str(size(zero_points,1))])
    end
    interpolated_values(i) = F(zero_points(i,:));
end

% replace zero values in the data
data(data <= 1) = interpolated_values;
data=imresize3(data,resamplingFactor);
data(~solid)=0;

writeTiffStackYDW('./domainKearneyInterp.tif',data);


data=imresize3(domainSunLinesAll,1/resamplingFactor,'nearest');
[nzX, nzY, nzZ] = ind2sub(size(data), find(data > 1));

% get x, y, z indices of zero elements
[zX, zY, zZ] = ind2sub(size(data), find(data <= 1));

% build vectors of x, y, z coordinates and corresponding values
nonzero_points = [nzX, nzY, nzZ];
nonzero_values = data(data >1);

% create scattered interpolant with inverse distance weighting
F = scatteredInterpolant(nonzero_points, double(nonzero_values), 'natural');

% now interpolate the zero values
zero_points = [zX, zY, zZ];
interpolated_values=zeros(size(zero_points,1),1);
for i=1:size(zero_points,1)
    if mod(i,1000)==0
        disp(['Interpolating Points: ', num2str(i),' of ',num2str(size(zero_points,1))])
    end
    interpolated_values(i) = F(zero_points(i,:));
end

% replace zero values in the data
data(data <= 1) = interpolated_values;
data=imresize3(data,resamplingFactor);
data(~solid)=0;

writeTiffStackYDW('./domainSunInterp.tif',data);




data=imresize3(domainDirectLinesAll,1/resamplingFactor,'nearest');
[nzX, nzY, nzZ] = ind2sub(size(data), find(data > 1));

% get x, y, z indices of zero elements
[zX, zY, zZ] = ind2sub(size(data), find(data <= 1));

% build vectors of x, y, z coordinates and corresponding values
nonzero_points = [nzX, nzY, nzZ];
nonzero_values = data(data >1);

% create scattered interpolant with inverse distance weighting
F = scatteredInterpolant(nonzero_points, double(nonzero_values), 'natural');

% now interpolate the zero values
zero_points = [zX, zY, zZ];
interpolated_values=zeros(size(zero_points,1),1);
for i=1:size(zero_points,1)
    if mod(i,1000)==0
        disp(['Interpolating Points: ', num2str(i),' of ',num2str(size(zero_points,1))])
    end
    interpolated_values(i) = F(zero_points(i,:));
end

% replace zero values in the data
data(data <= 1) = interpolated_values;
data=imresize3(data,resamplingFactor);
data(~solid)=0;

writeTiffStackYDW('./domainDirectInterp.tif',data);
