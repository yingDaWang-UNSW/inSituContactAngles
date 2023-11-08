addpath(genpath(pwd))
rockImg0=readTiffStackYDW('segmentedRockTestTheta30.tif');

%% for each rock image, calculate the contact angles using the 10% shinkrage cutoff on blobs larger than 8-cubed
warning off
geometricFlag=1;
topologicalFlag=1;
extendedTopologicalFlag=1;
clusterInd=2;
solidInd=1;
minRad=8;
maxRad=64;


thetaData0=computeMeshContactAnglesOnImage(rockImg0,geometricFlag,topologicalFlag,extendedTopologicalFlag,clusterInd,solidInd,minRad,maxRad);

save('thetaData0.mat','thetaData0','-v7.3')

