addpath(genpath(pwd))
rockImg0=readTiffStackYDW('segmentedRockTestTheta30.tif');

%% for each rock image, calculate the contact angles using the 10% shinkrage cutoff on blobs larger than 8-cubed
warning off
geometricFlag=1;
topologicalFlag=0;
extendedTopologicalFlag=1;
clusterInd=2;
solidInd=1;
minRad=8;
maxRad=64;


thetaData0=computeMeshContactAnglesOnImage(rockImg0,geometricFlag,topologicalFlag,extendedTopologicalFlag,clusterInd,solidInd,minRad,maxRad);

%%
names=fieldnames(thetaDataStruct);
angles=0;
domain={'Sandstone'};
allClustersPlot=figure('Name','3','units','normalized','outerposition',[0 0 0.4 0.5]);

i=1;


    data=thetaData0;
    thetaMeanVecSun=[];
    thetaMeanVecKearney=[];
    thetaMeanVecDirect=[];
    for n=1:numel(data)-1

        clusterData=data{n};
        if ~isempty(fieldnames(clusterData))


            thetaMacroSun=clusterData.thetaMacroSun;
            thetaMacroKearney=clusterData.thetaMacroKearney;
            directTheta=clusterData.directTheta;
            V=clusterData.V;
            Ff=clusterData.Ff;
            Fs=clusterData.Fs;
            Effs=clusterData.Effs;
            KCells=clusterData.KCells;
            numUpsamples=clusterData.numUpsamples;
            PixelIdxList=clusterData.domainInds;
            disp(['Loaded Blob ',num2str(n),' of ', num2str(numel(data)), '. EquivRadius: ', num2str((numel(PixelIdxList).^0.33))])

            thetaMeanVecSun=[thetaMeanVecSun;thetaMacroSun];
    %                 thetaMeanVecKearney=[thetaMeanVecKearney;thetaMacroKearney./thetaMacroKearney.*quantile(thetaMacroKearney,0.5)];
            thetaMeanVecKearney=[thetaMeanVecKearney;thetaMacroKearney];
            thetaMeanVecDirect=[thetaMeanVecDirect;directTheta];

        end

    end
%     clf;
    view(2);xlim([0 180]);hold on
    histogram(180-thetaMeanVecSun,'binwidth',9,'normalization','pdf','FaceColor','none','edgecolor',[0 0.4470 0.7410],'linewidth',2,'displaystyle','stairs');hold on
    histogram(180-thetaMeanVecKearney,'binwidth',10,'normalization','pdf','FaceColor','none','edgecolor',[0.8500 0.3250 0.0980],'linewidth',2,'displaystyle','stairs');hold on
    h=histogram(180-thetaMeanVecDirect,'binwidth',11,'normalization','pdf','FaceColor','none','edgecolor',[0.9290 0.6940 0.1250],'linewidth',2,'displaystyle','stairs');hold on
    plot([angles(i), angles(i)], [0 h.Parent.YLim(2)],'color',[0.6350 0.0780 0.1840],'linewidth',2)
    legend('ThetaMacro', 'ThetaMacroExtended', 'Direct', 'Surface Affinity')


    xlabel('Calculated Contact Angle (degrees)')
    ylabel('PDF')
    title([domain{i}, ', Surface Affinity = ',num2str(angles(i)),' (deg)'])

suptitle('Sandstone and Gas Diffusion Layer Calculated Wettabilities (All Clusters)')

% 
% print(allClustersPlot,'/home/user/Insync/thetaFigures/allClustersPlot','-dpng', '-r300')
% img3=imread('/home/user/Insync/thetaFigures/allClustersPlot.png');
% img3=autoCrop(img3);
% imwrite(img3, '/home/user/Insync/thetaFigures/allClustersPlot.png');
