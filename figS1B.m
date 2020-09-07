% figS1B
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
%
%   usage: figS1B()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: create Figure S1B, task-related response amplitude for control
%   experiment, with no stimulus
%
%

function[] = figS1B()

toZscore=1;%0 or 1
ConcatProj = 0;
curFolder = pwd;

dataFolder = '';
% dataFolder = '~/das/';

zScoreString = '';
if toZscore
    zScoreString = '_zscored';
end

ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'das_'  zScoreString  ConcatProjStr  '.mat'],   'subResponse', 'roiMeanTseries', ...
    'roiTC', ...
    'subFolders', 'roiNames','subTrialResponse',...
    'eccen','ang','areas','trialLength',...
    'subStd',...
    'voxTrials','meanVoxTrial',...
    'voxNum','concatInfo');

TR=1.5;
% trialLength=10;
plotColors = { [0 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
fontsize=9;
ROIs = 1:length(roiNames);
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
ntrials=15;
goodSubs = [1:length(subFolders)];
%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs

            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi}.tSeries(binVoxels,:),1);%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin},trialLength,[]);

                %average per subject
                subBinResponse(iSub,iRoi,ibin,:) = mean(subBinTrialResponse{iSub,iRoi,ibin},2);
                %trial-by-trial variability per subject
                subBinVar(iSub,iRoi,ibin) = mean(std(subBinTrialResponse{iSub,iRoi,ibin},0,2));
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin});
                subBinStdVar(iSub,iRoi,ibin) = std(singleTrialStd);
                subBinStdMean(iSub,iRoi,ibin) = mean(singleTrialStd);
            end
    end
end

subBinAmp = std(subBinResponse,0,4);%get mean timecourse amplitude per subject.(iRoi,ibin,rwd,1)
binMeanAmp = squeeze(mean(subBinAmp,1)); %mean over subjects. (iRoi,ibin,rwd)
binStdAmp = squeeze(std(subBinAmp,0,1)); %(iRoi,ibin,rwd)

binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)

%%
i=0;

rows=2;
cols = 9;
subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};
figure(1); clf
for r= 1:length(ROIs)
    %amplitude
    subplot(rows,cols,subplots{r});
    iRoi=ROIs(r);
    dsErrorsurface(binCenters, squeeze(binMeanAmp(iRoi,:)), squeeze(binStdAmp(iRoi,:))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{1},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binMeanAmp(iRoi,:)),'.','Color', plotColors{1},'linewidth',linewidth,'markersize',markersize);
end

%% Formatting
for r = 1:length(ROIs)
    for isubplot=[1:2]
        subplot(rows,cols,subplots{isubplot});
        switch isubplot
            case 1
                ylabel('response amplitude (std)');
            case 2
                ylabel('response amplitude (std)');

        end
        xticks([0:20:60]);
   
        xlabel('eccentricity (deg)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -16/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
            'labelFontSize',fontsize);
        legend off
        axis square
        
    end
    set(gcf,'position',[10 10 28 	14]);
end

for iSub=1:length(goodSubs)
    diff(concatInfo{goodSubs(iSub)}.runTransition)
end

print('-painters','-dpdf',['figS1B.pdf']);
