% fig3
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% Uses data saved by saveTaskData.m
%
%   usage: fig3()
%   by: zvi roth
%   date: 9/3/2020
%   purpose: create figure 3, task-related response at 3 different
%   eccentricity bins, for one subject
%

function[] = fig3()

dataFolder = '/Volumes/TaskDrive/min/';%folder containing all fMRI data

onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;

subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

randSeed = rng;
load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
TR=1.5;
nperms=10000;
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
end

ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString  ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT');
fmriUnits = '% change image intensity';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
linewidth = 1;
markersize=10;
fontsize = 9;
ROIs = 1:length(roiNames);
ROIs = [1 2];
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22

ntrials=15;


%% bin voxels within each ROI according to eccentricity, then average within bins
%eccentricity bin borders
binBorders = [0 1 10 100];
nbins = length(binBorders)-1;
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            if onlyCorrect ==1 %ONLY CORRECT
                goodTrials = trialCorrectnessVec==1;
            elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
            else % including all trials with a response
                goodTrials = trialResponseVec>0;
            end
            
            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshapedTrials;

                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                subBinStd(iSub,iRoi,ibin,rwd,:) = std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2);
            end           
        end
    end
end

%% PLOT timecourse and FFT for s0008
figure(1); clf
iSub=2;
cols=nbins;
rows=1;
iRoi=2;
%timecourse
for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        dsErrorsurface(TR*(0:trialLength-1), squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)), squeeze(subBinStd(iSub,iRoi,ibin,rwd,:))./sqrt(size(subBinTrialResponse{iSub,iRoi,ibin,rwd},2)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1), squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)),'-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    xlabel('time (sec)');
    t = {'fMRI response'; ['(' fmriUnits ')'] };
    if ibin==1
        ylabel(t);
    end
end

%%
subplot(rows,cols,3)
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 8/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',fontsize);
    legend off
    axis square
end
set(gcf,'position',[10 10 20.5 	14]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_' ConcatProjStr '.pdf']);


%%