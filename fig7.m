% fig7
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% Uses data saved by saveTaskData.m
%
%   usage: fig7()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: create figure 7, task-related timepoint variability

function[] = fig7()

onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%1=data is z-scored, 0=no z-scoring
ConcatProj = 1;%1=project out global mean. 0=don't project out mean
curFolder = pwd;
TR=1.5;
dataFolder = '';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
 dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
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
    'subRoiRuns','runMeanFFT');

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
fontsize=9;
ROIs = [1:4];

%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22


for iSub = 1:length(goodSubs)
     for iRoi=1:length(roiNames)
        for rwd=1:2
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
        %trial-by-trial variability per subject
           subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
        end
     end
end

%%
iRoi=2;%right EVC

groupMeanVar = squeeze(mean(subTimepointStd(:,iRoi,:,:)));
groupStdVar = squeeze(std(subTimepointStd(:,iRoi,:,:)));

subVarDiff = squeeze(subTimepointStd(:,iRoi,2,:) - subTimepointStd(:,iRoi,1,:));
meanVarDiff = mean(subVarDiff);
stdVarDiff = std(subVarDiff);
figure(1); clf
subplots = {1:3, 4:6, 7 , 9};
rows=2;
cols = 9;
subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface(TR*(0:trialLength-1), groupMeanVar(rwd,:), groupStdVar(rwd,:)./sqrt(size(subTimepointStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(TR*(0:trialLength-1),groupMeanVar(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
end
subplot(rows,cols,subplots{2})
dsErrorsurface(TR*(0:trialLength-1), meanVarDiff, stdVarDiff./sqrt(size(subTimepointStd,1)), [0 0 0],dsSurfaceAlpha);
hold on
plot(TR*(0:trialLength-1),meanVarDiff,'k-','linewidth',linewidth,'markersize',markersize);
hline(0);

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot});
    xlabel('time (sec)');
    if isubplot==1
        ylabel('timepoint variability (std)');
    else
        ylabel('\Delta timepoint variability (std)');
    end
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -18/64, 'xAxisMargin', 6/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',fontsize);
    axis square
    legend off
end
% set(gcf,'position',[10 10 25 	6]);
set(gcf,'position',[10 10 25 	14]);
print('-painters','-dpdf',['fig7.pdf']);
