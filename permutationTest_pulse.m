% permutationTest_pulse
%
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020).
% Task-related activity in human visual cortex.
%
%   uses file created by saveTaskData_localizer.m
%
%   usage: permutationTest_pulse()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: compute task-related response amplitude, timepoint variability, amplitude variability,
%   and temporal variability, and use permutations to test for significant
%   difference between high and low reward.
%
%

function[] = permutationTest_pulse()

dataFolder = '';

onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%1=data is z-scored, 0=no z-scoring
ConcatProj = 1;%1=project out global mean. 0=don't project out mean

load(fullfile(dataFolder,'physioRandSeed.mat'),'physioRandSeed');
rng(physioRandSeed);

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
load([dataFolder 'rwdTC_physio' onlyCorrectString zScoreString  ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'ecg','ecgPulseRate','interpPulseRate',...
    'respselect','resp',...
    'rwdPulseTC','rwdRvTC',...
    'designMatPulse','designMatRespPulse','designMatResp','deconvLength',...
    'allGoodTrials');

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);

goodSubs = [1:length(subFolders)]; 
% goodSubs = [2:length(subFolders)]; 

rwdString = {'H','L'};
ntrials=15;


%% PERMUTATIONS
for iSub = 1:length(goodSubs)%1:length(goodSubs)%length(subdirs)
    for rwd=1:2
        trPulse{iSub,rwd} = reshape(rwdPulseTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
        subMeanPulse(iSub,rwd,:) = nanmean(trPulse{iSub,rwd},2);
        subPulseStd(iSub,rwd) = std(subMeanPulse(iSub,rwd,:));%std amplitude of mean
        subPulseMean(iSub,rwd) = mean(subMeanPulse(iSub,rwd,:));%baseline defined as mean heart rate
    end
end

%% PERMUTATIONS

for iSub = 1:length(goodSubs)
    for rwd=1:2
        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!

    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    subPulse = [trPulse{iSub,1} trPulse{iSub,2}];

    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for rwd=1:2
            permPulse = subPulse(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permSubPulseStd(iSub,rwd,p) = std(nanmean(permPulse,2));%std amplitude of average pulse
            permSubPulseMean(iSub,rwd,p) = mean(nanmean(permPulse,2));%mean heart rate
        end 
    end
end

permSubPulseStdDiff = squeeze(permSubPulseStd(:,1,:) - permSubPulseStd(:,2,:));
permSubPulseMeanDiff = squeeze(permSubPulseMean(:,1,:) - permSubPulseMean(:,2,:));

%PULSE P VALUES
realSubPulseStdDiff = subPulseStd(:,1) - subPulseStd(:,2);
pval_pulseStd = sum(mean(permSubPulseStdDiff)>=mean(realSubPulseStdDiff))/nperms;
realSubPulseMeanDiff = subPulseMean(:,1) - subPulseMean(:,2);
pval_pulseMean = sum(mean(permSubPulseMeanDiff)>=mean(realSubPulseMeanDiff))/nperms;


pval_pulseStd
pval_pulseMean  

