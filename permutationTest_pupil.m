% permutationTest_pupil
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% uses file created by savePupilData.m
%
%   usage: permutationTest_pupil()
%   by: zvi roth
%   date: 9/6/2020
%   purpose: compute phasic and tonic pupil size, and test for
%   significantly larger pupil size for high reward using permutation test
%

function[] = permutationTest_pupil()

dataFolder = '/Volumes/TaskDrive/taskRelatedBehavioral/';

nperms=10000;%number of permutations


onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
tic
load([dataFolder 'behavioralData' onlyCorrectString '.mat'], 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','subRTvar','rwdLevel','numTrials',...
    'trialCorrectness' ,'trialResponse','trialRT','propCorrect','stairThresh');
%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh);
std(subThresh);

%%
goodSubs = 1:length(subFolders);
minLength=3000;%keeping first 3sec of each trial
baseT = 50;%baseline is first 50 ms of each trial
%%
clear permSubMeanTC
for iSub = 1:length(subFolders)%length(subdirs)
    for rwd=1:2
        numTrials(iSub,rwd) = size(rwdPupil{iSub,rwd},1);
    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    subTrials = [rwdPupil{iSub,1}(:,1:minLength); rwdPupil{iSub,2}(:,1:minLength)];
    subRT = [trialRT{iSub,1}; trialRT{iSub,2}];
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for rwd=1:2
            permTrials = subTrials(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1),:);% #trials X timepoints
            permSubMeanTC = nanmean(permTrials(:,1:minLength));
            permSubMeanStd(iSub,rwd,p) = nanstd(squeeze(permSubMeanTC));%std amplitude of phasic pupil response
            permSubMeanBaseline(iSub,rwd,p) = nanmean(squeeze(permSubMeanTC(1:baseT)));%baseline pupil size
        end
    end
    for rwd=1:2
        subRwdTrials = subTrials(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1,1:minLength);
        realSubMeanTC(iSub,rwd,:) =  nanmean(subRwdTrials);
        realSubStd(iSub,rwd) = nanstd(squeeze(realSubMeanTC(iSub,rwd,:)));%std amplitude of phasic pupil response
        realSubBaseline(iSub,rwd) = nanmean(squeeze(realSubMeanTC(iSub,rwd,1:baseT)));%baseline pupil size
        
        realSubRTvar(iSub,rwd) = std(trialRT{iSub,rwd});
        realSubRT(iSub,rwd) = mean(trialRT{iSub,rwd});
        
    end
    %high reward minus low reward
    realSubStdDiff(iSub) = realSubStd(iSub,1) - realSubStd(iSub,2);
    realSubBaselineDiff(iSub) = realSubBaseline(iSub,1) - realSubBaseline(iSub,2);
    permSubStdDiff(iSub,:) = permSubMeanStd(iSub,1,:) - permSubMeanStd(iSub,2,:);
    permSubBaselineDiff(iSub,:) = permSubMeanBaseline(iSub,1,:) - permSubMeanBaseline(iSub,2,:);
end

%%
%average over subjects
meanPermStdDiff = nanmean(permSubStdDiff);
meanRealStdDiff = nanmean(realSubStdDiff);
meanPermBaselineDiff =nanmean(permSubBaselineDiff);
meanRealBaselineDiff = nanmean(realSubBaselineDiff);

%get p-values
pVal_phasic = sum(meanPermStdDiff >= meanRealStdDiff)/nperms;
pVal_tonic = sum(meanPermBaselineDiff >= meanRealBaselineDiff)/nperms;

pVal_phasic
pVal_tonic
