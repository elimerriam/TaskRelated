% permutationTest
%
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020).
% Task-related activity in human visual cortex.
%
%   uses file created by saveTaskData_localizer.m
%
%   usage: permutationTest()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: compute task-related response amplitude, timepoint variability, amplitude variability,
%   and temporal variability, and use permutations to test for significant
%   difference between high and low reward.
%
%

function[] = permutationTest()

removeLocalizer = 0;%0=don't exclude voxels identified by the localizer. 1 = exclude localizer voxels.
locThresh=0.3;%localizer coherence threshold for excluding voxels that respond to the stimulus

nperms=10000;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%1=data is z-scored, 0=no z-scoring
ConcatProj = 1;%1=project out global mean. 0=don't project out mean

curFolder = pwd;
dataFolder = '/Users/rothzn/min/';

load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
rng(randSeed);


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
noLocStr = '';
if removeLocalizer~=0
    thrsh = num2str(locThresh,'%.2f');
    noLocStr = ['_noLoc' thrsh([1 3:4]) ];
    if removeLocalizer==-1
        noLocStr = ['_onlyLoc' thrsh([1 3:4]) ];
    end
    
    load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString ConcatProjStr noLocStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
        'roiTC', 'allTrials', ...
        'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
        'expName','stairThresh','eccen','ang','areas','trialLength',...
        'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
        'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
        'globalMean','regressBetasGlobal','runRwd',...
        'subRoiRuns','runMeanFFT',...
        'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
        'voxTrials','voxGoodTrials','meanVoxTrial',...
        'locThresh','locPercent','locCo','locPh','locAmp','origVoxNum','finalVoxNum',...
        'maxRT');
else
    load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString ConcatProjStr noLocStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
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
end

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
ROIs = 1:length(roiNames);

goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22

rwdString = {'H','L'};


i=0;
%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;

%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh);
std(subThresh);


%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

for iSub = 1:length(goodSubs)%length(subdirs)
    %     iSub = goodSubs(i)
    for iRoi=1:length(roiNames)
        for rwd=1:2
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};%10 timepoints X number of trials
            
            %trial-by-trial variability per subject
            subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
            temp = fft(reshapedTrials);
            roiFftAmp = abs(temp(2,:));%Fourier amplitude per trial
            roiFftPh = angle(temp(2,:));%Fourier phase per trial
            trialRoiFftPhVec{iSub,iRoi,rwd} = roiFftPh;
            trialRoiFftAmpVec{iSub,iRoi,rwd} = roiFftAmp;
            trialRoiAmpVec{iSub,iRoi,rwd} = std(reshapedTrials);%std amplitude per trial
            subFftAmpStd(iSub,iRoi,rwd) = std(roiFftAmp);%variability of Fourier amplitude
            subFftMeanAmp(iSub,iRoi,rwd) = mean(roiFftAmp);%Mean of Fourier amplitude
            subFftMeanPh(iSub,iRoi,rwd) = circ_mean(roiFftPh');%circular mean of Fourier phase
            subMeanStd(iSub,iRoi,rwd) = mean(std(reshapedTrials));%mean of amplitude measured by std
            subStdVar(iSub,iRoi,rwd) = std(std(reshapedTrials));%variability of amplitude measured by std
            subMedianStd(iSub,iRoi,rwd) = median(std(reshapedTrials));
            
            subFftPhStd(iSub,iRoi,rwd) = circ_std(roiFftPh');%circular std of phase
            
            if iSub==goodSubs(1)
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
            
        end
    end
end


for rwd=1:2
    for iRoi = 1:length(roiNames)
        meanResponse(iRoi,rwd,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(iRoi,rwd,:) = std(subResponse(goodSubs,iRoi,rwd,:));
        meanTimepointStd(iRoi, rwd,:) = mean(subTimepointStd(:,iRoi,rwd,:));%mean trial-by-trial timepoint variability
        stdTimepointStd(iRoi, rwd,:) = std(subTimepointStd(:,iRoi,rwd,:));%std of trial-by-trial timepoint variability
    end
end


%% PERMUTATIONS
for iSub = 1:length(goodSubs)
    for rwd=1:2
        temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);%exclude the first trial and the last trial
        trialCorrectVec{iSub,rwd} = temp(:); %include trials with no response
        subPropCorrect(iSub,rwd) = sum(trialCorrectVec{iSub,rwd})/length(trialCorrectVec{iSub,rwd});
        
        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
        trialResponseVec = temp(:);
        temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
        allTrialRTvec = temp(:);
        
        if onlyCorrect==0
            trialRTvec{iSub,rwd} = temp(trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);%only trials with a response, we're assuming onlyCorrect==0
        elseif onlyCorrect==1 %only correct
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==1 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        elseif onlyCorrect==2 %only incorrect
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==0 & trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        else
            trialRTvec{iSub,rwd} = temp(:);
        end
        
    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    
    for iRoi=ROIs
        %concatenate trials across reward level
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
    end
    subRT = [trialRTvec{iSub,1}; trialRTvec{iSub,2}];
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response
                
                %latency
                temp=fft(permSubMeanTC(iSub,iRoi,rwd,p,:));
                permSubMeanFftAmp(iSub,iRoi,rwd,p) = abs(temp(2));
                permSubMeanFftPh(iSub,iRoi,rwd,p) = angle(temp(2));
                
                permSubStd(iSub,iRoi,rwd,p,:) = std(permTrials,0,2);%variability timecourse over trials
                
                %trial-by-trial FFT amp and phase variability per subject
                temp = fft(permTrials);
                permFftAmp = abs(temp(2,:));
                permFftPh = angle(temp(2,:));
                permFftOther = sum(abs(temp([1 3:end],:)));
                permSubAmpStd(iSub,iRoi,rwd,p) = std(permFftAmp);
                permSubPhStd(iSub,iRoi,rwd,p) = circ_std(permFftPh');
                
                permSubMeanAmp(iSub,iRoi,rwd,p) = mean(permFftAmp);
                permSubMeanPh(iSub,iRoi,rwd,p) = circ_mean(permFftPh');
                permSubMeanStd(iSub,iRoi,rwd,p) = mean(std(permTrials));
                
                permSubStdVar(iSub,iRoi,rwd,p) = std(std(permTrials));
            end
        end
        for rwd=1:2
            %RT variability
            permRT = subRT(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1),:);
            permSubRTvar(iSub,rwd,p) = std(permRT);
            permSubRT(iSub,rwd,p) = mean(permRT);
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        %difference in amplitude
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
            temp = fft(squeeze(realSubMeanTC(iSub,iRoi,rwd,:)));
            realSubMeanFftAmp(iSub,iRoi,rwd) = abs(temp(2));
            realSubMeanFftPh(iSub,iRoi,rwd) = angle(temp(2));
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std
        permSubDiff(iSub,iRoi,:) = squeeze(permSubResp(iSub,iRoi,1,:) - permSubResp(iSub,iRoi,2,:));%difference between high std and low std
        pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) >= realSubDiff(iSub,iRoi))/nperms;
        
        %difference in latency
        realSubFftAmpDiff(iSub,iRoi) = realSubMeanFftAmp(iSub,iRoi,1) - realSubMeanFftAmp(iSub,iRoi,2);
        realSubFftPhDiff(iSub,iRoi) = circ_dist(realSubMeanFftPh(iSub,iRoi,2), realSubMeanFftPh(iSub,iRoi,1));
        permSubFftAmpDiff(iSub,iRoi,:) = permSubMeanFftAmp(iSub,iRoi,1,:) - permSubMeanFftAmp(iSub,iRoi,2,:);
        permSubFftPhDiff(iSub,iRoi,:) = circ_dist(permSubMeanFftPh(iSub,iRoi,2,:),permSubMeanFftPh(iSub,iRoi,1,:));
        pVal_fftPh_sub(iSub, iRoi) = sum(permSubFftPhDiff(iSub,iRoi,:) >= realSubFftPhDiff(iSub,iRoi))/nperms;
        
        %FFT variability
        realSubAmpStdDiff(iSub,iRoi) = subFftAmpStd(iSub,iRoi,2) - subFftAmpStd(iSub,iRoi,1);
        realSubPhStdDiff(iSub,iRoi) = subFftPhStd(iSub,iRoi,2) - subFftPhStd(iSub,iRoi,1);
        
        realSubMeanFftAmpDiff(iSub,iRoi) = subFftMeanAmp(iSub,iRoi,2) - subFftMeanAmp(iSub,iRoi,1);
        realSubMeanFftPhDiff(iSub,iRoi) = squeeze(circ_dist(subFftMeanPh(iSub,iRoi,2),subFftMeanPh(iSub,iRoi,1)));
        realSubMeanStdDiff(iSub,iRoi) = subMeanStd(iSub,iRoi,2) - subMeanStd(iSub,iRoi,1);
        
        realSubStdVarDiff(iSub,iRoi) = subStdVar(iSub,iRoi,2) - subStdVar(iSub,iRoi,1);
        
        permSubAmpStdDiff(iSub,iRoi,:) = squeeze(permSubAmpStd(iSub,iRoi,2,:) - permSubAmpStd(iSub,iRoi,1,:));
        permSubPhStdDiff(iSub,iRoi,:) = squeeze(permSubPhStd(iSub,iRoi,2,:) - permSubPhStd(iSub,iRoi,1,:));
        
        permSubMeanAmpDiff(iSub,iRoi,:) = squeeze(permSubMeanAmp(iSub,iRoi,2,:) - permSubMeanAmp(iSub,iRoi,1,:));
        permSubMeanPhDiff(iSub,iRoi,:) = squeeze(circ_dist(permSubMeanPh(iSub,iRoi,2,:),permSubMeanPh(iSub,iRoi,1,:)));
        permSubMeanStdDiff(iSub,iRoi,:) = squeeze(permSubMeanStd(iSub,iRoi,2,:) - permSubMeanStd(iSub,iRoi,1,:));
        
        permSubStdVarDiff(iSub,iRoi,:) = squeeze(permSubStdVar(iSub,iRoi,2,:) - permSubStdVar(iSub,iRoi,1,:));
    end
    for rwd=1:2
        %RT variability
        realSubRTvar(iSub,rwd) = std(trialRTvec{iSub,rwd});
        %RT
        realSubRT(iSub,rwd) = mean(trialRTvec{iSub,rwd});
    end
    
end

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealFftAmpDiff = squeeze(mean(realSubMeanFftAmpDiff));%mean over subjects
meanRealFftPhDiff = squeeze(mean(realSubMeanFftPhDiff));%mean over subjects
meanRealStdVarDiff = squeeze(mean(realSubStdVarDiff));%mean over subjects
meanRealMeanStdDiff = squeeze(mean(realSubMeanStdDiff));%mean trials amplitude, mean over subjects

meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermMeanAmpDiff = squeeze(mean(permSubMeanAmpDiff));%mean over subjects
meanPermMeanPhDiff = squeeze(mean(permSubMeanPhDiff));%mean over subjects
meanPermMeanStdDiff = squeeze(mean(permSubMeanStdDiff));%mean over subjects

meanPermStdVarDiff = squeeze(mean(permSubStdVarDiff));%mean over subjects


for iRoi= ROIs
    %Amplitude
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_stdAmp(iRoi) = sum(meanPermDiff(iRoi,:) >= meanRealDiff(iRoi))/nperms;
    
    %FFT Phase
    meanPermPhDiff(iRoi,:) = mean(permSubFftPhDiff(:,iRoi,:));%average over subjects
    meanRealPhDiff(iRoi) = mean(realSubFftPhDiff(:,iRoi));%average over subjects
    pVal_ph(iRoi) = sum(meanPermPhDiff(iRoi,:) >= meanRealPhDiff(iRoi))/nperms;
    meanPermAmpDiff(iRoi,:) = mean(permSubFftAmpDiff(:,iRoi,:));%average over subjects
    meanRealAmpDiff(iRoi) = mean(realSubFftAmpDiff(:,iRoi));%average over subjects
    pVal_amp(iRoi) = sum(meanPermAmpDiff(iRoi,:) >= meanRealAmpDiff(iRoi))/nperms;
    
    %variability
    permStd(iRoi,:,:,:) = mean(permSubStd(:,iRoi,:,:,:));%mean over subjects. permStd(iRoi,iRoi,rwd,p)
    for t=1:trialLength
        pVal_var_timecourse(iRoi,t) = sum( (permStd(iRoi,2,:,t)-permStd(iRoi,1,:,t)) >= (meanTimepointStd(iRoi,2,t)-meanTimepointStd(iRoi,1,t)))/nperms;
    end
    pVal_var(iRoi) = sum( mean(permStd(iRoi,2,:,:)-permStd(iRoi,1,:,:),4) >= mean(meanTimepointStd(iRoi,2,:)-meanTimepointStd(iRoi,1,:),3) )/nperms;
    
    %FFT variability
    pVal_ampVar(iRoi) = sum(meanPermAmpStdDiff(iRoi,:) >= meanRealAmpStdDiff(iRoi))/nperms;
    pVal_phVar(iRoi) = sum(meanPermPhStdDiff(iRoi,:) >= meanRealPhStdDiff(iRoi))/nperms;
    
    pVal_phMean(iRoi) = sum(meanPermMeanPhDiff(iRoi,:) >= meanRealFftPhDiff(iRoi))/nperms;
    pVal_ampMean(iRoi) = sum(meanPermMeanAmpDiff(iRoi,:) >= meanRealFftAmpDiff(iRoi))/nperms;
    pVal_stdMean(iRoi) = sum(meanPermMeanStdDiff(iRoi,:) >= meanRealMeanStdDiff(iRoi))/nperms;
    
    pVal_stdVar(iRoi) = sum(meanPermStdVarDiff(iRoi,:) >= meanRealStdVarDiff(iRoi))/nperms;
end
pVal_stdAmp;
pVal_ph;
pVal_amp;
pVal_var;

pVal_phVar;
pVal_ampMean; %mean trial amplitude measured by fft (per trial), low minus high
pVal_stdMean; %mean trial amplitude measured by std (per trial), low minus high
pVal_phMean;

%RT variability
permRTvar = squeeze(mean(permSubRTvar));
permRTvarDiff = permRTvar(2,:)-permRTvar(1,:);
realRTvar = mean(realSubRTvar);
realRTvarDiff = realRTvar(2) - realRTvar(1);
pVal_RTvar = sum(permRTvarDiff >= realRTvarDiff)/nperms;
%RT
permMeanRT = squeeze(mean(permSubRT));
permRTdiff = permMeanRT(2,:)-permMeanRT(1,:);
realMeanRT = mean(realSubRT);
realRTdiff = realMeanRT(2) - realMeanRT(1);
pVal_RT = sum(permRTdiff >= realRTdiff)/nperms;

pVal_RTvar
realSubRTvarDiff = realSubRTvar(:,2) - realSubRTvar(:,1);

pVal_RT
realSubRTdiff = realSubRT(:,2)-realSubRT(:,1);

pVal_stdVar; %variability of amplitude measured by std

%%
'mean accuracy'
mean(subPropCorrect)
'accuracy std'
std(subPropCorrect)
'mean RT'
mean(realSubRT)
'RT std'
std(realSubRT)

%%
%statistics of how many localizer voxels were excluded
if removeLocalizer
    mean(origVoxNum(goodSubs,:) - finalVoxNum(goodSubs,:));
    mean((origVoxNum(goodSubs,:) - finalVoxNum(goodSubs,:))./origVoxNum(goodSubs,:));
    percentExcluded = 100 - locPercent(goodSubs,:);
    mean(mean(percentExcluded,2))
    min(mean(percentExcluded,2))
    max(mean(percentExcluded,2))
end


%p-values
{'amplitude p-value', 'timepoint variability p-value','phase variability p-value','amplitude variability p-value'}
[pVal_stdAmp(2) pVal_var(2) pVal_phVar(2) pVal_stdVar(2)]

