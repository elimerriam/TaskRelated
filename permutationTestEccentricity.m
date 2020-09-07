


close all; clear all;

dataFolder = '/Users/rothzn/min/';
nperms=10000;%number of permutations


onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
ConcatProj=1;
curFolder = pwd;
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

load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString  ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd');

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
linewidth = 1;
markersize=10;
ROIs = 1:length(roiNames);
ROIs = [1 2];
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22


%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;


%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

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
                elseif onlyCorrect ==2 % ONLY INCORRECT
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
                temp = fft(reshapedTrials);
                roiBinFftAmp = abs(temp(2,:));
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshapedTrials;
               
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                %trial-by-trial variability per subject
                subBinStd(iSub,iRoi,ibin,rwd,:) = std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2);
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftPh = angle(temp(2,:));
                subBinFftAmpStd(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp);
                subBinFftPhStd(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh');
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                subBinStdMean(iSub,iRoi,ibin,rwd) = mean(singleTrialStd);

            end

        end
    end
end

subBinAmp = squeeze(std(subBinResponse,0,5));%get mean timecourse amplitude per subject
binAmp = squeeze(mean(subBinAmp,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binAmp(:,:,1) - binAmp(:,:,2));%(iRoi,ibin)

binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)

%mean trial-by-trial variability
binStdMean = squeeze(mean(subBinStd,1)); % mean variability across subjects. binStd(iRoi,ibin,rwd,timepoint)
binStdDiff = squeeze(binStdMean(:,:,2,:) - binStdMean(:,:,1,:));

%% PERMUTATIONS

for iSub = 1:length(goodSubs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        if rwd==1
            firstTrial(rwd)=1;
        else %rwd==2
            firstTrial(rwd)=numTrials(iSub,1)+1;
        end
    end
    %concatenate trials for both reward levels
    for iRoi=ROIs
        for ibin=1:nbins
            roiBinTrials{iRoi,ibin} = [subBinTrialResponse{iSub,iRoi,ibin,1} subBinTrialResponse{iSub,iRoi,ibin,2}];
        end
    end
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));%random order of all trials (low+high rwd)
        for iRoi= ROIs
            for ibin=1:nbins
                for rwd=1:2
                    permTrials = roiBinTrials{iRoi,ibin}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                    permSubMeanTC(iSub,iRoi,ibin,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                    permSubResp(iSub,iRoi,ibin,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,ibin,rwd,p,:)));%amplitude of task response
                    
                    
                    %trial-by-trial variability per subject
                    permSubStd(iSub,iRoi,ibin,rwd,p,:) = std(permTrials,0,2);%timecourse of variability over trials

                    %trial-by-trial FFT amp and phase variability per subject
                    temp = fft(permTrials);
                    permFftAmp = abs(temp(2,:));
                    permFftPh = angle(temp(2,:));
                    permSubAmpStd(iSub,iRoi,ibin,rwd,p) = std(permFftAmp);
                    permSubPhStd(iSub,iRoi,ibin,rwd,p) = circ_std(permFftPh');
                   
                    singleTrialStd = std(permTrials);
                    permSubStdMean(iSub,iRoi,ibin,rwd,p) = mean(singleTrialStd);
                end
            end
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for ibin=1:nbins
            for rwd = 1:2
                realSubMeanTC(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                
            end
            realSubDiff(iSub,iRoi,ibin) = std(realSubMeanTC(iSub,iRoi,ibin,1,:)) - std(realSubMeanTC(iSub,iRoi,ibin,2,:));%difference between high std and low std
            permSubDiff(iSub,iRoi,ibin,:) = squeeze(permSubResp(iSub,iRoi,ibin,1,:) - permSubResp(iSub,iRoi,ibin,2,:));%difference between high std and low std
            pVal_sub(iSub, iRoi,ibin) = sum(permSubDiff(iSub,iRoi,ibin,:) > realSubDiff(iSub,iRoi,ibin))/nperms;
            
            permSubStdDiff(iSub,iRoi,ibin,:,:) = squeeze(permSubStd(iSub,iRoi,ibin,2,:,:) - permSubStd(iSub,iRoi,ibin,1,:,:));%difference between high and low reward varability 

            realSubAmpStdDiff(iSub,iRoi,ibin) = subBinFftAmpStd(iSub,iRoi,ibin,2) - subBinFftAmpStd(iSub,iRoi,ibin,1);
            realSubPhStdDiff(iSub,iRoi,ibin) = subBinFftPhStd(iSub,iRoi,ibin,2) - subBinFftPhStd(iSub,iRoi,ibin,1);
            realSubStdMeanDiff(iSub,iRoi,ibin) = subBinStdMean(iSub,iRoi,ibin,2) - subBinStdMean(iSub,iRoi,ibin,1);

            permSubAmpStdDiff(iSub,iRoi,ibin,:) = squeeze(permSubAmpStd(iSub,iRoi,ibin,2,:) - permSubAmpStd(iSub,iRoi,ibin,1,:));
            permSubPhStdDiff(iSub,iRoi,ibin,:) = squeeze(permSubPhStd(iSub,iRoi,ibin,2,:) - permSubPhStd(iSub,iRoi,ibin,1,:));
            permSubStdMeanDiff(iSub,iRoi,ibin,:) = squeeze(permSubStdMean(iSub,iRoi,ibin,2,:) - permSubStdMean(iSub,iRoi,ibin,1,:));
        end

    end
end
permStdDiff = squeeze(mean(permSubStdDiff));%mean over subjects. permStd(iRoi,ibin,p,timepoint)

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealStdMeanDiff = squeeze(mean(realSubStdMeanDiff));%difference in single trial std

meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermStdMeanDiff= squeeze(mean(permSubStdMeanDiff));%mean over subjects
        
hypoth = pVal_sub<0.05;
for iRoi= ROIs
    for ibin=1:nbins
        %RFX
        meanPermDiff(iRoi,ibin,:) = mean(permSubDiff(:,iRoi,ibin,:));%average over subjects
        meanRealDiff(iRoi,ibin) = mean(realSubDiff(:,iRoi,ibin));%average over subjects
        pVal_rfx(iRoi,ibin) = sum(meanPermDiff(iRoi,ibin,:) > meanRealDiff(iRoi,ibin))/nperms;
        
        %variability
        for t=1:trialLength
            pVal_var_timecourse(iRoi,ibin,t) = sum( permStdDiff(iRoi,ibin,:,t) > binStdDiff(iRoi,ibin,t) )/nperms;
        end
        pVal_var(iRoi,ibin) = sum( mean(permStdDiff(iRoi,ibin,:,:),4) > mean(binStdDiff(iRoi,ibin,:),3) )/nperms;
        
        %FFT variability

        pVal_ampVar(iRoi,ibin) = sum(meanPermAmpStdDiff(iRoi,ibin,:) > meanRealAmpStdDiff(iRoi,ibin))/nperms;
        pVal_phVar(iRoi,ibin) = sum(meanPermPhStdDiff(iRoi,ibin,:) > meanRealPhStdDiff(iRoi,ibin))/nperms;

        %single trial STD amplitude mean
        pVal_stdMean(iRoi,ibin) = sum(meanPermStdMeanDiff(iRoi,ibin,:) > meanRealStdMeanDiff(iRoi,ibin))/nperms;
    end
end
pVal_rfx;%response amplitude
pVal_var;%timepoint variability
pVal_ampVar;%amplitude variability
pVal_phVar;%temporal/phase variability


%%
iroi=2; %right EVC
'timepoint variability p-value: '
pVal_var(iroi,:)

