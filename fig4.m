% fig4
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% Uses data saved by saveTaskData.m
%
%   usage: fig4()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: create figure 4, task-related response amplitude
%

function[] = fig4()


onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%1=data is z-scored, 0=no z-scoring
ConcatProj = 1;%1=project out global mean. 0=don't project out mean

curFolder = pwd;
dataFolder = '';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

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
    fmriUnits = 'std';
end
pscString = '';
if ~toZscore
    pscString = '_percSigCh';
    fmriUnits = 'percent';
end
ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
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
TR=1.5;

plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
fontsize=9;
ROIs = [2];%right EVC

dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;

ntrials=15;
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22


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
                    %                 numTrials = sum(trialCorrectnessVec);
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
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                %trial-by-trial variability per subject
                subBinVar(iSub,iRoi,ibin,rwd) = mean(std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2));
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftPh = angle(temp(2,:));
                subBinFftAmpVar(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp);
                subBinFftPhVar(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh');
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                subBinStdVar(iSub,iRoi,ibin,rwd) = std(singleTrialStd);
                subBinStdMean(iSub,iRoi,ibin,rwd) = mean(singleTrialStd);

            end
        end
    end
end

subBinAmp = std(subBinResponse,0,5);%get mean timecourse amplitude per subject.(iRoi,ibin,rwd,1)
binMeanAmp = squeeze(mean(subBinAmp,1)); %mean over subjects. (iRoi,ibin,rwd)
binStdAmp = squeeze(std(subBinAmp,0,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binMeanAmp(:,:,1) - binMeanAmp(:,:,2));%(iRoi,ibin)

subBinDiff = squeeze(subBinAmp(:,:,:,1) - subBinAmp(:,:,:,2));%(iSub,iRoi,ibin)
binMeanDiff = squeeze(mean(subBinDiff));
binStdDiff = squeeze(std(subBinDiff));
binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)




%%
i=0;

rows=2;
cols = 9;
subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};

for r= 1:length(ROIs)
    i=i+1; figure(i); clf; 
    %amplitude
    subplot(rows,cols,subplots{1});
    iRoi=ROIs(r);
    for rwd=1:2
        dsErrorsurface(binCenters, squeeze(binMeanAmp(iRoi,:,rwd)), squeeze(binStdAmp(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
    end
    for rwd=1:2
        plot(binCenters, squeeze(binMeanAmp(iRoi,:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    subplot(rows,cols,subplots{2});
    dsErrorsurface(binCenters, binMeanDiff(iRoi,:), binStdDiff(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binMeanAmp(iRoi,:,1) -binMeanAmp(iRoi,:,2)),'k.','linewidth',linewidth,'markersize',markersize);
    hline(0);

end

%% bar plot of subjects' amplitudes
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;

for r = 1:length(ROIs)
    figure(r);
    subplot(rows,cols,subplots{3})
    iROI = ROIs(r);
    clear smallSubAmp subPh
    for iSub=1:length(goodSubs)
        for rwd=1:2
            smallSubAmp(iSub,rwd) = std(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:)));
            f=angle(fft(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:))));
            subPh(iSub,rwd) = f(2);
        end
    end
    smallSubDiff = smallSubAmp(:,1) - smallSubAmp(:,2);
    minSubDiff = min(smallSubDiff);
    maxSubDiff= max(smallSubDiff);
    subjects = size(smallSubAmp,1);
    rewards = size(smallSubAmp,2);
    [rwdNum, subNum] = meshgrid(1:rewards, 1:subjects);
    scatterCmap=cool;
%     scatterCmap = scatterCmap(end:-1:1,:);%invert color map
    colormap(scatterCmap);

    l = size(scatterCmap,1);
    for iSub=1:length(goodSubs)
        subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
        plot(1:rewards,smallSubAmp(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),smallSubAmp(:),markerSize, [subColor; subColor] ,'filled');   
    % MEAN
    for iRwd=1:rewards
        scatter(iRwd,mean(smallSubAmp(:,iRwd)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(smallSubAmp),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
end



%% Formatting
if toZscore
    responseStr = ['response amplitude (std)'];
    deltaResponseStr = ['\Delta response amplitude (std)'];
else
    responseStr = {'response amplitude'; '(% change image intensity)'};
    deltaResponseStr = {'\Delta response amplitude';  '(% change image intensity)'};
end

for r = 1:length(ROIs)
    figure(r);
    for isubplot=[1:3]% 5:7]%1:length(subplots)
        subplot(rows,cols,subplots{isubplot});
        switch isubplot
            case 1
                ylabel(responseStr);
            case 2
                ylabel(deltaResponseStr);
            case 3
                ylabel(responseStr);
        end
        
        if isubplot==3
            xlabel('reward');
            drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -80/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                'xTick',[1 2], 'xTickLabel', {'high','low'},...
                'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
                'yAxisMajorTickLen',-4/32);
            legend off
        else
            
            xlabel('eccentricity (deg)');
            drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -16/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                'labelFontSize',fontsize);
            legend off
            axis square
        end
    end
    set(gcf,'position',[10 10 28 	14]);
end

print('-painters','-dpdf',['fig4_' ConcatProjStr  pscString '.pdf']);


