% fig2
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% Uses data saved by savePupilData.m, and saveTaskDataPhysio
%
%   usage: fig2()
%   by: zvi roth
%   date: 9/3/2020
%   purpose: create figure 2, mean pupil size for high- and low-reward
%   eye-tracking runs, mean heart-rate for high- and low-reward runs, and
%   pulse-to-BOLD kernel
%

function[] = fig2()

eyeDataFolder = '';%folder containing all eye-tracking data
dataFolder = '';%folder containing all fMRI data

onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%1=data is z-scored, 0=no z-scoring

iRoi=2;
rows=1;
cols=3;
subplots = {1,2,3,4};
linelength = 0.2;
linewidth = 1;
markersize = 10;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
fontsize=9;
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;

i=0; %figure counter
i=i+1; figure(i); clf
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
tic
load([eyeDataFolder 'behavioralData' onlyCorrectString '.mat'], 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','subRTvar','rwdLevel','numTrials',...
    'trialCorrectness' ,'trialResponse','trialRT','propCorrect','stairThresh');
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};

%%
minLength = 2000;
for rwd=1:2
    groupAll{rwd}(1,:) = meanPupil{1,rwd}(1:minLength);
    for iSub=2:length(subFolders)
        groupAll{rwd}(iSub,:) =  meanPupil{iSub,rwd}(1:minLength);
    end
    groupMean{rwd} = mean(groupAll{rwd});
    groupStd{rwd} = std(groupAll{rwd});
end
clear diffPupil groupAll
for iSub = 1:length(subFolders)
%     minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil(iSub,:) = meanPupil{iSub,1}(1:minLength) - meanPupil{iSub,2}(1:minLength);
end
subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface((1:minLength)*2, groupMean{rwd}(1:minLength), groupStd{rwd}(1:minLength)./sqrt(length(subFolders)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,groupMean{rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end
%% HEART RATE
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
zScoreString = '';
fmriUnits = '% change image intensity';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end

%ConcatProj = 1, project out global mean signal,
%ConcatProj = 0, do not project out global mean signal
for ConcatProj=0:1
    clear subPulseTrial numBinVoxels binMeanTserie subBinTrialResponse
    clear binPulseKernel binPulseResidualTC binPhysioKernel binPhysioResidualTC binRespKernel binRespResidualTC
    ConcatProjStr = '';
    if ConcatProj
        ConcatProjStr = 'ConcatProj';
    end
    load([dataFolder 'rwdTC_physio' onlyCorrectString zScoreString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
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
        'allGoodTrials','allGoodTRs');
    
    ROIs = 1:length(roiNames);
    eccMin = 0.2;
    eccMax = 70;
    nbins = 12;
    TR=1.5;
    binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
    nbins = length(binBorders)-1;
    for ibin=2:length(binBorders)
        binCenters(ibin-1) = (binBorders(ibin)+binBorders(ibin-1))/2;
    end
    goodSubs = [1:length(subFolders)];%already excludes s0022

    for iSub = 1:length(goodSubs)
        for iRoi=1:length(roiNames)
            for rwd=1:2
                %RESPIRATION
                trRV{iSub,rwd} = reshape(rwdRvTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
                subMeanRV(iSub,rwd,:) = nanmean(trRV{iSub,rwd},2);
                subRespStd(iSub,rwd) = std(subMeanRV(iSub,rwd,:));%std amplitude of mean
                subRespVar(iSub,rwd,:) = nanstd(trRV{iSub,rwd},0,2);%timepoint variability
                respStd = nanstd(trRV{iSub,rwd});%std amp per trial
                subRespStdVar(iSub,rwd) = std(respStd);
                f = fft(trRV{iSub,rwd});
                subRespPhVar(iSub,rwd) = nanstd(angle(f(2,:)));
                subRespAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
                %PULSE
                trPulse{iSub,rwd} = reshape(rwdPulseTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
                subMeanPulse(iSub,rwd,:) = nanmean(trPulse{iSub,rwd},2);
                subPulseStd(iSub,rwd) = std(subMeanPulse(iSub,rwd,:));%std amplitude of mean
                subPulseVar(iSub,rwd,:) = nanstd(trPulse{iSub,rwd},0,2);%timepoint variability
                pulseStd = nanstd(trPulse{iSub,rwd});%std amp per trial
                subPulseStdVar(iSub,rwd) = std(pulseStd);
                f = fft(trPulse{iSub,rwd});
                subPulsePhVar(iSub,rwd) = nanstd(angle(f(2,:)));
                subPulseAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
                
                %REGRESS OUT PHYSIO
                designMatPulse{goodSubs(iSub),rwd}= designMatPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                pulseKernel(iSub,iRoi,rwd,:) = designMatPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);%this is only good trials!
                pulseResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{goodSubs(iSub),rwd};
                for ilag=1:size(designMatPulse{goodSubs(iSub),rwd},1)
                    pulseBoldCorr{ConcatProj+1}(iSub,iRoi,rwd,ilag) = corr(subTrialResponse{goodSubs(iSub),iRoi,rwd}(:),designMatPulse{goodSubs(iSub),rwd}(ilag,:)');
                end
                designMatResp{goodSubs(iSub),rwd} = designMatResp{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                respKernel(iSub,iRoi,rwd,:) = designMatResp{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
                respResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{goodSubs(iSub),rwd};
                
                designMatRespPulse{goodSubs(iSub),rwd} = designMatRespPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                physioKernel(iSub,iRoi,rwd,:) = designMatRespPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
                physioResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{goodSubs(iSub),rwd};
                subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshape(physioResidualTC{iSub,iRoi,rwd}(:), trialLength, length(physioResidualTC{iSub,iRoi,rwd}(:))/trialLength);
                
                
            end
        end
    end

    iRoi = 2;%ipsilateral EVC

    %heart rate
    for iSub=1:length(subFolders)
        for rwd=1:2
            temp = reshape(rwdPulseTC{iSub,rwd},ecgTrial,[]);
            subPulseTrial(iSub,rwd,:) = nanmean(temp,2);
        end
    end
    subPulseTrial = subPulseTrial*ecgSampleRate;%changing to beats-per-second
    meanPulseTrial = squeeze(nanmean(subPulseTrial));
    stdPulseTrial = squeeze(nanstd(subPulseTrial));
    subPulseDiff = squeeze(subPulseTrial(:,1,:) - subPulseTrial(:,2,:));
    meanPulseDiff = squeeze(nanmean(subPulseDiff));
    stdPulseDiff = squeeze(nanstd(subPulseDiff));
    
    subplot(rows,cols,subplots{3})
    lineColor = [0 0 0];
    if ConcatProj==0
        linestyle = '-';
    else
        linestyle = '--';
    end
    meanRwdPulseKernel = squeeze(mean(pulseKernel,3));%average across rwd, per subject
    meanPulseKernel = squeeze(mean(meanRwdPulseKernel(:,iRoi,:)));%average across subjects
    stdPulseKernel = squeeze(std(meanRwdPulseKernel(:,iRoi,:)));
    dsErrorsurface(TR*(1:trialLength), squeeze(meanPulseKernel), squeeze(stdPulseKernel)./sqrt(size(pulseKernel,1)), dsSurfaceContrast*lineColor,dsSurfaceAlpha);
    hold on
    plot(TR*(1:trialLength), squeeze(meanPulseKernel),linestyle,'Color', lineColor,'linewidth',linewidth,'markersize',markersize);
    
    %KERNEL AMPLITUDE
    kernelAmp(ConcatProj+1) = std(squeeze(meanPulseKernel));
    kernelAmpSub(ConcatProj+1,:) = std(squeeze(meanRwdPulseKernel(:,iRoi,:))');
    
end
t = {'BOLD signal'; ['(' fmriUnits ')'] };
    ylabel(t);
kernelAmp(2)/kernelAmp(1);
kernelAmpSub(2,:)'./kernelAmpSub(1,:)';

% mean pulse rate
subplot(rows,cols,subplots{2})
for rwd=1:2
    dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanPulseTrial(rwd,:),stdPulseTrial(rwd,:)./sqrt(size(rwdPulseTC,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot((0:ecgTrial-1)/ecgSampleRate,meanPulseTrial(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
end
ylabel('pulse (beats/sec)');


%%

for isubplot=1%:2
   subplot(rows,cols,subplots{isubplot})
   xlabel('time (ms)');
   axis square
   if mod(isubplot,2)>0
       ylabel('pupil size (arb. units)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64,'labelFontSize',fontsize);
   else
       ylabel('\Delta pupil size (arb. units)');
       drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'labelFontSize',fontsize);
   end
   legend off
end
for isubplot=2:3
    subplot(rows,cols,subplots{isubplot})
    xlabel('time (s)');
    yticklabel = get(gca,'yTickLabel');
    ytick = get(gca,'yTick');
    if isubplot==3
        yticklabel = {'-0.015','0','0.03'};
        ytick = [-0.015 0 0.03];
    end
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -16/64, 'xAxisMargin', 8/64, 'yAxisMargin', 4/64,'xAxisMinMaxSetByTicks',1,...
        'labelFontSize',fontsize,'yTickLabel',yticklabel,'yTick',ytick);
    legend off
    axis square
end

%% SAVE FIGURE
set(gcf,'position',[10 10 21 8]);
print('-painters','-dpdf',['fig2.pdf']);

