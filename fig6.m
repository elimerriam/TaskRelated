% fig6
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
% Uses randSeed.mat
%
%   usage: fig6()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: create figure 6 and figure S6, simulations results
%
% left panel of fig6 uses hrfType=1, right panel uses hrfType=2.
%

function[] = fig6()

dataFolder = '';
load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
rng(randSeed);

toZscore=0;
upsampleFactor = 10;
TR = 1.5;
%set canonical HRF
modelParams = struct;
sampleDuration = TR/upsampleFactor;
sampleDelay=sampleDuration/2;
defaultParams=1;

hrfType=1;

if hrfType==1
    modelParams.x = 6;%9;
    modelParams.y = 16;
    modelParams.z = 6;
else
    % params that make the variability timecourse similar to temporal jitter
    modelParams.x = 12;
    modelParams.y = 16;
    modelParams.z = 6;
end

[modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
i=0;
plot(hrfModel);

%%
plotColorsSnr = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
numTrials = 16;
trialLength = 10;
runsPerRwd = 100;
T = numTrials*trialLength;
upT = T*upsampleFactor;
oneOverF(1,1,1,:) = (log10(1:T))*1./(1:T);
oneOverF(1,1,1,2:end) = oneOverF(1,1,1,1:end-1);%we only care about the second component onwards
oneOverF(1,1,1,T/2+1:end) = oneOverF(1,1,1,1+T/2:-1:2);%make it symmetric



%%

npoints = 10;%number of different noise levels
%noise levels for independent 1/f noise
snr = logspace(2, -0.5, npoints);
snr(1) = Inf;
%noise levels for temporal jitter
jitter = linspace(0, pi/2, npoints);
%noise levels for amplitude jitter
ampJitter = logspace(-2,0,npoints);
ampJitter(1) = 0;

for isnr=1:length(snr)
    plotColorsSnr{isnr} = [1 - (isnr-1)/(length(snr)-1), 0, (isnr-1)/(length(snr)-1)];
end
for ijitter=1:length(jitter)
    plotColorsJitter{ijitter} = [1 - (ijitter-1)/(length(jitter)-1), 0, (ijitter-1)/(length(jitter)-1)];
end
for iampJitter=1:length(ampJitter)
    plotColorsAmpJitter{iampJitter} = [1 - (iampJitter-1)/(length(ampJitter)-1), 0, (iampJitter-1)/(length(ampJitter)-1)];
end

taskAmp = 1;
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            for r=1:runsPerRwd
                runTC = zeros(1,upT);
                
                %temporal jitter
                taskTiming = 1:trialLength*upsampleFactor:upT;%temporal jitter will be added here, after upsampling
                noisyTiming = taskTiming + (jitter(ijitter)*upsampleFactor*trialLength/(2*pi))*randn(size(taskTiming));
                noisyTiming(noisyTiming<1) = 1;                
                noisyTiming(noisyTiming>T*upsampleFactor) = T;  
                
                %amplitude jitter
                runTC(ceil(noisyTiming)) = ones;
                runTC(ceil(noisyTiming)) = taskAmp*(1+tanh(ampJitter(iampJitter)*randn(numTrials,1)));
                
                temp = conv(runTC,hrfModel);
                runTC = temp(1:upT);%crop end
                rwdSignal(isnr,ijitter,iampJitter,:,r) = runTC(1:upsampleFactor:end);%downsample                
       
            end
            %1/f independent noise
            n = (taskAmp/snr(isnr))*randn(size(rwdSignal(isnr,ijitter,iampJitter,:,:)));
            rwdNoise(isnr,ijitter,iampJitter,:,:) = ifft(repmat(oneOverF,1,1,1,1,runsPerRwd).*fft(n));
        end
    end
end

%%
rwdTC = rwdSignal + rwdNoise;
rwdTC = rwdTC*upsampleFactor;
if toZscore
    rwdTC = zscore(rwdTC,0,4);
end
rwdTC = rwdTC(:,:,:,trialLength+1:end,:);%junk first cycle
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            for r=1:runsPerRwd
                %filter each run
                rwdTC(isnr,ijitter,iampJitter,:,r) = filterTC(squeeze(rwdTC(isnr,ijitter,iampJitter,:,r)));
            end
        end
    end
end

rwdTrials = reshape(rwdTC,length(snr),length(jitter),length(ampJitter),trialLength,[]);%(isnr,ijitter,iampjitter,trialLength,numTrials)
meanTrial = mean(rwdTrials,5);
meanRun = squeeze(mean(rwdTC,5));
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            f = squeeze(fft(rwdTrials(isnr,ijitter,iampJitter,:,:)));
            fftTrialAmp(isnr,ijitter,iampJitter,:) = abs(f(2,:));
            fftTrialPh(isnr,ijitter,iampJitter,:) = angle(f(2,:));
            fftRun(isnr,ijitter,iampJitter,:) = abs(fft(meanRun(isnr,ijitter,iampJitter,:)));
        end
    end
end

nframes = length(fftRun);
%from plotMeanFourierAmp.m
fftRun = fftRun / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
fftRun = fftRun(:,:,:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));



%%
trialsAmp = squeeze(std(rwdTrials,0,4));
ampMeanTrial = std(meanTrial,0,4);

%% 
linewidth = 0.1;
markersize = 10;
fontsize=9;

i=i+1; figure(i) ;clf
rows=4; cols=2;
linewidth = 1;
markersize = 10;
%HRF % derivative
subplot(rows,cols,1)
hrf = hrfModel(1:upsampleFactor:end);
hrf = zscore(hrf);
plot(TR*(0:length(hrf)-1),hrf,'k')
hold on
hrfDiff = diff(hrf);
plot(TR*(0:length(hrf)-2),abs(hrfDiff),'k--')
xlabel('time (sec)');
xlim([0 length(hrf)-1]);
% convolved hrf & derivative
subplot(rows,cols,2)

threeTrials = zeros(3*trialLength,1);
threeTrials(1:trialLength:end) = ones;
threeTrials = conv(threeTrials,hrf);
threeTrialsDiff = diff(threeTrials);
oneTrial = threeTrials(trialLength+1:2*trialLength);
oneTrial = zscore(oneTrial);
oneTrialDiff = threeTrialsDiff(trialLength+1:2*trialLength);

plot(TR*(0:length(oneTrial)-1),oneTrial,'k')
hold on
plot(TR*(0:length(oneTrialDiff)-1),abs(oneTrialDiff),'k--')
xlabel('time (sec)');
xlim([0 TR*trialLength]);
%SNR
subplot(rows,cols,cols+1)
ijitter=1; iampJitter=1;
for isnr=1:length(snr)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,cols+2)
for isnr=1:length(snr)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end


%amplitude variability
subplot(rows,cols,2*cols+1)
isnr=1; ijitter=1;
for iampJitter=1:length(ampJitter)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2*cols+2)
for iampJitter=1:length(ampJitter)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end


%temporal variability
subplot(rows,cols,3*cols+1)
isnr=1; iampJitter=1;
for ijitter=1:length(jitter)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,3*cols+2)
for ijitter=1:length(jitter)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end
%%
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    xlabel('time (sec)');
    if isubplot<3
        ylabel('response (arb units)');
    elseif mod(isubplot,2)==1 
        ylabel('response (arb units)');
    else
        ylabel('variability (std)');
    end
    if toZscore
        drawPublishAxis('xLabelOffset', -9/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',fontsize);
    else
        drawPublishAxis('xLabelOffset', -9/64,'yLabelOffset', -16/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',fontsize);
    end
    axis square
    legend off
end
set(gcf,'position',[10 10 9.7 21]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig6_hrfType' num2str(hrfType) '.pdf']);

%% FIGURE S6
%plot amplitude as function of noise level
i=i+1; figure(i) ;clf
rows=1;
cols=1;
%SNR
subplot(rows,cols,1)
ijitter=1; iampJitter=1;
baseAmp = ampMeanTrial(1,ijitter,iampJitter);
plot(squeeze(ampMeanTrial(:,ijitter,iampJitter)./baseAmp));
hold on
isnr=1; ijitter=1;
baseAmp = ampMeanTrial(isnr,ijitter,1);
plot(squeeze(ampMeanTrial(isnr,ijitter,:))./baseAmp);
isnr=1; iampJitter=1;
baseAmp = ampMeanTrial(isnr,1,iampJitter);
plot(squeeze(ampMeanTrial(isnr,:,iampJitter))./baseAmp);
xlabel('noise level');
ylabel('normalized amplitude (std)');
xlim([1 npoints]);
axis square

drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
        'yAxisMinMaxSetByTicks',0,'labelFontSize',fontsize,'xTick',[1 npoints],'xTickLabel',{'low','high'},'yTick',[0 1.2],...
        'yAxisMin',0,'yAxisMax',1.2);
    if hrfType==1
        h=legend({'independent noise','amplitude jitter','temporal jitter'},'Location','southwest');
    else
        legend off
    end
set(gcf,'position',[10 10 7 5]);
print('-painters','-dpdf',['figS6_simAmplitude_hrfType' num2str(hrfType) '.pdf']);
end
%%
function filteredTC = filterTC(tc)
filt = ones(1,length(tc));
filtNonzeros = [0,0,0.117503097415405,0.588887709492812,0.904365555167461,0.988891003461758,0.999355620179431];
filt(2:8) = filtNonzeros;
filt(end-6:end) = filtNonzeros;
filteredTC = real(ifft(fft(tc) .* filt' ))';
end