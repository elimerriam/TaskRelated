% figS3
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
%   uses file created by saveTaskEyeData.m
%
%   usage: figS3()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: create Figure S3, analysis of eye movements in eye tracking
%   experiment
%
%

function[] = figS3()

dataFolder = '';
fontsize=9;

ampMax=10;
ampMin = 0.1;
durationMax = 100;
peakVelMax = 1000;
peakVelMin = 10;

sacTime = 1.00;
            
load([dataFolder 'eyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'numSubs', 'pupilRun', 'numRuns',...
    'allMicrosacc','imageWidth','screenWidth','imageHeight','screenHeight','allSac','allFix','startTimes','allMgl');

plotColors = {'b','k','r','g','y','c','m','b','k','r','g','y','c','m','b','k','r','g','y','c','m'};
numSubs = length(subFolders);

rows=2;
cols=numSubs;
[targetX targetY] = pol2cart(2*pi - pi/4, 5);
colorLength = 64;
myColor = cool(colorLength);
fov = 8;

allTh=[];
allR = [];
allSacStart=[];
allPeakVel = [];
allSacDur = [];

colormap('cool');

[targetX targetY] = pol2cart(2*pi - pi/4, 5);
colorLength = 64;
myColor = cool(colorLength);
fov = 8;

allX = [];
allY = [];
allFixDur = [];
allTrialTime = [];
allSubNum=[];

for iSub=1:numSubs
    
    for rwd=1:2
        for r=1:numRuns(iSub,rwd)
            fixStartTime = allFix{iSub,rwd,r}.startTime;%in ms
            fixEndTime = allFix{iSub,rwd,r}.endTime;
            fixH = allFix{iSub,rwd,r}.aveH;
            fixV = allFix{iSub,rwd,r}.aveV;
            fixH = imageWidth(iSub,rwd,r)*(fixH-screenWidth(iSub,rwd,r)/2)/screenWidth(iSub,rwd,r);
            fixV = imageHeight(iSub,rwd,r)*(fixV-screenHeight(iSub,rwd,r)/2)/screenHeight(iSub,rwd,r);
            fixDuration = fixEndTime - fixStartTime;

            %compute saccade time from trial start
            fixStartTime = fixStartTime-fixStartTime(1);
            fixStartTime = fixStartTime./1000;
            trialStartTime = startTimes{iSub,rwd,r};
            trialLength = diff(trialStartTime);
            trialStartTime(end+1) = trialStartTime(end) + min(trialLength);%last 'trial' is actually the 15 sec message screen
            for ifix=1:length(fixStartTime)
                temp = fixStartTime(ifix) - trialStartTime;
                %assign this saccade to a trial number
                if max(temp)>=0
                    [fixTrialTime{iSub,rwd,r}(ifix) fixTrialNum{iSub,rwd,r}(ifix)] = min(temp(temp>=0));
                else
                    fixTrialTime{iSub,rwd,r}(ifix)=0;
                    fixTrialNum{iSub,rwd,r}(ifix)=0;
                end
            end
            fixColor = fixTrialTime{iSub,rwd,r};
            maxFixTime(iSub,rwd,r) = max(fixColor);
            %only plot fixations that are before the end of the experiment
            beforeEnd = fixTrialNum{iSub,rwd,r}<length(trialStartTime);
            afterEnd = fixTrialNum{iSub,rwd,r}>=length(trialStartTime);
            
            allX = [allX fixH(beforeEnd)];
            allY = [allY fixV(beforeEnd)];
            allFixDur = [allFixDur fixDuration(beforeEnd)];
            allTrialTime = [allTrialTime fixTrialTime{iSub,rwd,r}(beforeEnd)];

            sacStart = allSac{iSub,rwd,r}.startTime;
            sacEnd = allSac{iSub,rwd,r}.endTime;
            sacStart = sacStart./1000;
            sacStartTime = sacStart - sacStart(1);
            
            clear sacTrialTime
            for isac=1:length(sacStartTime)
                temp = sacStartTime(isac) - trialStartTime;
                %assign this saccade to a trial number
                if max(temp)>=0
                    [sacTrialTime{iSub,rwd,r}(isac) sacTrialNum{iSub,rwd,r}(isac)] = min(temp(temp>=0));
                else
                    sacTrialTime{iSub,rwd,r}(isac)=0;
                    sacTrialNum{iSub,rwd,r}(isac)=0;
                end
            end
            peakVel = allSac{iSub,rwd,r}.peakVel;

            sacDur = allSac{iSub,rwd,r}.endTime - allSac{iSub,rwd,r}.startTime;%in ms
            x = allSac{iSub,rwd,r}.endH - allSac{iSub,rwd,r}.startH;%in pixels
            y = allSac{iSub,rwd,r}.endV - allSac{iSub,rwd,r}.startV;
            x = imageWidth(iSub,rwd,r)*x/screenWidth(iSub,rwd,r);
            y = imageHeight(iSub,rwd,r)*y/screenHeight(iSub,rwd,r);

            [theta rho] = cart2pol(x,y);%in pixels

            %good saccades are before the end of the experiment, exclude
            %the first saccade from each run, and aren't incredibly large
            %or fast
            goodSac = sacTrialNum{iSub,rwd,r}>1 & sacTrialNum{iSub,rwd,r}<length(trialStartTime) & sacTrialTime{iSub,rwd,r}>0 ...
                &  sacDur<durationMax & peakVel<peakVelMax & peakVel>peakVelMin;
            
            allTh = [allTh theta(goodSac)];
            allR = [allR rho(goodSac)];%in pixels
            allSacStart = [allSacStart sacTrialTime{iSub,rwd,r}(goodSac)];
            allPeakVel = [allPeakVel peakVel(goodSac)];
            allSacDur = [allSacDur sacDur(goodSac)]; %in ms
            allSubNum = [allSubNum iSub*ones(size(sacDur(goodSac)))];
            
        end
        axis square
    end
end
%%
figure(1)
clf
rows=2;
cols=2;
isubplots={1,2,3,4};
% isubplots={1:3,5:8,9:11, 13:16};

fontsize=11;
markersize = 1;
edgeAlpha = 0.07;% 0.15;
markerColor = [0 0 0];
subplot(rows,cols,isubplots{1})

s=scatter(allR(allSacStart<sacTime), allPeakVel(allSacStart<sacTime),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([peakVelMin peakVelMax])
xlim([ampMin ampMax])
xticks([0.1 1 10]);
xticklabels([0.1 1 10]);
yticks([10 100 1000]);
yticklabels([10 100 1000]);
grid on
axis square
xlabel('amplitude (deg)','fontangle','italic','fontsize',fontsize);
ylabel('peak velocity (deg/s)','fontangle','italic','fontsize',fontsize);
box on


subplot(rows,cols,isubplots{4})
h=polarhistogram(allTh(allSacStart<sacTime),40,'FaceColor',[0.4 0.4 0.4], 'FaceAlpha',0.5,'normalization','count');
subplot(rows,cols,isubplots{2})
h=histogram(allR(allSacStart<sacTime));
h.FaceColor = [0.5 0.5 0.5];
xlim([0 ampMax])
% axis tight
axis square
xlabel('amplitude (deg)','fontangle','italic','fontsize',fontsize);
ylabel('# saccades','fontangle','italic','fontsize',fontsize);

% scatter plot of all saccades
subplot(rows,cols,isubplots{3})
[x y] = pol2cart(allTh(allSacStart<sacTime),allR(allSacStart<sacTime));
s=scatter(x, y,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
colormap cool
xlimit = 5;
xlim([-xlimit xlimit]);
ylim([-xlimit xlimit]);
hline(0);
vline(0);
axis square
xlabel('x displacement (deg)','fontangle','italic','fontsize',fontsize);
ylabel('y displacement (deg)','fontangle','italic','fontsize',fontsize);
box on

set(gcf,'position',[120 80 450 450]);
print('-painters','-dpdf',['figS3.pdf']);

