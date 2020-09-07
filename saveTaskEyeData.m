% saveTaskEyeData
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
%   creates file used by figS3.m
%
%   usage: saveTaskEyeData()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: save eye tracking data from eye tracking experiment
%
%

function[] = saveTaskEyeData()

VFAC = 8;
MINDUR = 7;

dataFolder = '/Users/rothzn/min/';
saveFolder = dataFolder;
samplerate=500;%in ms
subFolders = {'s000520180126','s000820171116','s001620171103','s002220171122','s003220180105','s003920180220',...
    's004020180221','s004120180223','s004320180306','s004520180308','s004620180323','s004720180326',...
    's004920180404'};

numSubs = length(subFolders);
curFolder = pwd;
onlyCorrect=0;
cd(dataFolder);
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};

%smooth out the microsaccade data
L = 100;
filter = ones(L,1);
filter = filter./sum(filter);
allSac = cell(2,1);
for iSub = 1:numSubs
    cd(subFolders{iSub});
    
    stimfiles = dir('*.mat');
    listRwdRuns{iSub,1} = [];
    listRwdRuns{iSub,2} = [];
    for iRun=1:length(stimfiles)
        stimfile=stimfiles(iRun).name;
        allStimfiles{iRun}=load(stimfile);
        
        expName{iSub,iRun} = allStimfiles{iRun}.task{1}.taskFilename;

        rwdType = allStimfiles{iRun}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            keyboard
        end
        concatRwdTypeNum(iSub,iRun) = rwd;
        listRwdRuns{iSub,rwd} = [listRwdRuns{iSub,rwd}; iRun];
    end
    for rwd=1:2
        clear s
        for r=1:length(listRwdRuns{iSub,rwd})%going through all stimfiles for this reward level
            iRun = listRwdRuns{iSub,rwd}(r);
            s{r} = allStimfiles{iRun};
            imageWidth(iSub,rwd,r) = s{r}.myscreen.imageWidth;
            screenWidth(iSub,rwd,r) = s{r}.myscreen.screenWidth;
            imageHeight(iSub,rwd,r) = s{r}.myscreen.imageHeight;
            screenHeight(iSub,rwd,r) = s{r}.myscreen.screenHeight;

            trialCorrectness{iSub,rwd}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
            trialResponse{iSub,rwd}(r,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response
            trialRT{iSub,rwd}(r,:) = [s{r}.stimulus.trialRT NaN(1,17-size(s{r}.stimulus.trialRT,2))];
            propCorrect{iSub,rwd}(r) = s{r}.stimulus.percentCorrect;
            stairThresh{iSub,rwd}(r,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
                
            stimfile = stimfiles(iRun).name;
            e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=10');%e{r}.saccades, e{r}.fixations
            allSac{iSub,rwd,r} = e{r}.saccades;%
            allFix{iSub,rwd,r} = e{r}.fixations;%from mglPrivateEyelinkEDFRead.m
            allMgl{iSub,rwd,r} = e{r}.mgl;
            startTimes{iSub,rwd,r} = e{r}.complete.startTime;
            eyetrackerTime{iSub,rwd,r} = e{r}.complete.time;
            runEyeSize{rwd}(r,:) = size(e{r}.eye.pupil); 
            startTimes{iSub,rwd,r} = e{r}.trialTime-e{r}.trialTime(1);
        end
        trialLengthEye(rwd) = max(runEyeSize{rwd}(:,2));
        numTrials(rwd) = sum(runEyeSize{rwd}(:,1));
        sacRwd{iSub,rwd} = zeros(numTrials(rwd), trialLengthEye(rwd));
        trialCounter=0;
        numRuns(iSub,rwd) = length(s);
        for r=1:length(s)
            
            eyePos = [e{r}.complete.xPos' e{r}.complete.yPos'];% time X 2, x and y
            if sum(~isnan(eyePos(:)))~=0
                eyevel{r} = vecvel(eyePos,samplerate, 3);
                sacOutput = microsacc(eyePos, eyevel{r}, VFAC, MINDUR);
                
                allMicrosacc{iSub,rwd,r} = sacOutput;
                
                sacOnset = sacOutput(:,1);               
                sacRun{iSub,rwd,r} = sacOnset;%saccade onset
                sacRunTC = zeros(size(eyePos,1),1);
                sacRunTC(sacOnset) = ones;
                sacRunSmooth{iSub,rwd,r} = conv(sacRunTC, filter, 'same');
                pupilRun{iSub,rwd,r} = e{r}.complete.pupil;%e{r}.eye.pupil is longer because of zero padding trials;

                for itrial=1:length(e{r}.complete.startTime)
                    trialCounter=trialCounter+1;
                    sac{r,itrial} = sacOnset(sacOnset >= e{r}.complete.startTime(itrial) & sacOnset<=e{r}.complete.endTime(itrial));%saccade onsets for this trial
                    sac{r,itrial} = sac{r,itrial} - e{r}.complete.startTime(itrial) + 1;%time relative to trial
                    sacTrial = zeros(1,trialLengthEye(rwd));
                    sacTrial(sac{r,itrial}) = ones;
                    sacRwd{iSub,rwd}(trialCounter,:) = sacTrial;
                    
                    sacRunSmoothTrial{r,itrial} = sacRunSmooth{iSub,rwd,r}(e{r}.complete.startTime(itrial):e{r}.complete.endTime(itrial));
                    sacTrial = NaN(1,trialLengthEye(rwd));
                    sacTrial(1:length(sacRunSmoothTrial{r,itrial})) = sacRunSmoothTrial{r,itrial};
                    sacRunSmoothTrialZeroFilled{iSub,rwd}(trialCounter,:) = sacTrial;
                end
            end
        end

        binnedSac{iSub,rwd} = sum(sacRwd{iSub,rwd});
        rwdPupil{iSub,rwd} = NaN(numTrials(rwd), trialLengthEye(rwd));
        nullTrials{iSub,rwd} = NaN(numTrials(rwd),1);
        trialCounter=0;
        for r=1:length(s)         
            rwdPupil{iSub,rwd}(trialCounter+1:trialCounter+runEyeSize{rwd}(r,1), 1:runEyeSize{rwd}(r,2)) = e{r}.eye.pupil;
            trialCounter = trialCounter + runEyeSize{rwd}(r,1);
        end
        %plot mean pupil size
        meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';  
        
        smoothSac{iSub,rwd} = nanmean(sacRunSmoothTrialZeroFilled{iSub,rwd});
    end
    
    %% return to home directory
    cd(dataFolder);
end

%%
save([saveFolder 'eyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil',...
    'sacRun','pupilRun', 'numRuns','eyetrackerTime',...
    'listRwdRuns', 'concatRwdTypeNum','expName',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L',...
    'allMicrosacc','imageWidth','screenWidth','imageHeight','screenHeight','allSac','allFix','startTimes','allMgl',...
    'VFAC','MINDUR');

