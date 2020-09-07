% saveControlData
% 
% associated with the following publication: Roth ZN, Ryoo M, and Merriam EP (2020). 
% Task-related activity in human visual cortex. 
%
%
%   usage: saveControlData()
%   by: zvi roth
%   date: 9/4/2020
%   purpose: save data from control experiment, with no stimulus
%
%
function[] = saveControlData()

toZscore=1;%0 or 1
ConcatProj = 0;
curFolder = pwd;
dataFolder = '/Volumes/TaskDrive/controlExp/';
saveFolder = dataFolder;
subFolders = {'000920170810', '001020170815','001520170809','002120170814','002220170817','002420170828'};
numSubs=length(subFolders);
roiNames = { 'leftBenson', 'rightBenson'};

cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
trialLength=10;
rows=2;
cols=ceil(numSubs/rows);
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2],[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};

for iSub = 1:numSubs
    iSub
    cd(subFolders{iSub});
    v=newView;
    
    % Find the correct Concatenation group
    if ConcatProj
        concatGroupNum = viewGet(v,'groupNum','ConcatenationProj'); %concatenation of multiple sessions
    else
        concatGroupNum = viewGet(v,'groupNum','Concatenation');
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);
    iScan=1;
    concatInfo{iSub} = viewGet(v, 'concatInfo', iScan);
    for iRoi = 1:length(roiNames)
        roiTC{iSub,iRoi} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
        roiTC{iSub,iRoi}.tSeries = 100*(roiTC{iSub,iRoi}.tSeries-1);%percent signal change
        voxNum(iSub,iRoi) = size(roiTC{iSub,iRoi}.tSeries,1);
        
        if toZscore
            roiTC{iSub,iRoi}.tSeries = zscoreConcat(roiTC{iSub,iRoi}.tSeries, concatInfo{iSub});
        end
        roiMeanTseries{iSub,iRoi}(:) = nanmean(roiTC{iSub,iRoi}.tSeries);%mean across voxels
        subTrialResponse{iSub,iRoi} = reshape(roiMeanTseries{iSub,iRoi}(:), trialLength, length(roiMeanTseries{iSub,iRoi}(:))/trialLength);
        %average per subject
        subResponse(iSub,iRoi,:) = mean(subTrialResponse{iSub,iRoi},2);
        subStd(iSub,iRoi,:) = std(subTrialResponse{iSub,iRoi},0,2);
        %single voxel responses
        voxTrials{iSub,iRoi} = reshape(roiTC{iSub,iRoi}.tSeries, size(roiTC{iSub,iRoi}.tSeries,1), trialLength, length(roiMeanTseries{iSub,iRoi}(:))/trialLength);
        meanVoxTrial{iSub,iRoi} = squeeze(nanmean(voxTrials{iSub,iRoi},3));
        
    end

    %load benson eccentricity maps
    v = viewSet(v, 'curGroup', 'templates');
    templateGroup = viewGet(v,'curGroup');
    v = loadAnalysis(v, 'mrDispOverlayAnal/templateRet.mat');
    for iRoi = 1:length(roiNames)
        bensonData = loadROIbensonMatching(v,roiNames{iRoi},1,templateGroup,1,concatGroupNum);
        eccen{iSub,iRoi} = bensonData{1}.eccen;
        ang{iSub,iRoi} = bensonData{1}.ang;
        areas{iSub,iRoi} = bensonData{1}.areas;
    end

    deleteView(v);
    cd ..
end

%%

zScoreString = '';
if toZscore
    zScoreString = '_zscored';
end

ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end

%%
save([saveFolder 'das_'  zScoreString  ConcatProjStr  '.mat'],   'subResponse', 'roiMeanTseries', ...
    'roiTC', ...
    'subFolders', 'roiNames','subTrialResponse',...
    'eccen','ang','areas','trialLength',...
    'subStd',...
    'voxTrials','meanVoxTrial',...
    'voxNum','concatInfo');
end

%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
    newConcat = [newConcat zscore(thisRun,0,2)];
end
end