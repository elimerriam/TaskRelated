% figS2
% 
% associated with the following publication: Roth, ZN, Ryoo, M, and Merriam, EP (2020). 
% Task-related activity in human visual cortex. 
%
%   usage: figS2()
%   by: zvi roth
%   date: 9/3/2020
%   purpose: create figure S2, compare task-related activity to simulus
%   activity using inverted encoding pRF model
%


function[] = figS2()

dataFolder = '~/taskRelatedData/';%change this to the folder containing all subjects' fMRI data
loadData = 0;%1 = load saved data, or 0 = retrieve it anew
axRange = 20;%range of view of the visual field
xRes = 50;%resolution


dx=axRange/xRes;
if loadData
    load('voxVisualField.mat',...
        'dataFolder','subFolders','roiNames','goodSubs','roiNums','rwdScans','locScan','prfGroupName',...
        'matchGroup','matchScan',...
        'allRwdCo','allRwdPh','allRwdAmp','allEcc','allAngle','allAreas','numVox','allSubs');
else
    subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212',...
        '003220180105', '0034i20180209',  '003520180328', '0039i20180220', '004020180328','004120180320', ...
        '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
    roiNames = { 'leftBenson', 'rightBenson','Lviz_localizer','Rviz_localizer','lh_17Networks_16','rh_17Networks_16','leftCerebellarCortex','rightCerebellarCortex'};
    numSubs=length(subFolders);%number of subjects
    goodSubs = [1:3 5:length(subFolders)];%excluding 1 subject
    
    avgGroupName = 'Averages';
    roiNums = 1:4;%1:2
    rwdScans = 1:2;%first 2 scans are high-reward and low-reward
    locScan = 3;%third scan is localizer
    prfGroupName = 'templates';
    matchGroup = 'Averages';
    matchScan = rwdScans(1);
    cd(dataFolder);
    
    %initialize arrays
    for iRoiNum=1:length(roiNums)
        for rwd=1:3
            allRwdCo{iRoiNum,rwd} =[];
            allRwdPh{iRoiNum,rwd} =[];
            allRwdAmp{iRoiNum,rwd} =[];
        end
        allEcc{iRoiNum}=[];
        allAreas{iRoiNum}=[];
        allAngle{iRoiNum}=[];
        allSubs{iRoiNum} = [];
    end

    %get data, insert into arrays
    for iGoodSub=1:length(goodSubs)
        
        iSub = goodSubs(iGoodSub);
        cd(subFolders{iSub});
        v=newView;

        % load task-related data
        for iRoiNum=1:length(roiNums)
            iRoi = roiNums(iRoiNum);
            for iscan=1:2
                rwdCorrData = loadROIcoranalMatching(v, roiNames{iRoi}, rwdScans(iscan), avgGroupName, matchScan, matchGroup);
                %figure out what reward type it is
                s = viewGet(v, 'stimfile', iscan);
                if ~isfield(s{1}, 'stimulus')
                    s=s{1};
                end
                rwdType = s{1}.stimulus.rewardVal;
                if strcmp(rwdType, 'H')
                    rwd = 1;
                elseif strcmp(rwdType, 'L')
                    rwd = 2;
                end
                
                rwdCo{iGoodSub,iRoiNum,rwd} = rwdCorrData{1}.co;%coherence
                rwdPh{iGoodSub,iRoiNum,rwd} = rwdCorrData{1}.ph;%phase
                allRwdCo{iRoiNum,rwd} = [allRwdCo{iRoiNum,rwd}; rwdCorrData{1}.co];
                allRwdPh{iRoiNum,rwd} = [allRwdPh{iRoiNum,rwd}; rwdCorrData{1}.ph];
                allRwdAmp{iRoiNum,rwd} = [allRwdAmp{iRoiNum,rwd}; rwdCorrData{1}.amp];
            end
            locCorrData = loadROIcoranalMatching(v, roiNames{iRoi}, locScan, avgGroupName, matchScan, matchGroup);
            allRwdCo{iRoiNum,3} = [allRwdCo{iRoiNum,3}; locCorrData{1}.co];
            allRwdPh{iRoiNum,3} = [allRwdPh{iRoiNum,3}; locCorrData{1}.ph];
            allRwdAmp{iRoiNum,3} = [allRwdAmp{iRoiNum,3}; locCorrData{1}.amp];
        end
        %load benson prf data
        v = viewSet(v, 'curGroup', 'templates');
        templateGroup = viewGet(v,'curGroup');
        v = loadAnalysis(v, 'mrDispOverlayAnal/templateRet.mat');
        for iRoiNum=1:length(roiNums)
            iRoi = roiNums(iRoiNum);
            bensonData = loadROIbensonMatching(v,roiNames{iRoi},1,templateGroup,matchScan,matchGroup);
            eccen{iGoodSub,iRoiNum} = bensonData{1}.eccen;
            ang{iGoodSub,iRoiNum} = bensonData{1}.ang;%0 = top. 180 = bottom.
            areas{iGoodSub,iRoiNum} = bensonData{1}.areas;
            
            
            allEcc{iRoiNum} = [allEcc{iRoiNum}; bensonData{1}.eccen];
            allAngle{iRoiNum} = [allAngle{iRoiNum}; bensonData{1}.ang];
            allAreas{iRoiNum} = [allAreas{iRoiNum}; bensonData{1}.areas];
            
            numVox(iSub,iRoiNum) = length(bensonData{1}.eccen);
            allSubs{iRoiNum} = [allSubs{iRoiNum}; iSub*ones(size(bensonData{1}.eccen))];
        end
        
        deleteView(v);
        cd(dataFolder);
    end
    save('voxVisualField.mat',...
        'dataFolder','subFolders','roiNames','goodSubs','roiNums','rwdScans','locScan','prfGroupName',...
        'matchGroup','matchScan',...
        'allRwdCo','allRwdPh','allRwdAmp','allEcc','allAngle','allAreas','numVox','allSubs');
end


%now combine left and right ROIs
newRoiNum = 3;
roiNums(newRoiNum+1:end) = [];%keep only left and right EVC (benson atlas)
allAngle{newRoiNum} = [allAngle{1}; allAngle{2}];
allEcc{newRoiNum} = [allEcc{1}; allEcc{2}];
allAreas{newRoiNum} = [allAreas{1}; allAreas{2}];
for rwd=1:3
    allRwdCo{newRoiNum,rwd} = [allRwdCo{1,rwd}; allRwdCo{2,rwd}];
    allRwdPh{newRoiNum,rwd} = [allRwdPh{1,rwd}; allRwdPh{2,rwd}];
    allRwdAmp{newRoiNum,rwd} = [allRwdAmp{1,rwd}; allRwdAmp{2,rwd}];
end
roiNames{newRoiNum} = 'benson';



% get the Cartesian coordinates for each voxel
for iRoiNum=1:length(roiNums)
    polarAngle = allAngle{iRoiNum}*pi/180;
    if mod(iRoiNum,2)==1%left ROI
        polarAngle = mod(pi - polarAngle + (3*pi/2), 2*pi);
    else
        polarAngle = polarAngle+pi/2;
    end
    [allX{iRoiNum} allY{iRoiNum}] = pol2cart(polarAngle, allEcc{iRoiNum});
end
allX{newRoiNum} = [allX{1}; allX{2}];
allY{newRoiNum} = [allY{1}; allY{2}];

fontsize=11;
figure(1); clf;

myCols = circshift(hsv(64), [-16 1]);
myCols_amp = hot(64);
[gX,gY] = meshgrid(-axRange:dx:axRange,-axRange:dx:axRange);
rows=1;%length(roiNums);
cols=3;
for iRoiNum=3:length(roiNums)
    for rwd=1:3
        tic
        voxlist = find(allAreas{iRoiNum}>0);%voxels in v1,v2,v3

        gCo = zeros([size(gX) length(voxlist)]);
        gPh = zeros(size(gCo));
        gAmp = zeros(size(gCo));

        rfSizes = (allEcc{iRoiNum}.^1).*(allAreas{iRoiNum}.^0.7)/5;

        for iVox = 1:length(voxlist)
            thisvox = voxlist(iVox);
            rfSize = rfSizes(thisvox);
            p(1) = allRwdCo{iRoiNum,rwd}(thisvox);
            p(2) = allX{iRoiNum}(thisvox);
            p(3) = allY{iRoiNum}(thisvox);
            p(4) = rfSize;
            p(5) = rfSize;
            p(6) = 0;
            p(7) = 0;
            gCo(:,:,iVox) = eventRelatedGauss(p,gX,gY);
            gPh(:,:,iVox) = ones(size(gX)) .* allRwdPh{iRoiNum,rwd}(thisvox);
            gAmp(:,:,iVox) = gCo(:,:,iVox)*allRwdAmp{iRoiNum,rwd}(thisvox)/allRwdCo{iRoiNum,rwd}(thisvox);
        end
        
        
        gAmp(isnan(gAmp)) = 0;
        gCo(isnan(gCo)) = 0;
        gPh(isnan(gPh)) = 0;
        
        %coherence and phase
        gCpx = gCo .* (cos(gPh) + 1i*sin(gPh));
        ph = angle(mean(gCpx,3));
        co = abs(mean(gCpx,3));
        co = flipud(co)./.01;
        ph = flipud(ph);
        minCo(iRoiNum,rwd) = min(co(:));
        maxCo(iRoiNum,rwd) = max(co(:));
        allCo{iRoiNum,rwd} = co;
        
        figure(1)
        subplot(rows,cols,rwd);
        imshow(ph, 'colormap', myCols_amp, 'InitialMag', 'fit')
        colormap(myCols)
        xticks([1 xRes 2*xRes]);
        xticklabels([-axRange 0 axRange]);
        yticks([1 xRes 2*xRes]);
        yticklabels([axRange 0 -axRange]);
        overlay = cat(3, zeros(size(co)),zeros(size(co)), zeros(size(co)));
        hold on
        hCo(iRoiNum,rwd) = imshow(overlay);
        hold off
        %         set(h, 'AlphaData', 1-co);
        axis on
        caxis([-pi pi]);
        
        
        toc
    end
end

%%

figure(1)
set(gcf,'position',[100 100 500 500]);
for iRoiNum=3:length(roiNums)
    subplot(rows,cols,1)
    ylabel('y (deg)','fontangle','italic','fontsize',fontsize);
end
iRoiNum=3;
for rwd=1:3
    subplot(rows,cols,rwd)
    xlabel('x (deg)','fontangle','italic','fontsize',fontsize);
end


figure(1)
for iRoiNum=3:length(roiNums)
    for rwd=1:3
        subplot(rows,cols,rwd)
        c= allCo{iRoiNum,rwd};
        c = c-min(minCo(:));
        c = c./(max(maxCo(:)) - min(minCo(:)));
        set(hCo(iRoiNum,rwd), 'AlphaData', 1-sqrt(c));
    end
end


figure(1)
print('-painters','-dpdf',['figS2_' num2str(axRange) '.pdf']);
