%%% -----------------------------------------------------------------------
%%% This script will generate plots associated to a given state vector
%%% (figure 2 d-g)
%%% -----------------------------------------------------------------------
%%% Requires `main.m` `locaNMF` and `computeStateVectors` to be run first

%% Load project data

clearvars;
addpath(genpath(pwd));
projectFolder = ['..' filesep 'widefieldChoice' filesep 'data' filesep 'excMice'];
load([projectFolder filesep 'projectData.mat']);


%% Preload animal data

dsetAnimal = cell(Ndatasets, 1);
dataLocaAnimal = cell(Ndatasets, 1);
dataSVDAnimal = cell(Ndatasets, 1);
dataBehaviorAnimal = cell(Ndatasets, 1);
baseCompAnimal = cell(Ndatasets, 1);
for itAnimal = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{itAnimal});
  dset.rootFolder = projectFolder; % Overwrite for machine swaps
  dataLoca = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_locaNMF_hdf5_output.mat', dset.name)));
  dataSVD = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_SVD_minimal.mat', dset.name)));
  %dataSVDfull = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_SVD.mat', dset.name)));
  dataBehavior = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_jointSessionsBehavior.mat', dset.name)));
  baseComp = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_baseComponents.mat', dset.name)));
  dsetAnimal{itAnimal} = dset;
  dataLocaAnimal{itAnimal} = dataLoca;
  dataSVDAnimal{itAnimal} = dataSVD;
  dataBehaviorAnimal{itAnimal} = dataBehavior;
  baseCompAnimal{itAnimal} = baseComp;
end
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
regionOrder = [1 2 3 4 6 8 9 5 7 10];
regionColor = [2 3 3 4 4 4 5 6 6 7];
regionPosition = [3 5 6 8 9 10 12 14 16];
areaGroupIDs = {{[2 ,3], [4, 5, 6, 7, 8], [9]}, {[2, 3], [4, 6, 8], [9]}, {[1], [10], [5, 7]}};
areaGroupNames = {{'dorsal', 'parietalSS', 'ventral'}, {'dorsal', 'parietal', 'ventral'}, {'V1', 'RS', 'SS'}};
areaGroupNames2 = {{'PM,AM', 'A,AL,RL,SSt,SSb', 'L'}, {'PM,AM', 'A,RL,AL', 'L'}, {'V1', 'RS', 'SSt,SSb'}};

%% Compute everything
% Paramters
singleTimePoint = 0.3; % With respect to the movement time
kfoldN = 5; % For the cross-valiation
wSize = 3;
zeroComps = true;
reverseTrainTest = false;
zscoreComp = true;
movWindow = 31:150;
noMovWindow = 15:30;
%%% PupA definitions
presFrames = 1:30;
OLFrames = 31:75;
vectorModes = {'mov','no mov'}; % Vector will be V1 - V2

% Initialization
projMeanFull1 = cell(Ndatasets, 1);
projSEMFull1 = cell(Ndatasets, 1);
projMeanFull2 = cell(Ndatasets, 1);
projSEMFull2 = cell(Ndatasets, 1);
projDprimeFull = cell(Ndatasets, 1);
projDprimeCompFull = cell(Ndatasets, 1);
wheelposFull1 = cell(Ndatasets, 1);
wheelposFull2 = cell(Ndatasets, 1);
pupilPosFull1 = cell(Ndatasets, 1);
pupilPosFull2 = cell(Ndatasets, 1);
projDprimeTimeCompFull = cell(Ndatasets, 1);
projDprimeTimeCompPerAreaFull = cell(Ndatasets, 1);
sVecFull = cell(Ndatasets, 1);
sVecMatFull = cell(Ndatasets, 1);
regionIdxFull = cell(Ndatasets, 1);
trialSplits = zeros(Ndatasets, 2);

singleProjMeanFull1 = cell(Ndatasets, 1);
singleProjSEMFull1 = cell(Ndatasets, 1);
singleProjMeanFull2 = cell(Ndatasets, 1);
singleProjSEMFull2 = cell(Ndatasets, 1);
singleProjDprimeFull = cell(Ndatasets, 1);
singleProjDprimeTimeCompPerAreaFull = cell(Ndatasets, 1);
singleProjDprimeTimeCompPerAreaGroupFull = cell(Ndatasets, 1);
projDprimeTimeCompPerAreaGroupFull = cell(Ndatasets, 1);

for itAnimal = 1:Ndatasets
  dset = dsetAnimal{itAnimal};
  dataLoca = dataLocaAnimal{itAnimal};
  dataSVD = dataSVDAnimal{itAnimal};
  %dataSVDfull = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_SVD.mat', dset.name)));
  dataBehavior = dataBehaviorAnimal{itAnimal};
  baseComp = baseCompAnimal{itAnimal};
  slog.info(dset.name);
  regionMasks = baseComp.a_dict.region_mats;
  % Because there is no SSt in dataset 6 (A15352)
  if((length(unique(dset.data.areasAllenLocaNMF.*dataSVD.mask))-1) ~= length(regionNames))
    slog.warning('Weird number of regions on dset %s', dset.name);
    dataLoca.region_ranks = [dataLoca.region_ranks(1:5) {int64(0)} dataLoca.region_ranks(6:end)];
    regionMasks = [regionMasks(1:4, :); zeros(1, size(regionMasks, 2)); regionMasks(5:end, :)];
  end

  regionCompSizes = cell2mat(dataLoca.region_ranks(2:end));
  regionIdx = [];
  for it = 1:length(regionCompSizes)
    regionIdx = [regionIdx, ones(1, regionCompSizes(it))*it];
  end

  v1 = var(dataLoca.compData);
  vexp  = zeros(length(v1), 1);
  for it0 = 1:length(v1)
    vexp(it0) = sum(dataLoca.compXY(it0, :).^2*v1(it0));
  end
  if(zscoreComp)
    locaZS = zscore(dataLoca.C, 0, 2);
  end
  %%% Reshape loca
  trialID = dataBehavior.trialID;
  locaReshaped = nan(length(dataBehavior.time), length(trialID), size(dataLoca.compData, 2));
  %locaReshaped = zscore(locaReshaped, 0, 3);
  for it = 1:length(trialID)
    curTrial = trialID(it);
    validIdx = find(dataSVD.trialID == curTrial);
    frames = find(dataSVD.frames(validIdx));
    %locaReshaped(frames, it, :) = dataLoca.compData(validIdx, :);
    if(zscoreComp)
      locaReshaped(frames, it, :) = locaZS(validIdx, :);  
    else
      locaReshaped(frames, it, :) = dataLoca.C(validIdx, :);  
    end
  end
  
  %%% Generate the state vectors for the stimulus
  vectorModesFrames = {find(dataSVD.time >= 0 & dataSVD.time < 2.5), find(dataSVD.time >= 0 & dataSVD.time < 2.5)}; % Frames used to average the modes
  DprimeFrames = find(dataSVD.time >= 1 & dataSVD.time < 1.5); % Frames used to average the projection and compute d prime
  
  leftOri = dset.data.eventTimesCondition.orientation(dataBehavior.trialID, 1);
  rightOri = dset.data.eventTimesCondition.orientation(dataBehavior.trialID, 2);
  sigPupA = max(dataBehavior.pupilArea(OLFrames, :)-mean(dataBehavior.pupilArea(presFrames, :)))'; % Max OL - mean prestim

  trialDiff = abs(abs(leftOri)-abs(rightOri));

  choices = dset.data.trialsSummary.choicedir(dataBehavior.trialID);
  choices2 = nan(size(dataBehavior.wheelMovementsSz, 2), 1);
  for it1 = 1:length(choices2)
    idx = find(~isnan(dataBehavior.wheelMovementsSz(:, it1)));
    firstIdx = idx(find(dataSVD.time(idx) > dataSVD.time(movWindow(1)) & dataSVD.time(idx) <= dataSVD.time(movWindow(end)), 1, 'first')); % Looking for the 1st mov on the window
    lastIdxNoMov = idx(find(dataSVD.time(idx) > dataSVD.time(noMovWindow(1)) & dataSVD.time(idx) <= dataSVD.time(noMovWindow(end)), 1, 'last')); % Looking for the last mov on the no window
    if(dataSVD.time(firstIdx) - dataSVD.time(lastIdxNoMov) < 0.5) % At least half a sec separation
      choices2(it1) = 2;
      continue;
    end
    if(sum(~isnan(dataBehavior.saccades(firstIdx-15:firstIdx, it1))) > 0) % Half a sec without saccades before movement
      choices2(it1) = 3;
      continue;
    end
    idxSac = find(~isnan(dataBehavior.saccades(:, it1)));
    interferingSacc = idxSac(find(dataSVD.time(idxSac) > 1 & dataSVD.time(idxSac) <= 1.5, 1, 'first'));  % Half a sec without saccades after movement
    if(~isempty(interferingSacc))
      choices2(it1) = 3;
      continue;
    end
    if(dataBehavior.wheelMovementsSz(firstIdx, it1) > 0)
      choices2(it1) = -1;
    elseif(dataBehavior.wheelMovementsSz(firstIdx, it1) < 0)
      choices2(it1) = 1;
    end
  end
  outcome = dset.data.trialsSummary.outcome(dataBehavior.trialID);

  selTrial1 = find(abs(choices2) == 1);
  selTrial2 = find(abs(choices2) ~= 1);
  
  movTime1 = zeros(size(selTrial1));
  movTime2 = zeros(size(selTrial2));
  % Let's do 1 sec before and after
  eventTimeData1 = nan(60, length(movTime1), size(locaReshaped, 3));
  eventTimeDataWheel1 = nan(60, length(movTime1));
  eventTimeDataPupilPos1 = nan(60, 2, length(movTime1));
  %dt = dataSVD.time(2)-dataSVD.time(1);
  eventTimeDataT = dataSVD.time(1:60)-1+dataSVD.time(2);
  for it0 = 1:length(selTrial1)
    movTime1(it0) = movWindow(1)-1+find(~isnan(dataBehavior.wheelMovementsSz(movWindow, selTrial1(it0))), 1, 'first');
    %movTime1(it0) = 60-1+find(dataBehavior.wheelMovements(movWindow, selTrial1(it0)) >= 1, 1, 'first');
    firstFrame = max(movTime1(it0)-29, 1);
    lastFrame = min(movTime1(it0)+30, size(dataBehavior.wheelMovementsSz, 1));
    eventTimeData1((firstFrame:lastFrame)-movTime1(it0)+30, it0, :) = locaReshaped(firstFrame:lastFrame, selTrial1(it0), :);
    eventTimeDataWheel1((firstFrame:lastFrame)-movTime1(it0)+30, it0) = dataBehavior.wheelP(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.wheelP(firstFrame, selTrial1(it0));
    eventTimeDataPupilPos1((firstFrame:lastFrame)-movTime1(it0)+30, 1, it0) = dataBehavior.pupilX(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.pupilX(firstFrame, selTrial1(it0));
    eventTimeDataPupilPos1((firstFrame:lastFrame)-movTime1(it0)+30, 2, it0) = dataBehavior.pupilY(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.pupilY(firstFrame, selTrial1(it0));
  end
  eventTimeData2 = nan(60, length(movTime2), size(locaReshaped, 3));
  eventTimeDataWheel2 = nan(60, length(movTime2));
  eventTimeDataPupilPos2 = nan(60, 2, length(movTime2));
  %dt = dataSVD.time(2)-dataSVD.time(1);
  %eventTimeDataT = dataSVD.time(1:60)-1+dataSVD.time(2);
  for it0 = 1:length(selTrial2)
    %movTime2(it0) = 60-1+find(dataBehavior.wheelMovements(movWindow, selTrial2(it0)) >= 1, 1, 'first');
    %movTime2(it0) = movWindow(1)-1+find(~isnan(dataBehavior.wheelMovementsSz(movWindow, selTrial2(it0))), 1, 'first');
    movTime2(it0) = movWindow(1)-1+randi(length(movWindow)); % Center these trials randomly on the same interval
    firstFrame = max(movTime2(it0)-29, 1);
    lastFrame = min(movTime2(it0)+30, size(dataBehavior.wheelMovementsSz, 1));
    eventTimeData2((firstFrame:lastFrame)-movTime2(it0)+30, it0, :) = locaReshaped(firstFrame:lastFrame, selTrial2(it0), :);
    eventTimeDataWheel2((firstFrame:lastFrame)-movTime2(it0)+30, it0) = dataBehavior.wheelP(firstFrame:lastFrame, selTrial2(it0))-dataBehavior.wheelP(firstFrame, selTrial2(it0));
    eventTimeDataPupilPos2((firstFrame:lastFrame)-movTime2(it0)+30, 1, it0) = dataBehavior.pupilX(firstFrame:lastFrame, selTrial2(it0))-dataBehavior.pupilX(firstFrame, selTrial2(it0));
    eventTimeDataPupilPos2((firstFrame:lastFrame)-movTime2(it0)+30, 2, it0) = dataBehavior.pupilY(firstFrame:lastFrame, selTrial2(it0))-dataBehavior.pupilY(firstFrame, selTrial2(it0));
  end
  
  [~, singleTimePointFrame] = min(abs(eventTimeDataT-singleTimePoint));

  projMean1 = zeros(size(eventTimeDataT, 1), kfoldN);
  projSEM1 = zeros(size(eventTimeDataT, 1), kfoldN);
  projMean2 = zeros(size(eventTimeDataT, 1), kfoldN);
  projSEM2 = zeros(size(eventTimeDataT, 1), kfoldN);
  projDprime = zeros(size(eventTimeDataT, 1), kfoldN);
  projDprimeComp = zeros(size(locaReshaped, 3), kfoldN); % Estimate of the D prime for each component, i.e.,  weight of SV
  projDprimeTimeComp = zeros(size(eventTimeDataT, 1), size(locaReshaped, 3), kfoldN); % Estimate of the D prime for each component, i.e.,  weight of SV
  projDprimeTimeCompPerArea = zeros(size(eventTimeDataT, 1), length(regionNames), kfoldN); % Estimate of the D prime for each component, i.e.,  weight of SV
  
  partitionSplit1 = cvpartition(length(selTrial1), 'KFold', kfoldN);
  partitionSplit2 = cvpartition(length(selTrial2), 'KFold', kfoldN);
  
  sVecMat = zeros(size(eventTimeDataT, 1)*(size(eventTimeDataT, 1)-1)/2, kfoldN);
  sVecFold = zeros(size(eventTimeDataT, 1), size(locaReshaped, 3), kfoldN);
  
  singleProjMean1 = zeros(size(eventTimeDataT, 1), kfoldN);
  singleProjSEM1 = zeros(size(eventTimeDataT, 1), kfoldN);
  singleProjMean2 = zeros(size(eventTimeDataT, 1), kfoldN);
  singleProjSEM2 = zeros(size(eventTimeDataT, 1), kfoldN);
  singleProjDprime = zeros(size(eventTimeDataT, 1), kfoldN);
  singleProjDprimeTimeComp = zeros(size(eventTimeDataT, 1), size(locaReshaped, 3), kfoldN);
  singleProjDprimeTimeCompPerArea = zeros(size(eventTimeDataT, 1), length(regionNames), kfoldN);
  singleProjDprimeTimeCompPerAreaGroup = zeros(size(vectorModesFrames{1}, 1), 3, 3, kfoldN); % subareas / region definitions
  projDprimeTimeCompPerAreaGroup = zeros(size(vectorModesFrames{1}, 1), 3, 3, kfoldN); % subareas / region definitions

  for curFold = 1:kfoldN
    if(~reverseTrainTest)
      trainSet1 = training(partitionSplit1, curFold);
      testSet1 = test(partitionSplit1, curFold);
      trainSet2 = training(partitionSplit2, curFold);
      testSet2 = test(partitionSplit2, curFold);
    else
      trainSet1 = test(partitionSplit1, curFold);
      testSet1 = training(partitionSplit1, curFold);
      trainSet2 = test(partitionSplit2, curFold);
      testSet2 = training(partitionSplit2, curFold);
    end
%     
    data1 = eventTimeData1(:, trainSet1, :);
    data2 = eventTimeData2(:, trainSet2, :);
    if(zeroComps)
      data1 = data1 - data1(1, :, :);
      data2 = data2 - data2(1, :, :);
    end
    if(wSize > 1)
      data1 = filter((1/wSize)*ones(1,wSize), 1, data1);
      data2 = filter((1/wSize)*ones(1,wSize), 1, data2);
    end
  
    
    mu1 = squeeze(nanmean(data1, 2)); %  Avg across trials
    mu2 = squeeze(nanmean(data2, 2)); %  Avg across trials
    s12 = zeros(size(mu1));
    for it1 = 1:size(data1, 1)
      n1 = sum(~isnan(sum(squeeze(data1(it1, :, :)),2)));
      n2 = sum(~isnan(sum(squeeze(data2(it1, :, :)),2)));
      %s12(it1, :) = sqrt((squeeze(nansum(data1(it1, :, :).^2,2))+squeeze(nansum(data2(it1, :, :).^2,2))-n1*squeeze(mu1(it1, :)').^2-n2*squeeze(mu2(it1, :)').^2)/(n1+n2));
      s12(it1, :) = sqrt(0.5*(nanvar(data1(it1, :, :), [], 2)+nanvar(data2(it1, :, :), [], 2)));
    end
    sVecOrigRaw = (mu1-mu2); % Classical definition
    sVecOrig = normalize(sVecOrigRaw, 'norm');
    sVecRaw = squeeze((mu1-mu2)./s12); % The one we will be using
    sVec = normalize(sVecRaw', 'norm')';
    sVecFold(:, :, curFold) = sVec;
    
    singleSVec = sVec(singleTimePointFrame, :);
    
    %Mdl = fitclinear([squeeze(data1(singleTimePointFrame, :, :)); squeeze(data2(singleTimePointFrame, :, :))], [zeros(size(data1,2), 1); ones(size(data2,2), 1)]);
    % HACK to fix sVec
    %sVec = repmat(sVec(29, :), [size(sVec, 1), 1]);
    
    % Hack to orthorgonalize from movement
%     selFrames = [35 35];
%     svMov = nanmean(A.sVecFull{itAnimal},3);
%     svMov = normalize(svMov', 'norm')';
%     svMov = svMov(selFrames(1), :);
%     for it2 = 1:size(sVec, 1)
%       svChoice = sVec(it2, :);
%       [q, r] = qr([svMov' svChoice'], 0);
%       sVec(it2, :) = q(:,2)';
%     end
    
    sVecMat(:, curFold) = pdist(sVec, 'corr');
    projDprimeTimeComp(:, :, curFold) = sVecRaw;
  %   [~, idx] = sort(abs(sVec), 'descend');
  %   sVec(idx(1:111)) = 0;
    % Now do the projections - Stim is special case, so data1 and data2 is
    % the same
    data1 = eventTimeData1(:, testSet1, :);
    data2 = eventTimeData2(:, testSet2, :);
    
    if(zeroComps)
      data1 = data1 - data1(1, :, :);
      data2 = data2 - data2(1, :, :);
    end
    if(wSize > 1)
      data1 = filter((1/wSize)*ones(1,wSize), 1, data1);
      data2 = filter((1/wSize)*ones(1,wSize), 1, data2);
    end
    projFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
    projFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
    singleProjFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
    singleProjFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
    for it0 = 1:size(eventTimeDataT, 1)
      for it1 = 1:size(data1, 2)
        projFull1(it0, it1) = squeeze(data1(it0, it1, :))'*sVec(it0, :)';
        singleProjFull1(it0, it1) = squeeze(data1(it0, it1, :))'*singleSVec(1,:)';
      end
      for it1 = 1:size(data2, 2)
        projFull2(it0, it1) = squeeze(data2(it0, it1, :))'*sVec(it0, :)';
        singleProjFull2(it0, it1) = squeeze(data2(it0, it1, :))'*singleSVec(1,:)';
      end
    end
    projMean1(:, curFold) = nanmean(projFull1, 2);
    projSEM1(:, curFold) = nanstd(projFull1, [], 2)./sqrt(sum(~isnan(projFull1),2));
    projMean2(:, curFold) = nanmean(projFull2, 2);
    projSEM2(:, curFold) = nanstd(projFull2, [], 2)./sqrt(sum(~isnan(projFull2),2));
    
    singleProjMean1(:, curFold) = nanmean(singleProjFull1, 2);
    singleProjSEM1(:, curFold) = nanstd(singleProjFull1, [], 2)./sqrt(sum(~isnan(singleProjFull1),2));
    singleProjMean2(:, curFold) = nanmean(singleProjFull2, 2);
    singleProjSEM2(:, curFold) = nanstd(singleProjFull2, [], 2)./sqrt(sum(~isnan(singleProjFull2),2));
    for it0 = 1:size(eventTimeDataT, 1)
      p1 = projFull1(it0, :);
      p2 = projFull2(it0, :);
      projDprime(it0, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');
      
      p1 = singleProjFull1(it0, :);
      p2 = singleProjFull2(it0, :);
      singleProjDprime(it0, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');
    end
    % Now the dprime per area
    for itA = 1:length(regionNames)
      validComp = find(regionIdx == itA);
      projFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
      projFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
      singleProjFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
      singleProjFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
      for it0 = 1:size(eventTimeDataT, 1)
        for it1 = 1:size(data1, 2)
          projFull1(it0, it1) = squeeze(data1(it0, it1, validComp))'*sVec(it0, validComp)';
          singleProjFull1(it0, it1) = squeeze(data1(it0, it1, validComp))'*singleSVec(1, validComp)';
        end
        for it1 = 1:size(data2, 2)
          projFull2(it0, it1) = squeeze(data2(it0, it1, validComp))'*sVec(it0, validComp)';
          singleProjFull2(it0, it1) = squeeze(data2(it0, it1, validComp))'*singleSVec(1, validComp)';
        end
      end

      for it0 = 1:size(eventTimeDataT, 1)
        p1 = projFull1(it0, :);
        p2 = projFull2(it0, :);
        projDprimeTimeCompPerArea(it0, itA, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');
        
        p1 = singleProjFull1(it0, :);
        p2 = singleProjFull2(it0, :);
        singleProjDprimeTimeCompPerArea(it0, itA, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');
      end
    end
    % Now the d prime per region
    %singleProjDprimeTimeCompPerAreaGroup = zeros(size(vectorModesFrames{1}, 1), 3, 3, kfoldN); % subareas / region definitions
    for itR = 1:length(areaGroupNames)
      for itA = 1:length(areaGroupNames{itR})
        validRegions = areaGroupIDs{itR}{itA};
        validComp = find(ismember(regionIdx, validRegions));
        projFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
        projFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
        singleProjFull1 = nan(size(eventTimeDataT, 1), size(data1, 2));
        singleProjFull2 = nan(size(eventTimeDataT, 1), size(data2, 2));
        for it0 = 1:size(eventTimeDataT, 1)
          for it1 = 1:size(data1, 2)
            projFull1(it0, it1) = squeeze(data1(it0, it1, validComp))'*sVec(it0, validComp)';
            singleProjFull1(it0, it1) = squeeze(data1(it0, it1, validComp))'*singleSVec(1, validComp)';
          end
          for it1 = 1:size(data2, 2)
            projFull2(it0, it1) = squeeze(data2(it0, it1, validComp))'*sVec(it0, validComp)';
            singleProjFull2(it0, it1) = squeeze(data2(it0, it1, validComp))'*singleSVec(1, validComp)';
          end
        end

        for it0 = 1:size(eventTimeDataT, 1)
          p1 = projFull1(it0, :);
          p2 = projFull2(it0, :);
          projDprimeTimeCompPerAreaGroup(it0, itR, itA, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');

          p1 = singleProjFull1(it0, :);
          p2 = singleProjFull2(it0, :);
          singleProjDprimeTimeCompPerAreaGroup(it0, itR, itA, curFold) = dPrimeCalc(p1(~isnan(p1))', p2(~isnan(p2))');
        end
      end
    end
    %projDprimeTimeCompPerArea
  end
  projMeanFull1{itAnimal} = projMean1;
  projSEMFull1{itAnimal} = projSEM1;
  projMeanFull2{itAnimal} = projMean2;
  projSEMFull2{itAnimal} = projSEM2;
  projDprimeFull{itAnimal} = projDprime;
  sVecMatFull{itAnimal} = sVecMat;
  projDprimeCompFull{itAnimal} = projDprimeComp;
  projDprimeTimeCompFull{itAnimal} = projDprimeTimeComp;
  wheelposFull1{itAnimal} = eventTimeDataWheel1;
  wheelposFull2{itAnimal} = eventTimeDataWheel2;
  pupilPosFull1{itAnimal} = eventTimeDataPupilPos1;
  pupilPosFull2{itAnimal} = eventTimeDataPupilPos2;
  sVecFull{itAnimal} = sVecFold;
  projDprimeTimeCompPerAreaFull{itAnimal} = projDprimeTimeCompPerArea;
  regionIdxFull{itAnimal} = regionIdx;
  
  singleProjMeanFull1{itAnimal} = singleProjMean1;
  singleProjSEMFull1{itAnimal} = singleProjSEM1;
  singleProjMeanFull2{itAnimal} = singleProjMean2;
  singleProjSEMFull2{itAnimal} = singleProjSEM2;
  singleProjDprimeFull{itAnimal} = singleProjDprime;
  singleProjDprimeTimeCompPerAreaFull{itAnimal} = singleProjDprimeTimeCompPerArea;
  singleProjDprimeTimeCompPerAreaGroupFull{itAnimal} = singleProjDprimeTimeCompPerAreaGroup;
  projDprimeTimeCompPerAreaGroupFull{itAnimal} = projDprimeTimeCompPerAreaGroup;
  
end
save('compMovements.mat', 'projMeanFull1', 'projMeanFull2', 'projSEMFull1', 'projSEMFull2', 'projDprimeFull', 'sVecMatFull', 'sVecFull', 'projDprimeCompFull', ...
     'wheelposFull1', 'wheelposFull2', 'projDprimeTimeCompPerAreaFull', 'projDprimeTimeCompFull', 'regionIdxFull', 'pupilPosFull1', 'pupilPosFull2', ...
     'singleProjMeanFull1', 'singleProjSEMFull1', 'singleProjMeanFull2', 'singleProjSEMFull2', 'singleProjDprimeFull', 'singleProjDprimeTimeCompPerAreaFull', 'singleProjDprimeTimeCompPerAreaGroupFull', 'projDprimeTimeCompPerAreaGroupFull',...
     'eventTimeDataT', 'areaGroupIDs', 'areaGroupNames', 'areaGroupNames2', 'singleTimePointFrame'); 
slog('Wheel mov done!', 'w');

%% Plot the projection split across animals

%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 18, 'height', 12);
p = panel();
p.margin = 10;
p.pack(4, 2);

for itAnimal = 1:Ndatasets
%for itAnimal = 1
  c = ceil(itAnimal/4);
  r = mod(itAnimal-1,4)+1;
  p(r ,c).select();
  cmap = lines(8);
  hold on;
  h1 = [];
  h2 = [];
  projMean = projMeanFull1{itAnimal};
  projSEM = projSEMFull1{itAnimal};
  time = eventTimeDataT;
  y = nanmean(projMean, 2);
  ysem = 2*nanmean(projSEM,2);
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  %h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
  h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{1}, trialSplits(itAnimal, 1)))];
  
  projMean = projMeanFull2{itAnimal};
  projSEM = projSEMFull2{itAnimal};
  y = nanmean(projMean, 2);
  ysem = 2*nanmean(projSEM,2);
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1+1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  h2 = [h2; plot(time(valid)', y(valid), '-', 'Color', cmap(1+1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{2}, trialSplits(itAnimal, 2)))];

  yl = ylim;
  plot([0 0], yl, 'k-');

 % title(datasetList.dataset{itAnimal});
  xlabel('time (s)');
  ylabel('Projection (a.u.)');
  xlim(time([1 end]));
  l = legend([h1 h2]);
  l.Box = 'off';
  l.Location = 'best';
  offsetAxes();
  spaceOutAxes();
end

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
%copygraphics(hFig, 'ContentType', 'vector');

%% Plot the dprime in time split
slmRegression = true;
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 18, 'height', 12);
p = panel();
p.margin = 10;
p.pack(4, 2);
slmFull = cell(Ndatasets+1);
slopeList = zeros(Ndatasets, 3);
kinkList = zeros(Ndatasets, 1);
for itAnimal = 1:Ndatasets
%for itAnimal = 1
  c = ceil(itAnimal/4);
  r = mod(itAnimal-1,4)+1;
  p(r ,c).select();
  cmap = lines(8);
  hold on;
  h1 = [];
  h2 = [];
  projMean = abs(projDprimeFull{itAnimal});

  time = eventTimeDataT+0.033;
  y = nanmean(projMean, 2);
  ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  %h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
  h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{1}, trialSplits(itAnimal, 1)))];
  plot(time(valid)', y(valid), 'o');
  if(slmRegression)
    t = time(valid)';
    y = y(valid);
    %[~, idx] = max(y);
    idx = find(y >= prctile(y, 80), 1, 'first');
    t = t(1:idx);
    y = y(1:idx);

    [slm, xp, yp] = slmengine(t, y,'degree',1,'plot','off','knots',3,'interiorknots','free', 'verbosity',0);
    %plot(t, slmeval(t, slm),'-ko');
    range2 = find(t >= slm.knots(1) & t < slm.knots(2));
    a = fit(t(range2)', y(range2), 'poly1');
    c = coeffvalues(a);
    b = confint(a);
    label = sprintf('%.2f (%.2f, %.2f)', c(1), b(1,1), b(2,1));
    h1 = [h1; plot(xp, yp, 'k-', 'DisplayName', label, 'LineWidth', 1)];
    slmFull{itAnimal} = slm;
    slopeList(itAnimal, :) = [c(1) b(1,1) b(2,1)];
    kinkList(itAnimal) = slm.knots(2);
  end
  yl = ylim;
  plot([1 1], yl, 'k-');
  xl = xlim;
  l = legend(h1);
  l.Box = 'off';
  l.Location = 'best';
 % title(datasetList.dataset{itAnimal});
  xlabel('time (s)');
  ylabel('d^\prime');
  xlim(eventTimeDataT([1 end]));
  yl = ylim;
  ylim([0 yl(2)]);
  offsetAxes();
  spaceOutAxes();
end


c = 2;
r = 4;
p(r ,c).select();
cmap = lines(8);
hold on;
h1 = [];
h2 = [];
pd = projDprimeFull([1 2 3 5 6 7]);
projMean = cell2mat(cellfun(@(x)mean(abs(x),2), pd(~cellfun(@isempty,pd)), 'UniformOutput', false)');

time = eventTimeDataT+0.033;
y = nanmean(projMean, 2);
ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(2, :), 'DisplayName', 'average')];
plot(time(valid)', y(valid), 'o');
if(slmRegression)
  t = time(valid)';
  y = y(valid);
  %[~, idx] = max(y);
  idx = find(y >= prctile(y, 80), 1, 'first');
  t = t(1:idx);
  y = y(1:idx);

  [slm, xp, yp] = slmengine(t, y,'degree',1,'plot','off','knots',3,'interiorknots','free', 'verbosity',0);
  %plot(t, slmeval(t, slm),'-ko');
  range2 = find(t >= slm.knots(1) & t < slm.knots(2));
  a = fit(t(range2)', y(range2), 'poly1');
  c = coeffvalues(a);
  b = confint(a);
  label = sprintf('%.2f (%.2f, %.2f)', c(1), b(1,1), b(2,1));
  h1 = [h1; plot(xp, yp, 'k-', 'DisplayName', label, 'LineWidth', 1)];
  slmFull{end} = slm;
end
yl = ylim;
plot([1 1], yl, 'k-');
xl = xlim;
l = legend(h1);
l.Box = 'off';
l.Location = 'best';
% title(datasetList.dataset{itAnimal});
xlabel('time (s)');
ylabel('d^\prime');
xlim(eventTimeDataT([1 end]));
yl = ylim;
ylim([0 yl(2)]);
offsetAxes();
spaceOutAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

%exportgraphics(hFig, 'movSlopes.pdf');

%% Plot the single projection split across animals

%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 18, 'height', 12);
p = panel();
p.margin = 10;
p.pack(4, 2);

for itAnimal = 1:Ndatasets
%for itAnimal = 1
  c = ceil(itAnimal/4);
  r = mod(itAnimal-1,4)+1;
  p(r ,c).select();
  cmap = lines(8);
  hold on;
  h1 = [];
  h2 = [];
  projMean = singleProjMeanFull1{itAnimal};
  projSEM = singleProjSEMFull1{itAnimal};
  time = eventTimeDataT+0.033;
  y = nanmean(projMean, 2);
  ysem = 2*nanmean(projSEM,2);
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  %h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
  h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{1}, trialSplits(itAnimal, 1)))];
  
  projMean = singleProjMeanFull2{itAnimal};
  projSEM = singleProjSEMFull2{itAnimal};
  y = nanmean(projMean, 2);
  ysem = 2*nanmean(projSEM,2);
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1+1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  h2 = [h2; plot(time(valid)', y(valid), '-', 'Color', cmap(1+1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{2}, trialSplits(itAnimal, 2)))];

  yl = ylim;
  plot([0 0], yl, 'k-');

 % title(datasetList.dataset{itAnimal});
  xlabel('time (s)');
  ylabel('Projection (a.u.)');
  xlim(eventTimeDataT([1 end]));
  l = legend([h1 h2]);
  l.Box = 'off';
  l.Location = 'best';
  offsetAxes();
  spaceOutAxes();
end

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
copygraphics(hFig, 'ContentType', 'vector');

%% Plot the single dprime in time split
slmRegression = false;
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 18, 'height', 12);
p = panel();
p.margin = 10;
p.pack(4, 2);
slmFull = cell(Ndatasets+1);
slopeList = zeros(Ndatasets, 3);
kinkList = zeros(Ndatasets, 1);
time = eventTimeDataT;
%time = eventTimeDataT+0.033;
for itAnimal = 1:Ndatasets
%for itAnimal = 1
  c = ceil(itAnimal/4);
  r = mod(itAnimal-1,4)+1;
  p(r ,c).select();
  cmap = lines(8);
  hold on;
  h1 = [];
  h2 = [];
  projMean = abs(singleProjDprimeFull{itAnimal});

  
  y = nanmean(projMean, 2);
  ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));
  upper = y + ysem;
  lower = y - ysem;
  valid = find(~isnan(y));
  fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
  %h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
  h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s (%d)', vectorModes{1}, trialSplits(itAnimal, 1)))];
  %plot(time(valid)', y(valid), 'o');
  if(slmRegression)
    t = time(valid)';
    y = y(valid);
    %[~, idx] = max(y);
    idx = find(y >= prctile(y, 80), 1, 'first');
    t = t(1:idx);
    y = y(1:idx);

    [slm, xp, yp] = slmengine(t, y,'degree',1,'plot','off','knots',3,'interiorknots','free', 'verbosity',0);
    %plot(t, slmeval(t, slm),'-ko');
    range2 = find(t >= slm.knots(1) & t < slm.knots(2));
    a = fit(t(range2)', y(range2), 'poly1');
    c = coeffvalues(a);
    b = confint(a);
    label = sprintf('%.2f (%.2f, %.2f)', c(1), b(1,1), b(2,1));
    h1 = [h1; plot(xp, yp, 'k-', 'DisplayName', label, 'LineWidth', 1)];
    slmFull{itAnimal} = slm;
    slopeList(itAnimal, :) = [c(1) b(1,1) b(2,1)];
    kinkList(itAnimal) = slm.knots(2);
  end
  yl = ylim;
  plot([1 1], yl, 'k-');
  xl = xlim;
  l = legend(h1);
  l.Box = 'off';
  l.Location = 'best';
 % title(datasetList.dataset{itAnimal});
  xlabel('time (s)');
  ylabel('d^\prime');
  xlim(eventTimeDataT([1 end]));
  yl = ylim;
  ylim([0 yl(2)]);
  offsetAxes();
  spaceOutAxes();
end


c = 2;
r = 4;
p(r ,c).select();
cmap = lines(8);
hold on;
h1 = [];
h2 = [];
pd = singleProjDprimeFull([1 2 3 5 6 7]);
projMean = cell2mat(cellfun(@(x)mean(abs(x),2), pd(~cellfun(@isempty,pd)), 'UniformOutput', false)');

y = nanmean(projMean, 2);
ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(2, :), 'DisplayName', 'average')];
%plot(time(valid)', y(valid), '-');
if(slmRegression)
  t = time(valid)';
  y = y(valid);
  %[~, idx] = max(y);
  idx = find(y >= prctile(y, 80), 1, 'first');
  t = t(1:idx);
  y = y(1:idx);

  [slm, xp, yp] = slmengine(t, y,'degree',1,'plot','off','knots',3,'interiorknots','free', 'verbosity',0);
  %plot(t, slmeval(t, slm),'-ko');
  range2 = find(t >= slm.knots(1) & t < slm.knots(2));
  a = fit(t(range2)', y(range2), 'poly1');
  c = coeffvalues(a);
  b = confint(a);
  label = sprintf('%.2f (%.2f, %.2f)', c(1), b(1,1), b(2,1));
  h1 = [h1; plot(xp, yp, 'k-', 'DisplayName', label, 'LineWidth', 1)];
  slmFull{end} = slm;
end
yl = ylim;
plot([1 1], yl, 'k-');
xl = xlim;
l = legend(h1);
l.Box = 'off';
l.Location = 'best';
% title(datasetList.dataset{itAnimal});
xlabel('time (s)');
ylabel('d^\prime');
xlim(eventTimeDataT([1 end]));
yl = ylim;
ylim([0 yl(2)]);
offsetAxes();
spaceOutAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

exportgraphics(hFig, 'movSlopes.pdf');

%% Prepare drawings to plot our 10 allen-based areas

R = imref2d([67, 67]);
dsetAllen = wf.dataset.load(projectFolder, 'allen');
areasAllen = dsetAllen.data.areas;
retData = dsetAllen.data.areas;  
%areasAllen = round(imwarp(retData, dsetAllen.data.alignment.transform, 'nearest', 'OutputView', R));
% Save also our redefinition of Allen areas
areasAllenLocaNMF = zeros(size(areasAllen));
areasNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RSP'};
areasID = {14, 29, 25, 26, 4, 27, 21, 11, [7 12 19 32 34], [28 31]};
for it1 = 1:length(areasNames)
  for it2 = 1:length(areasID{it1})
    areasAllenLocaNMF(areasAllen == areasID{it1}(it2)) = it1;
  end
end

%%% Bit to be able to plot the allen data on top of our figures
areasAllen = dsetAllen.data.areas;
areasNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RSP'};
areasID = 1:length(areasNames);

contourList = cell(size(areasID));
centerList = cell(size(areasID));
labelCenterList = cell(size(areasID));
contourListTransform = cell(size(areasID));
labelCenterListTransform = cell(size(areasID));
for it0 = 1:length(areasID)
  valid = find(areasAllenLocaNMF == areasID(it0));
  areaBoundaries = zeros(size(areasAllenLocaNMF));
  areaBoundaries(valid) = 1;
  prop = regionprops(areaBoundaries);

  center = round(prop.Centroid);
  firstPixel = find(areaBoundaries(center(2), :) == 1, 1, 'first');

  contour = bwtraceboundary(areaBoundaries, [center(2) firstPixel], 'W', 8);
  contourList{it0} = [contour(:,2) contour(:,1)];
  %xn = contour(:,2);
  %yn = contour(:,1);
  [xn,yn] = transformPointsForward(dsetAllen.data.alignment.transform, contour(:,2), contour(:,1));
  contourListTransform{it0} = [xn yn];
  centerList{it0} = center;
  %if(it0 <= length(areasNames))
    labelCenterList{it0} = center;
    [xn,yn] = transformPointsForward(dsetAllen.data.alignment.transform, center(1), center(2));
    labelCenterListTransform{it0} = [xn yn];
    %labelCenterListTransform{it0} = center;
  %end  
end


%% The figure plot
load('compMovements.mat');
%
vectorModes = {'mov', 'no mov'};
wSize = 3;
eventTimeDataT = dataSVD.time(1:60)-1+dataSVD.time(2);
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.1 0.35 0.35 0.2]);

selAnimal = 1;

itAnimal = selAnimal;

%%% The projection
p(1 ,2).select();
cmap = cblindmap();
cmap = cmap(1:end,:);
hold on;
h1 = [];
h2 = [];
projMean = singleProjMeanFull1{itAnimal};
projSEM = singleProjSEMFull1{itAnimal};
time = eventTimeDataT+0.033;
y = nanmean(projMean, 2);
ysem = 2*nanmean(projSEM,2);
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}))];

projMean = singleProjMeanFull2{itAnimal};
projSEM = singleProjSEMFull2{itAnimal};
y = nanmean(projMean, 2);
ysem = 2*nanmean(projSEM,2);
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1+1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(time(valid)', y(valid), '-', 'Color', cmap(1+1, :), 'DisplayName', sprintf('%s', vectorModes{2}))];
yl = ylim;
ylim([0 yl(2)*1.2]);
plot([0 0], ylim, 'k--');


x = [time(singleTimePointFrame-wSize+1) time(singleTimePointFrame)];
y = ylim;
hp1 = fill([x(1) x(1) x(end) x(end)], [y(2)*0.95 y(2) y(2) y(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');

xlabel('time (s)');
ylabel('Projection (a.u.)');
%xlim(eventTimeDataT([1 end]));
xlim([-1 1]);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();


  
%%% The d prime

p(1, 3).select();

hold on;
h1 = [];
h2 = [];
%pd = singleProjDprimeFull([1 2 3 5 6 7]);
pd = singleProjDprimeFull;
projMean = cell2mat(cellfun(@(x)mean(abs(x),2), pd(~cellfun(@isempty,pd)), 'UniformOutput', false)');

y = nanmean(projMean, 2);
ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));

upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(3, :), 'DisplayName', 'average')];
%plot(time(valid)', y(valid), '-');
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
[a,b]=max(y);
[y(b) ysem(b)]
x = [time(singleTimePointFrame-wSize+1) time(singleTimePointFrame)];
y = ylim;
hp1 = fill([x(1) x(1) x(end) x(end)], [y(2)*0.95 y(2) y(2) y(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');

l = legend(h1);
l.Box = 'off';
l.Location = 'best';
% title(datasetList.dataset{itAnimal});
xlabel('time (s)');
ylabel('d''');
%xlim(eventTimeDataT([1 end]));
xlim([-1 1]);
offsetAxes();
spaceOutAxes();

% Now the d prime spatial map
p(1, 4).select();

componentsPixels = cell(length(regionNames), 1);
dprimeTrajs = nan(length(eventTimeDataT), length(regionNames));
for it2 = 1:length(regionNames)
  componentsPixels{it2} = find(dset.data.areasAllenLocaNMF == it2);
  %pd = singleProjDprimeTimeCompPerAreaFull([1 2 3 5 6 7]);
  pd = singleProjDprimeTimeCompPerAreaFull;
  dprimeData = cell2mat(cellfun(@(x)x(:, it2), cellfun(@(x)mean(x,3), pd, 'UniformOutput', false), 'UniformOutput', false)');
  dprimeTrajs(:, it2) = nanmean(dprimeData, 2);
end
[~, bestFrame] = max(nanmean(dprimeTrajs,2));
dprimeSet = dprimeTrajs(bestFrame, :);

movData = nan(size(dset.data.areasAllenLocaNMF));
for it2 = 1:length(regionNames)
  movData(componentsPixels{it2}) = dprimeSet(it2);
end

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
%c.Location = 'west';
caxis([0 max(movData(:))]);
caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
colormap(parula(10));
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
%exportgraphics(hFig, 'fig2Mov.pdf', 'ContentType', 'vector');
%copygraphics(hFig, 'ContentType', 'vector');


%% The d prime split univar

%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.2 0.35 0.35 0.1]);


% Now the d prime spatial map
p(1, 1).select();

componentsPixels = cell(length(regionNames), 1);
dprimeTrajs = nan(length(eventTimeDataT), length(regionNames));
for it2 = 1:length(regionNames)
  componentsPixels{it2} = find(dset.data.areasAllenLocaNMF == it2);
  %pd = singleProjDprimeTimeCompPerAreaFull([1 2 3 4 5 6 7]);
  pd = singleProjDprimeTimeCompPerAreaFull;
  dprimeData = cell2mat(cellfun(@(x)x(:, it2), cellfun(@(x)mean(x,3), pd, 'UniformOutput', false), 'UniformOutput', false)');
  dprimeTrajs(:, it2) = nanmean(dprimeData, 2);
end
[~, bestFrame] = max(nanmean(dprimeTrajs,2));
dprimeSet = dprimeTrajs(bestFrame, :);

movData = nan(size(dset.data.areasAllenLocaNMF));
for it2 = 1:length(regionNames)
  movData(componentsPixels{it2}) = dprimeSet(it2);
end

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
%c.Location = 'west';
%caxis([0 max(movData(:))]);
caxis([0 1.2]);
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
colormap(parula(8));
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

% Now the d prime spatial map
p(1, 2).select();

data = cell2mat(cellfun(@(x)x(bestFrame, :), cellfun(@(x)nanmean(x,3), singleProjDprimeTimeCompPerAreaFull, 'UniformOutput', false),'UniformOutput', false));
dataG = cell2mat(cellfun(@(x)x(bestFrame, :), cellfun(@(x)nanmean(x,2), singleProjDprimeFull, 'UniformOutput', false),'UniformOutput', false));
data = [dataG, data];
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Compression', 10, 'PointSize', 4);
set(gca, 'XTick', 1:length(regionNames)+1);
set(gca, 'XTickLabel', {'global', regionNames{:}});
set(gca, 'XTickLabelRotation', 90);
% for it1 = 1:length(hList)
%   hList(it1).SizeData = 8;
% end


p(1, 3).select();
time = eventTimeDataT(2:end);
cmap = divergingBlueRedCmap(256);
sVec = mean(cell2mat(cellfun(@(x)mean(x,2),sVecMatFull,'UniformOutput',false)'),2);
mm = squareform(sVec);
imagesc(1-mm(2:end,2:end), 'XData', time, 'YData', time);
xlim(time([1 end]));
ylim(time([1 end]));
hold on;
plot([0 0], ylim, 'k--');
plot(xlim, [0 0], 'k--');
xlabel('time (s)');
ylabel('time (s)');
axis square ij;
caxis([-1 1]);
colormap(gca,cmap);
colorbar;


set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
exportgraphics(hFig, 'fig2MovED.pdf', 'ContentType', 'vector');

%%
y = nanmean(data);
ys = nanstd(data)./sqrt(sum(~isnan(data)));
[~, nidx] = sort(y, 'descend');
figure;
errorbar(y(nidx), ys(nidx));
names = {'global', regionNames{:}};
set(gca, 'XTick', 1:11);
set(gca, 'XTickLabels', names(nidx));
set(gca, 'XTickLabelRotation', 90);

%%
figure;
errorbar(y(nidx), ys(nidx));
names = {'global', regionNames{:}};
set(gca, 'XTick', 1:11);
set(gca, 'XTickLabels', names(nidx));
set(gca, 'XTickLabelRotation', 90);
selSource = 1;
pList = zeros(length(nidx), 1);
for it1 = 2:length(nidx)
  %[p, h] = signrank(data(:, selSource), data(:, nidx(it1)));
  [h, p] = ttest(data(:, selSource), data(:, nidx(it1)));
  pList(it1) = p;
  if(p < 0.05)
    sigstar([selSource it1], p);
  end
end
%%
figure;
hold on;
nidxList = zeros(size(data));
for it1 = 1:size(data, 1)
  [~, nidx] = sort(data(it1, :), 'descend');
  plot(nidx, 1:11, 'o-');
  nidxList(it1, nidx) = 1:11;
end

%%
figure;
errorbar(mean(nidxList,1),std(nidxList)/sqrt(7),'o')

%%

%% The figure plot
load('compMovements.mat');
%
vectorModes = {'mov', 'no mov'};
wSize = 3;
eventTimeDataT = dataSVD.time(1:60)-1+dataSVD.time(2);
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.1 0.35 0.35 0.2]);

selAnimal = 1;

itAnimal = selAnimal;

%%% The projection
p(1 ,2).select();
cmap = cblindmap();
cmap = cmap(1:end,:);
hold on;
h1 = [];
h2 = [];
projMean = projMeanFull1{itAnimal};
projSEM = projSEMFull1{itAnimal};
time = eventTimeDataT+0.033;
y = nanmean(projMean, 2);
ysem = 2*nanmean(projSEM,2);
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}))];

projMean = projMeanFull2{itAnimal};
projSEM = projSEMFull2{itAnimal};
y = nanmean(projMean, 2);
ysem = 2*nanmean(projSEM,2);
upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(1+1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(time(valid)', y(valid), '-', 'Color', cmap(1+1, :), 'DisplayName', sprintf('%s', vectorModes{2}))];
yl = ylim;
ylim([0 yl(2)*1.2]);
plot([0 0], ylim, 'k--');


x = [time(singleTimePointFrame-wSize+1) time(singleTimePointFrame)];
y = ylim;
hp1 = fill([x(1) x(1) x(end) x(end)], [y(2)*0.95 y(2) y(2) y(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');

xlabel('time (s)');
ylabel('Projection (a.u.)');
%xlim(eventTimeDataT([1 end]));
xlim([-1 1]);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();


  
%%% The d prime

p(1, 3).select();

hold on;
h1 = [];
h2 = [];
%pd = singleProjDprimeFull([1 2 3 5 6 7]);
pd = projDprimeFull;
projMean = cell2mat(cellfun(@(x)mean(abs(x),2), pd(~cellfun(@isempty,pd)), 'UniformOutput', false)');

y = nanmean(projMean, 2);
ysem = nanstd(projMean, [], 2)/sqrt(size(projMean,2));

upper = y + ysem;
lower = y - ysem;
valid = find(~isnan(y));
fill([time(valid)' fliplr(time(valid)')],[upper(valid)' fliplr(lower(valid)')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(1, :), 'DisplayName', datasetList.dataset{itAnimal})];
h1 = [h1; plot(time(valid)', y(valid), '-', 'Color', cmap(3, :), 'DisplayName', 'average')];
%plot(time(valid)', y(valid), '-');
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
[a,b]=max(y);
[y(b) ysem(b)]
x = [time(singleTimePointFrame-wSize+1) time(singleTimePointFrame)];
y = ylim;
hp1 = fill([x(1) x(1) x(end) x(end)], [y(2)*0.95 y(2) y(2) y(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');

l = legend(h1);
l.Box = 'off';
l.Location = 'best';
% title(datasetList.dataset{itAnimal});
xlabel('time (s)');
ylabel('d''');
%xlim(eventTimeDataT([1 end]));
xlim([-1 1]);
offsetAxes();
spaceOutAxes();

% Now the d prime spatial map
p(1, 4).select();

componentsPixels = cell(length(regionNames), 1);
dprimeTrajs = nan(length(eventTimeDataT), length(regionNames));
for it2 = 1:length(regionNames)
  componentsPixels{it2} = find(dset.data.areasAllenLocaNMF == it2);
  %pd = singleProjDprimeTimeCompPerAreaFull([1 2 3 5 6 7]);
  pd = singleProjDprimeTimeCompPerAreaFull;
  dprimeData = cell2mat(cellfun(@(x)x(:, it2), cellfun(@(x)mean(x,3), pd, 'UniformOutput', false), 'UniformOutput', false)');
  dprimeTrajs(:, it2) = nanmean(dprimeData, 2);
end
[~, bestFrame] = max(nanmean(dprimeTrajs,2));
dprimeSet = dprimeTrajs(bestFrame, :);

movData = nan(size(dset.data.areasAllenLocaNMF));
for it2 = 1:length(regionNames)
  movData(componentsPixels{it2}) = dprimeSet(it2);
end

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
%c.Location = 'west';
caxis([0 max(movData(:))]);
caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
colormap(parula(10));
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
exportgraphics(hFig, 'fig2Mov.pdf', 'ContentType', 'vector');
%copygraphics(hFig, 'ContentType', 'vector');


%% The d prime split univar

%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.2 0.35 0.35 0.1]);


% Now the d prime spatial map
p(1, 1).select();

componentsPixels = cell(length(regionNames), 1);
dprimeTrajs = nan(length(eventTimeDataT), length(regionNames));
for it2 = 1:length(regionNames)
  componentsPixels{it2} = find(dset.data.areasAllenLocaNMF == it2);
  %pd = singleProjDprimeTimeCompPerAreaFull([1 2 3 4 5 6 7]);
  pd = singleProjDprimeTimeCompPerAreaFull;
  dprimeData = cell2mat(cellfun(@(x)x(:, it2), cellfun(@(x)mean(x,3), pd, 'UniformOutput', false), 'UniformOutput', false)');
  dprimeTrajs(:, it2) = nanmean(dprimeData, 2);
end
[~, bestFrame] = max(nanmean(dprimeTrajs,2));
dprimeSet = dprimeTrajs(bestFrame, :);

movData = nan(size(dset.data.areasAllenLocaNMF));
for it2 = 1:length(regionNames)
  movData(componentsPixels{it2}) = dprimeSet(it2);
end

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
%c.Location = 'west';
%caxis([0 max(movData(:))]);
caxis([0 1.2]);
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
colormap(parula(8));
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

% Now the d prime spatial map
p(1, 2).select();

data = cell2mat(cellfun(@(x)x(bestFrame, :), cellfun(@(x)nanmean(x,3), singleProjDprimeTimeCompPerAreaFull, 'UniformOutput', false),'UniformOutput', false));
dataG = cell2mat(cellfun(@(x)x(bestFrame, :), cellfun(@(x)nanmean(x,2), singleProjDprimeFull, 'UniformOutput', false),'UniformOutput', false));
data = [dataG, data];
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Compression', 10, 'PointSize', 4);
set(gca, 'XTick', 1:length(regionNames)+1);
set(gca, 'XTickLabel', {'global', regionNames{:}});
set(gca, 'XTickLabelRotation', 90);
% for it1 = 1:length(hList)
%   hList(it1).SizeData = 8;
% end


p(1, 3).select();
time = eventTimeDataT(2:end);
cmap = divergingBlueRedCmap(256);
sVec = mean(cell2mat(cellfun(@(x)mean(x,2),sVecMatFull,'UniformOutput',false)'),2);
mm = squareform(sVec);
imagesc(1-mm(2:end,2:end), 'XData', time, 'YData', time);
xlim(time([1 end]));
ylim(time([1 end]));
hold on;
plot([0 0], ylim, 'k--');
plot(xlim, [0 0], 'k--');
xlabel('time (s)');
ylabel('time (s)');
axis square ij;
caxis([-1 1]);
colormap(gca,cmap);
colorbar;


set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%export_fig('2ext_AreaMasks.pdf','-nocrop');
exportgraphics(hFig, 'fig2MovED.pdf', 'ContentType', 'vector');