%%% -----------------------------------------------------------------------
%%% This script will compute the state vectors used throughout the manuscript
%%% See main.m
%%% -----------------------------------------------------------------------
%%% We provide the computation of state vectors for saccades as a coding sample.
%%% Rest of the scripts will be available at publication time.

%% Load project data

clearvars;
addpath(genpath(pwd));
projectFolder = ['..' filesep 'widefieldChoice' filesep 'data' filesep 'excMice'];
load([projectFolder filesep 'projectData.mat']);


%% Preload data

dsetAnimal = cell(Ndatasets, 1);
dataLocaAnimal = cell(Ndatasets, 1);
dataSVDAnimal = cell(Ndatasets, 1);
dataBehaviorAnimal = cell(Ndatasets, 1);
baseCompAnimal = cell(Ndatasets, 1);
for itAnimal = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{itAnimal});
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


%% Compute state vectors
%%% -----------------------------------------------------------------------
%%% Saccades
%%% -----------------------------------------------------------------------

vectorModes = {'sacc','no sacc'}; % Vector will be V1 - V2
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
    idx = find(~isnan(dataBehavior.saccades(:, it1)));
    firstIdx = idx(find(dataSVD.time(idx) > dataSVD.time(movWindow(1)) & dataSVD.time(idx) <= dataSVD.time(movWindow(end)), 1, 'first')); % Looking for the 1st mov on the window
    lastIdxNoMov = idx(find(dataSVD.time(idx) > dataSVD.time(noMovWindow(1)) & dataSVD.time(idx) <= dataSVD.time(noMovWindow(end)), 1, 'last')); % Looking for the last mov on the no window
    if(dataSVD.time(firstIdx) - dataSVD.time(lastIdxNoMov) < 0.5) % At least half a sec separation
      choices2(it1) = 2;
      continue;
    end
    if(sum(~isnan(dataBehavior.wheelMovementsSz(firstIdx-15:firstIdx, it1))) > 0) % Half a sec without mov before saccades
      choices2(it1) = 3;
      continue;
    end
    idxSac = find(~isnan(dataBehavior.wheelMovementsSz(:, it1)));
    interferingSacc = idxSac(find(dataSVD.time(idxSac) > 1 & dataSVD.time(idxSac) <= 1.5, 1, 'first'));  % Half a sec without mov after sacc
    if(~isempty(interferingSacc))
      choices2(it1) = 3;
      continue;
    end
    if(dataBehavior.saccades(firstIdx, it1) ~= 0)
      choices2(it1) = 1;
    elseif(dataBehavior.saccades(firstIdx, it1) == 0)
      choices2(it1) = 0;
    end
  end
  outcome = dset.data.trialsSummary.outcome(dataBehavior.trialID);

  selTrial1 = find(abs(choices2) == 1);
  selTrial2 = find(abs(choices2) ~= 1);
  [length(selTrial1) length(selTrial2)]
  
  movTime1 = zeros(size(selTrial1));
  movTime2 = zeros(size(selTrial2));
  % Let's do 1 sec before and after
  eventTimeData1 = nan(60, length(movTime1), size(locaReshaped, 3));
  eventTimeDataWheel1 = nan(60, length(movTime1));
  eventTimeDataPupilPos1 = nan(60, 2, length(movTime1));
  
  eventTimeDataT = dataSVD.time(1:60)-1+dataSVD.time(2);
  for it0 = 1:length(selTrial1)
    movTime1(it0) = movWindow(1)-1+find(~isnan(dataBehavior.saccades(movWindow, selTrial1(it0))), 1, 'first');
    
    firstFrame = max(movTime1(it0)-29, 1);
    lastFrame = min(movTime1(it0)+30, size(dataBehavior.saccades, 1));
    eventTimeData1((firstFrame:lastFrame)-movTime1(it0)+30, it0, :) = locaReshaped(firstFrame:lastFrame, selTrial1(it0), :);
    eventTimeDataWheel1((firstFrame:lastFrame)-movTime1(it0)+30, it0) = dataBehavior.wheelP(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.wheelP(firstFrame, selTrial1(it0));
    eventTimeDataPupilPos1((firstFrame:lastFrame)-movTime1(it0)+30, 1, it0) = dataBehavior.pupilX(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.pupilX(firstFrame, selTrial1(it0));
    eventTimeDataPupilPos1((firstFrame:lastFrame)-movTime1(it0)+30, 2, it0) = dataBehavior.pupilY(firstFrame:lastFrame, selTrial1(it0))-dataBehavior.pupilY(firstFrame, selTrial1(it0));
  end
  eventTimeData2 = nan(60, length(movTime2), size(locaReshaped, 3));
  eventTimeDataWheel2 = nan(60, length(movTime2));
  eventTimeDataPupilPos2 = nan(60, 2, length(movTime2));
  
  for it0 = 1:length(selTrial2)
    movTime2(it0) = movWindow(1)-1+randi(length(movWindow)); % Center these trials randomly on the same interval
    firstFrame = max(movTime2(it0)-29, 1);
    lastFrame = min(movTime2(it0)+30, size(dataBehavior.saccades, 1));
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
      s12(it1, :) = sqrt(0.5*(nanvar(data1(it1, :, :), [], 2)+nanvar(data2(it1, :, :), [], 2)));
    end
    sVecOrigRaw = (mu1-mu2); % Classical definition
    sVecOrig = normalize(sVecOrigRaw, 'norm');
    sVecRaw = squeeze((mu1-mu2)./s12); % The one we will be using
    sVec = normalize(sVecRaw', 'norm')';
    sVecFold(:, :, curFold) = sVec;
    
    singleSVec = sVec(singleTimePointFrame, :);
    
    sVecMat(:, curFold) = pdist(sVec, 'corr');
    projDprimeTimeComp(:, :, curFold) = sVecRaw;

    % Now do the projections
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
save('compSacc.mat', 'projMeanFull1', 'projMeanFull2', 'projSEMFull1', 'projSEMFull2', 'projDprimeFull', 'sVecMatFull', 'sVecFull', 'projDprimeCompFull', ...
     'wheelposFull1', 'wheelposFull2', 'projDprimeTimeCompPerAreaFull', 'projDprimeTimeCompFull', 'regionIdxFull', 'pupilPosFull1', 'pupilPosFull2', ...
     'singleProjMeanFull1', 'singleProjSEMFull1', 'singleProjMeanFull2', 'singleProjSEMFull2', 'singleProjDprimeFull', 'singleProjDprimeTimeCompPerAreaFull', 'singleProjDprimeTimeCompPerAreaGroupFull', 'projDprimeTimeCompPerAreaGroupFull',...
     'eventTimeDataT', 'areaGroupIDs', 'areaGroupNames', 'areaGroupNames2', 'singleTimePointFrame'); 
slog('Saccades done!', 'w');

%%
