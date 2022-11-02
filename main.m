%%% -----------------------------------------------------------------------
%%% Main script to create and process a new (or existing) project
%%% A project consists of a list of datasets (animals)
%%% Used data should already been preprocessed as in: 
%%% https://doi.org/10.1016/j.celrep.2021.109377
%%% Attention separates sensory and motor signals in the mouse visual cortex
%%% M. Abdolrahmani, D. R. Lyamzin, R. Aoki, A. Benucci, Cell Reports 2021
%%%
%%% From the prev reference, each animal data is stored in a single structure,
%%% L, that contains information of every trial.
%%% L relevant contents (N trials in S sessions):
%%% .sess [Nx1] - session ID of each trial (1 to T)
%%% .sessname [Tx1] - name of each session
%%% .trial [Nx1] - trial idx within a session
%%% .skip [Nx1] - wether that trial should be skipped (bool)
%%% .outcome [Nx1] - -1/0/1 for incorrect / time out /correct
%%% .choicedir [Nx1] - -1/1 for left and right choices
%%% .aux - metadata
%%% .ori [Nx2] - orientation in degrees of the 2 presentedf stimulii
%%% .fps [Nx2] - framerate of imaging and eye-tracking data
%%%
%%% Each trial data (imaging and behavior) is stored in a separate structure, found
%%% on separate files (usually dset.info.trialsSummaryFolder/tile/trialNumber.mat )
%%% Relevant contents:
%%% .fps [1x2] - framerate of imaging and eye-tracking data
%%% .sessnum - session ID the trial belongs to
%%% .trialnum - trial ID within the session
%%% .time - frames timestamps
%%% .resppix - preprocessing metadata
%%% .resp - [1xPxT] - imaging data (dF/F) for P pixels in T frames
%%% .wheelpos [Kx2] - wheelposition for K samples (time (aligned to imaging times), encoding value)
%%% .eventtimes - struct containing trial temporal data (time of stimulus presentation, closed loop start, closed loop end, etc.)
%%% .ppos [2xR] - pupil position for R frames (framerate might be different from imaging one), (X and Y)
%%% .parea [1xR] - measured pupil area on each frame
%%%
%%% -----------------------------------------------------------------------

%% Create a new project
clearvars;
addpath(genpath(pwd));

%%% External locations

% L structure location
LfileRoot = 'F:\dataWF\'; % Full path should be: [LfileRoot datasetList.dataset{it} '\trials\L_tile.mat']
% Retinotopies location - dset.info.retinotopy is read from the datasetLit (table below)
retinotopyRootFolder = '\\172.17.150.231\DATA\MOUSE\IMAGING\GCAMP\'; % Full path should be: [retinotopyRootFolder dset.info.retinotopy{1} '\ANALYZED\' dset.info.retinotopy{2} '\AnalyzedRet.mat']

projectFolder = [pwd filesep 'data' filesep 'excMice'];
% Generate a table with some metadata from each animal: Animal name,
% invalid sessions, retinotopy file, Numeric IDs of retinotopic areas used
% for alignment (V1, PM, RL, AM) and trained reference orientation (either
% 0 or 90 degrees)
% Invalid sessions come from pasing the TO < 0.2 and PerfNoTO > 0.6 and consistent imaging frames criteria (no gaps)
datasetList = {'A15098', [4 5 6 13 14], {'M151125_MA' '1'}, 3, 2, 7, 4, 0;
               'A15309', [5 8 10 13 15:18 20:22], {'M160415_RA' '1'}, 5, 3, 7, 4, 0;
               'A15100', [1 2 5 7 9 10], {'M160215_MA' '1'}, 3, 4, 6, 5, 0;
               'A15312', [6 7 10 13:15 26 28], {'M160415_RA' '2'}, 5, 3, 6, 4, 0;
               'A15301', [1:5 7 13:15 25 27], {'M160413_RA' '2'}, 4, 1, 5, 3, 90;
               'A15352', [1 2 4 6 7 10 13:16 18:21 24 25 27 32], {'M160511_MA' '1'}, 6, 4, 7, 3, 90;
               'A16032', [2:7 20 24:28 36:41 51 54:57 60:64],   {'M160923_RA' '1'}, 6, 3, 7, 1, 90};
datasetList = table(datasetList(:,1), datasetList(:,2), datasetList(:,3), datasetList(:,4), datasetList(:,5), datasetList(:,6), datasetList(:,7), datasetList(:,8), ...
                    'variableNames', {'dataset', 'invalidSessions', 'retinotopy', 'V1', 'PM', 'RL', 'AM', 'reference'}); 
datasetList = sortrows(datasetList,1);
Ndatasets = size(datasetList,1);

save([projectFolder filesep 'projectData.mat'], 'LfileRoot', 'retinotopyRootFolder', 'projectFolder', 'datasetList', 'Ndatasets');

%% Create the datasets - only need to do this once - Change to invalid sessions

for it = 1:Ndatasets
  dset = wf.dataset(datasetList.dataset{it}, 'rootFolder', projectFolder, 'comments', 'exc mice main analysis', 'tags', {'exc'});
  % Assign the session list
  dset.info.invalidSessions = datasetList.invalidSessions{it};
  dset = dset.resetValidSessions();
  dset = dset.setValidSessions(setxor(dset.getValidSessions(), dset.info.invalidSessions));
  dset.save();
end

%% Assign the trials data

% Needs the summary file - File containing the L structure
for it = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it});
  dset.rootFolder = [pwd filesep 'data' filesep 'excMice']; % Overwrite the root folder in case we are working on a new machine (Dropbox sync issues)
  dset = dset.loadTrialsSummary([LfileRoot datasetList.dataset{it} '\trials\L_tile.mat']);
  % Save the retinotopy in the correct place
  dset.data.areas = dset.data.trialsSummary.aux.visarea_ed;
  dset.save();
end
 
%% Load and store retinotopy data and visual field sign

for it = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it});
  dset.info.retinotopy = datasetList.retinotopy{it};
  dset.data.retinotopy = load([retinotopyRootFolder dset.info.retinotopy{1} '\ANALYZED\' dset.info.retinotopy{2} '\AnalyzedRet.mat']);
  [VFS, VFS_thr] = getVisualSign(dset.data.retinotopy);
  dset.data.visualFieldSign = VFS;
  dset.data.visualFieldSignThreshold = VFS_thr;
  dset.save();
end

%% Plot all the retinotopies

hFig = createCenteredFigure('height', 10);
p = panel();
p.margin = 5;
p.pack(2, Ndatasets);
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  p(1, it1).select();
  retData =  dset.data.visualFieldSign;
  imagesc(retData);
  colormap(gca,parula);
  setImageAxis();
  title(dset.name);
  p(2, it1).select();
  retData = dset.data.areas;
  colormap(gca, parula);
  imagesc(retData);
  setImageAxis();
end

%% Store area centroid information

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  retData = dset.data.areas;
  areaIdx = [datasetList.V1{it1} datasetList.PM{it1} datasetList.RL{it1} datasetList.AM{it1}];
  areaNames = {'V1' 'PM' 'RL' 'AM'};
  
  dset.data.retinotopyCoordinates = struct;
  dset.data.retinotopyCoordinates.areaID = struct;
  for it2 = 1:length(areaIdx)
    [r, c] = find(retData == areaIdx(it2));
    dset.data.retinotopyCoordinates.(areaNames{it2}) = [mean(c) mean(r)];
    dset.data.retinotopyCoordinates.areaID.(areaNames{it2}) = areaIdx(it2);
  end
  dset.data.retinotopyCoordinates.alignmentTriangle = polyshape([dset.data.retinotopyCoordinates.V1; dset.data.retinotopyCoordinates.PM; dset.data.retinotopyCoordinates.RL]);
  dset.save();
end

%% Plot the alignment triangles

hFig = createCenteredFigure('height', 10);
p = panel();
p.margin = 5;
p.pack(2, Ndatasets);
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  p(1, it1).select();
  retData =  dset.data.visualFieldSign;
  imagesc(retData);
  colormap(gca,parula);
  setImageAxis();
  title(dset.name);
  
  p(2, it1).select();
  retData = dset.data.areas;
  imagesc(retData);
  colormap(gca, parula);
  hold on;
  h = plot(dset.data.retinotopyCoordinates.alignmentTriangle, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 1);
  setImageAxis();
end

%% Generate a new dataset from the Allen data (for retinotopic alignment)

allen = wf.dataset('allen', 'rootFolder', projectFolder, 'comments', 'allen reference dataset');
allen.rootFolder = [pwd filesep 'data' filesep 'excMice']; % Overwrite for machine swaps
  
allen.data.areas = imread('refData/allenRG.tif');
allen.data.labels = imread('refData/allenR.tif');
allenImgRet = imread('refData/allenAvgRetinotopyNew.png');
allenImgRet = imresize(allenImgRet, size(allen.data.areas, [1 2]));
allen.data.visualFieldSign = allenImgRet(:, end:-1:1, :);

% Recolor the areas
retData = double(allen.data.areas);
uc = unique(retData);
uc(uc == 255) = [];
ucorder = randperm(length(uc));
for it2 = 1:length(uc)
  retData(retData==uc(it2)) = -ucorder(it2);
end
retData = abs(retData);
retData(retData == 255) = 0;
allen.data.areas = uint8(retData);

%V1 PM RL
coordsV1 = [3483 3617];
coordsPM = [3217 3364];
coordsRL = [3820 3350];
colV1 = allen.data.areas(coordsV1(2), coordsV1(1));
colPM = allen.data.areas(coordsPM(2), coordsPM(1));
colRL = allen.data.areas(coordsRL(2), coordsRL(1));

areaIdx = [colV1 colPM colRL];
areaNames = {'V1' 'PM' 'RL'};
allen.data.retinotopyCoordinates = struct;
for it2 = 1:length(areaIdx)
  [r, c] = find(allen.data.areas == areaIdx(it2));
  allen.data.retinotopyCoordinates.(areaNames{it2}) = [mean(c) mean(r)];
end
allen.data.retinotopyCoordinates.alignmentTriangle = polyshape([allen.data.retinotopyCoordinates.V1; allen.data.retinotopyCoordinates.PM; allen.data.retinotopyCoordinates.RL]);
allen.save();

%% Do the alignment transformations
% This is the reference animal (dataset ID) all other retinotopies will be aligned to
refAnimal = 5;

for it1 = 1:(Ndatasets+1)
  if(it1 <= Ndatasets)
    dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  else
    dset = wf.dataset.load(projectFolder, 'allen');
  end
  dset = dset.alignToReference(datasetList.dataset{refAnimal});
  dset.save();
end


%% Plot the transformations
% Reference image
R = imref2d([67, 67]);

hFig = createCenteredFigure('height', 12);
p = panel();
p.margin = 5;
p.pack(3, Ndatasets+1);
for it1 = 1:(Ndatasets+1)
  if(it1 <= Ndatasets)
    dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  else
    dset = wf.dataset.load(projectFolder, 'allen');
  end
  p(1, it1).select();
  retData =  dset.data.visualFieldSign;
  imagesc(retData);
  colormap(gca,parula);
  setImageAxis();
  title(dset.name);
  
  p(2, it1).select();
  
  retData = dset.data.areas;
  
  imagesc(retData);
  colormap(gca, parula);
  setImageAxis();
  hold on;
  plot(dset.data.retinotopyCoordinates.alignmentTriangle, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 1);
  p(3, it1).select();
  
  retData = dset.data.areas;
  
  B = imwarp(retData, dset.data.alignment.transform, 'nearest', 'OutputView', R);
  B = round(B);
  colormap(gca, parula);
  imagesc(B);
  setImageAxis();
  hold on;
  plot(dset.data.alignment.triangle, 'EdgeColor', [0 0 0], 'FaceColor', [0 0 0], 'LineWidth', 1);
  plot(dset.data.alignment.referenceTriangle, 'EdgeColor', [1 0 0], 'FaceColor', [1 0 0], 'LineWidth', 0.5);
  
  xlim(R.XWorldLimits);
  ylim(R.YWorldLimits);
end

%% Plot V1, RL, PM on top of the allen reference
hFig = createCenteredFigure();

allen = wf.dataset.load(projectFolder, 'allen');

img = allen.data.labels(:, :, 2);
img(img < 255) = 0;
imshow(double(~~img));
hold on;

%cols = lines(4);
h =[];
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  
  curTransform = dset.data.alignment.transform;
  % Transform to the Allen reference
  newT = curTransform.T*inv(allen.data.alignment.transform.T);
  newT(:, 3) = [0 0 1]; % needed due to rounding errors in the prev transform
  curTransform.T = newT;
  
  retData = dset.data.areas;
  fnames = fieldnames(dset.data.retinotopyCoordinates.areaID);
  ax = gca;
  ax.ColorOrderIndex = 1;
  for it2 = 1:length(fnames)
    img = zeros(size(retData));
    img(retData == dset.data.retinotopyCoordinates.areaID.(fnames{it2})) = 1;
    B = bwboundaries(img);
    [xn,yn] = transformPointsForward(curTransform, B{1}(:,2), B{1}(:,1));
    %hb = plot(xn, yn, 'DisplayName', fnames{it2});
    hb = plot(smooth(xn), smooth(yn), 'DisplayName', fnames{it2});
    if(it1 == 1)
      h = [h; hb];
    end
  end
end
setImageAxis();
axis(gca, 'equal', 'tight', 'ij');
xlim([2600 4600]);
ylim([2600 4400]);

l = legend(h);
l.Location = 'SouthEast';

%% Save the allen areas into each dataset
R = imref2d([67, 67]);
dsetAllen = wf.dataset.load(projectFolder, 'allen');
retData = dsetAllen.data.areas;  
areasAllen = round(imwarp(retData, dsetAllen.data.alignment.transform, 'nearest', 'OutputView', R));
% Save also our redefinition of Allen areas
areasAllenLocaNMF = zeros(size(areasAllen));
areasNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RSP'};
areasID = {14, 29, 25, 26, 4, 27, 21, 11, [7 12 19 32 34], [28 31]};
for it1 = 1:length(areasNames)
  for it2 = 1:length(areasID{it1})
    areasAllenLocaNMF(areasAllen == areasID{it1}(it2)) = it1;
  end
end

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  dset.data.areasAllen = areasAllen;
  dset.data.areasAllenLocaNMF = areasAllenLocaNMF;
  dset.data.areasAllenLocaNMFnames = areasNames;
  dset.save();
  %[hFig, dset] = dset.plotAverageImage('session', sessList(4), 'superimposeAreas', true, 'registered', true, 'registeredRef', imref2d([67 67]), 'export', 'pdf');
end

%%
%%% -----------------------------------------------------------------------
%%% After all the alignments we can prepare the datasets for locaNMF
%%% -----------------------------------------------------------------------

%% Generate main time series
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  % Input values are the time points of: prestimulus start, stimulus presentation, closed loop start, last time point (here 5 sec of closed loop)
  dset = dset.generateTimeSeries(0, 1, 2.5, 7.5);
end


%% Concatenate event times to harmonize the dataset
fnames = {'trialStartedTime';'intermissionStartedTime';'intermissionEndedTime';'stimulusBackgroundStartedTime';'stimulusCueStartedTime';'OpenLoopStartedTime';'OpenLoopEndedTime';'interactiveStartedTime';'inputThresholdCrossedTime';'interactiveEndedTime';'responseMadeTime';'feedbackStartedTime';'feedbackType';'feedbackPositiveStartedTime';'responseMadeID';'inputThresholdCrossedID';'feedbackPositiveEndedTime';'feedbackEndedTime';'stimulusCueEndedTime';'stimulusBackgroundEndedTime';'trialEndedTime';'feedbackNegativeStartedTime';'feedbackNegativeEndedTime'};
fnamesCond = {'visCueContrast';'responseForThreshold';'feedbackForResponse';'orientation';'distBetweenTargets';'targetWidth';'targetThreshold';'cyclesPerDegree';'gaussianSigmaP';'phase';'targetOrientation';'horizontalTarget';'useMosaic'};

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  Ntrials = size(dset.data.trialsSummary.trial,1);
  eventTimes = struct;
  eventTimesCondition = struct;
  for it2 = 1:length(fnames)
    eventTimes.(fnames{it2}) = nan(Ntrials, 1);
  end
  eventTimes.lastImagingTime = nan(Ntrials, 1);
  eventTimes.imagingFrames = nan(Ntrials, 1);
  eventTimes.lastWheelTime = nan(Ntrials, 1);
  eventTimes.wheelFrames = nan(Ntrials, 1);
  eventTimes.pupilFrames = nan(Ntrials, 1);
  for it2 = 1:length(fnamesCond)
    switch fnamesCond{it2}
      case {'visCueContrast', 'responseForThreshold', 'feedbackForResponse', 'orientation', 'phase'}
        eventTimesCondition.(fnamesCond{it2}) = nan(Ntrials, 2);
      otherwise
        eventTimesCondition.(fnamesCond{it2}) = nan(Ntrials, 1);
    end
  end
  ncbar('Storing event times');
  for it3 = 1:Ntrials
    trialData = dset.pullSingleTrial(it3);
    for it2 = 1:length(fnames)
      if(isfield(trialData.eventtimes, fnames{it2}) && ~isempty(trialData.eventtimes.(fnames{it2})))
        eventTimes.(fnames{it2})(it3) = trialData.eventtimes.(fnames{it2});
      end
    end
    eventTimes.lastImagingTime(it3) = trialData.time(end);
    eventTimes.imagingFrames(it3) = length(trialData.time);
    eventTimes.lastWheelTime(it3) = trialData.wheelpos(end, 2);
    eventTimes.firstWheelTime(it3) = trialData.wheelpos(1, 2);
    eventTimes.wheelFrames(it3) = size(trialData.wheelpos, 1);
    eventTimes.pupilFrames(it3) = length(trialData.parea);
    for it2 = 1:length(fnamesCond)
      if(isfield(trialData.eventtimes.condition, fnamesCond{it2}) && ~isempty(trialData.eventtimes.condition.(fnamesCond{it2})))
        eventTimesCondition.(fnamesCond{it2})(it3, :) = trialData.eventtimes.condition.(fnamesCond{it2});
      end
    end
    ncbar.update(it3/Ntrials);
  end
  ncbar.close();
  dset.data.eventTimes = eventTimes;
  dset.data.eventTimesCondition = eventTimesCondition;
  dset.save();
end

%% Generate the pixel mask across sessions and the final aligned version of imaging and behavioral data
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  
  sessList = dset.getValidSessions();
  pixelVar1 = zeros(numel(dset.data.areas), length(sessList));
  pixelVar2 = zeros(numel(dset.data.areas), length(sessList));
  for it2 = 1:length(sessList)
    data = dset.pullSessionTrials(sessList(it2), 'fps', 30, 'bar', false, 'force', true, 'maxTime', 5);
    for it3 = 1:size(pixelVar1, 1)
      tmp = data.DFF(:,it3,:);
      pixelVar1(it3, it2) = nanstd(tmp(:));
      pixelVar2(it3, it2) = nanmean(tmp(:));
    end
  end

  %%% This should be the joint sessions mask
  pixelVarS = min(pixelVar1')./mean(pixelVar1');

  figure;
  imagesc(reshape(pixelVarS, size(dset.data.areas)));
  % Removing the pixels with v low variance (background)
  pixelVarS(pixelVarS < 0.1) = 0;
  dset.data.sessAvgStdImg = reshape(mean(pixelVar1,2), size(dset.data.areas));
  figure;
  imagesc(reshape(~~pixelVarS, size(dset.data.areas)));
  dset.data.areasPixelMask = reshape(~~pixelVarS, size(dset.data.areas));
  dset.save();
end

%% Plot 1 session average image for consistency checks

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  sessList = dset.getValidSessions();
  [hFig, dset] = dset.plotAverageImage('session', sessList(1), 'superimposeAreas', 'allenLocaNMF', 'registered', true, 'registeredRef', imref2d([67 67]));
end

%% Check the valid pixel set to use for everything else

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  sessList = dset.getValidSessions();
  [hFig, dset] = dset.plotAverageImage('session', sessList(1), 'superimposeAreas', 'allenLocaNMF', 'registered', true, 'registeredRef', imref2d([67 67]));
  regMask = imwarp(dset.data.areasPixelMask, dset.data.alignment.transform, 'nearest', 'OutputView', R);
  hFig.Children(2).Children(end).CData = hFig.Children(2).Children(end).CData.*regMask;
end

%% Do some data cleaning
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  dset.rootFolder = [pwd filesep 'data' filesep 'excMice']; % Overwrite for machine swaps
  sessList = dset.getValidSessions();
  for it2 = 1:length(sessList)
    dset.cleanupSessionTrials(sessList(it2), 'removeFnans', true, 'realignPrestim', 1, 'lastTimePoint', 4.999, 'removeSimSaccades', true, 'simSaccadesTimeThreshold', 3, 'openLoopDuration', 1.5);
  end
end

%% Concatenate all data

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  data = dset.joinCleanSessions();
end

%% Compress split and register the joint data

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  dset.splitAndCompressJointSessions();
  dset.registerCompressedSessions(imref2d([67 67]));
end


%% Optional: Bit to z-score the pupil area information (since it might not be)

for itAnimal = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{itAnimal});

  dataBehavior = load(fullfile(dset.rootFolder, dset.name, dset.dataFolder, sprintf('%s_jointSessionsBehavior.mat', dset.name)));
  figure;
  hold on;
  for it2 = dset.getValidSessions()
    validTrials = find(dataBehavior.sessionID == it2);
    data = dataBehavior.pupilArea(:, validTrials);
    [f, xi] = ksdensity(data(:));
    plot(xi, f);
    dataM = nanmean(data(:));
    dataS = nanstd(data(:));
    data = (data-dataM)./dataS;
    dataBehavior.pupilArea(:, validTrials) = data;
  end
  save(fullfile(dset.rootFolder, dset.name, dset.dataFolder, sprintf('%s_jointSessionsBehaviorZScored.mat', dset.name)), '-struct', 'dataBehavior');
end

%% Do some data cleaning based on consistency checks: Remvoe trials with inconsistent times or that used a slightly different protocol

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  sessList = dset.getValidSessions();
  for it2 = 1:length(sessList)
    dset.cleanupSessionTrials(sessList(it2), 'removeFnans', true, 'realignPrestim', 1, 'lastTimePoint', 4.999, 'removeSimSaccades', true, 'simSaccadesTimeThreshold', 3, 'openLoopDuration', 1.5);
  end
end


%% Perform SVD on all the datasets

for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  data = dset.SVDregisteredSessions('targetVar', 99.9);
end

%% Save SVD metadata in a different file (for faster loading)
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  data = load(fullfile(dset.rootFolder, dset.name, dset.dataFolder, sprintf('%s_SVD.mat', dset.name)));
  data = rmfield(data, 'coeff');
  data = rmfield(data, 'score');
  save(fullfile(dset.rootFolder, dset.name, dset.dataFolder, sprintf('%s_SVD_minimal.mat', dset.name)), '-struct' ,'data');
  save(fullfile(dset.rootFolder, 'data', sprintf('%s_SVD_minimal.mat', dset.name)), '-struct' ,'data');
end

%% Now you run locaNMF (python version)


