%%% -----------------------------------------------------------------------
%%% This script should be run after locaNMF has run on all animals
%%% See main.m
%%% -----------------------------------------------------------------------

%% Load project data

clearvars;
addpath(genpath(pwd));
projectFolder = [pwd filesep 'data' filesep 'excMice'];
load([projectFolder filesep 'projectData.mat']);

%% Load the LocaNMF decomp for the extended data plots

% Variance explained threshold of SVD and locaNMF bits (for normalization purposes)
SVDvarianceCap = 0.999;
locaNMFvarianceCap = 0.985;

locaVarExplained = cell(Ndatasets, 1);
SVDvarExplained = cell(Ndatasets, 1);
sortedVar = cell(Ndatasets, 1);
sortedVarIdx = cell(Ndatasets, 1);
sortedVarRaw = cell(Ndatasets, 1);
sortedVarIdxRaw = cell(Ndatasets, 1);

regionSizesAreas = cell(Ndatasets, 1);
comPerRegion = cell(Ndatasets, 1);
comPerRegionNorm = cell(Ndatasets, 1);
for it1 = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{it1});
  dataLoca = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_locaNMF_hdf5_output.mat', dset.name)));
  dataSVD = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_SVD_minimal.mat', dset.name)));
  baseComp = load(fullfile(dset.rootFolder, dset.dataFolder, sprintf('%s_baseComponents.mat', dset.name)));

  regionMasks = baseComp.a_dict.region_mats;
  regionCompSizes = cell2mat(dataLoca.region_ranks(2:end));
  regionIdx = [];
  for it = 1:length(regionCompSizes)
    regionIdx = [regionIdx, ones(1, regionCompSizes(it))*it];
  end
  regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'VL', 'RS'};
  regionSizes = sum(regionMasks==0,2);
  %%
  v1 = var(dataLoca.compData);
  vexp  = zeros(length(v1), 1);
  for it0 = 1:length(v1)
    vexp(it0) = sum(dataLoca.compXY(it0, :).^2*v1(it0));
  end
  maxVar = SVDvarianceCap*locaNMFvarianceCap*sum(sort(vexp,'descend'));

  vexpArea = zeros(size(regionSizes));
  for it = 1:length(vexpArea)
    vexpArea(it) = sum(vexp(regionIdx == it));
  end
  % Variables needed here
  locaVarExplained{it1} = SVDvarianceCap*locaNMFvarianceCap*cumsum(sort(vexp,'descend'))/max(cumsum(sort(vexp,'descend')));
  SVDvarExplained{it1} = cumsum(dataSVD.latent);
  [sortedVar{it1}, sortedVarIdx{it1}] = sort(vexpArea./regionSizes/maxVar, 'descend');
  [sortedVarRaw{it1}, sortedVarIdxRaw{it1}] = sort(vexpArea/maxVar, 'descend');

  regionSizesAreas{it1} = regionSizes(sortedVarIdxRaw{it1})/max(regionSizes);
  comPerRegion{it1} = double(regionCompSizes);
  comPerRegionNorm{it1} = double(regionCompSizes')./regionSizes;
  it1
end

%% Variance explained plots

figure;
hold on;
h1 = plot(locaVarExplained,'o-', 'DisplayName', 'loca');
h2 = plot(SVDvarExplained, 'o-', 'DisplayName', 'SVD');
title('total variance explained');
xlabel('Num components');
xlim([0 25]);
legend([h1 h2]);

%% Sorted var per region with size normalization

figure;
plot(1:10, sortedVar, 'o-')
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames(sortedVarIdx));
xlabel('Region idx');
ylabel('Norm var per pixel');
title('Localized var explained');


%% Sorted var per region w/o normalization

figure;
plot(1:10, sortedVarRaw, 'o-')
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames(sortedVarIdxRaw));
xlabel('Region idx');
ylabel('Norm var');
yyaxis right;
plot(1:10, regionSizesAreas,'o-');
ylabel('Normalized area');
title('Localized var explained raw');

%% Num components per region

figure;
bar(comPerRegion);
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames);
ylabel('Num components');

%% Num components per region normalized

figure;
bar(comPerRegionNorm);
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames);
ylabel('Num components per pixel');

%% Norm variance and area plots

figure;
[sortedVar, sortedVarIdx] = sort(vexp./regionSizes(regionIdx)/maxVar, 'descend');
[sortedVarRaw, sortedVarIdxRaw] = sort(vexp/maxVar, 'descend');

regionIdx(sortedVarIdx(1:10))
plot(1:10, sortedVar(1:10), 'o-')
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames(regionIdx(sortedVarIdx(1:10))));
xlabel('Region idx'); 
ylabel('Norm var per pixel');
yyaxis right;
plot(1:10, regionSizes(regionIdx(sortedVarIdx(1:10)))/max(regionSizes),'o-');
ylabel('Normalized area');


%% More plots

vexpArea = zeros(size(regionSizes));
for it = 1:length(vexpArea)
  vexpArea(it) = sum(vexp(regionIdx == it));
end
figure;
[sortedVar, sortedVarIdx] = sort(vexpArea./regionSizes/maxVar, 'descend');
[sortedVarRaw, sortedVarIdxRaw] = sort(vexpArea/maxVar, 'descend');

regionIdx(sortedVarIdx)
plot(1:10, sortedVar, 'o-')
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames(sortedVarIdx));
xlabel('Region idx');
ylabel('Norm var per pixel');
yyaxis right;
plot(1:10, regionSizes(sortedVarIdx)/max(regionSizes),'o-');
ylabel('Normalized area');
title('Localized var explained');


%% More plots

vexpArea = zeros(size(regionSizes));
for it = 1:length(vexpArea)
  vexpArea(it) = sum(vexp(regionIdx == it));
end
figure;
[sortedVar, sortedVarIdx] = sort(vexpArea./regionSizes/maxVar, 'descend');
[sortedVarRaw, sortedVarIdxRaw] = sort(vexpArea/maxVar, 'descend');

regionIdx(sortedVarIdx)
plot(1:10, sortedVarRaw, 'o-')
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames(sortedVarIdxRaw));
xlabel('Region idx');
ylabel('Norm var');
yyaxis right;
plot(1:10, regionSizes(sortedVarIdxRaw)/max(regionSizes),'o-');
ylabel('Normalized area');
title('Localized var explained raw');

%% More plots

figure;
bar(double(regionCompSizes));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames);
ylabel('Num components');

%% More plots

figure;
bar(double(regionCompSizes')./regionSizes);
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', regionNames);
ylabel('Num components per pixel');


