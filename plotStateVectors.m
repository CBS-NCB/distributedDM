%% Load project data

clearvars;
addpath(genpath(pwd));
projectFolder = [pwd filesep 'data' filesep 'excMice'];
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


%%
%%% -----------------------------------------------------------------------
%%% Saccades
%%% -----------------------------------------------------------------------
load('compSacc.mat');


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

exportgraphics(hFig, 'stimSlopes.pdf');

%% Plot the single projection split across animals

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
load('compSacc.mat');

wSize = 3;

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
%caxis([0 1.2]);

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

exportgraphics(hFig, 'fig2Sacc.pdf', 'ContentType', 'vector');


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
%caxis([0 1.2]);
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

exportgraphics(hFig, 'fig2SaccED.pdf', 'ContentType', 'vector');

