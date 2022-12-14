% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs3.mat');

regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
XPositionsOrig = [ 1     3     5     6     8     9    10    12    14    15    17];
regionColor = [2     3     3     4     4     4     5     6     6     7];
regionOrder = [1     2     3     4     6     8     9     5     7    10];
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 29.7);
p = panel();
p.margin = 10;
p.pack([0.15 0.15 0.15 0.15 -1], [0.1 0.2 0.45 0.25]);


for itType = 1:4
  switch itType
    case 1
      titleName = 'stimulus';
      matRange = [-0.5 0.5];
      linePos = 0;
      matLabel = 'time from stimulus (s)';
    case 2
      titleName = 'wheel movements';
      matRange = [-1 1];
      linePos = 0;
      matLabel = 'time from movement (s)';
    case 3
      titleName = 'saccades';
      matRange = [-1 1];
      linePos = 0;
      matLabel = 'time from saccade (s)';
    case 4
      titleName = 'attention';
      matRange = [0 2.5];
      linePos = 1;
      matLabel = 'trial time (s)';
  end

  % Now the d prime spatial map
  p(itType, 2).select();


  movData = ddata{itType};
  x = imagesc(movData);
  x.AlphaData = ~isnan(movData);
  setImageAxis();
  axp = get(gca, 'Position');
  c = colorbar;
  c.Label.String = 'd''';
  %c.Location = 'west';
  %caxis(gca, [0 max(movData(:))]);
  caxis(gca, [0 1.2]);
  c.Ticks = 0:0.2:1.2;
  %caxis([0 1.2]);
  hold on;
  for it3 = 1:length(contourListTransform)
    plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
  end
  for it3 = 1:length(contourListTransform)
    if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
      text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), regionNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
    end
  end
  colormap(gca, parula(128));
  axis off;
  xlim([1 67]);
  ylim([1 67]);  
  set(gca, 'Position', axp-[0.035 0 0 0]);
  title(titleName);

  % Now the d prime spatial map
  p(itType, 3).select();

  
  cmapOrig = [0 0 0; cblindmap];
  cmap = cmapOrig([1 regionColor], :);
  
  data = univardata{itType};
  
  [xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Width', 1, 'Whiskers', 'CI', 'Compression', 10, 'PointSize', 4, 'XPositions', XPositionsOrig, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
  set(gca, 'XTick', XPositionsOrig);
  set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
  set(gca, 'XTickLabelRotation', 90);
  
  ylim([0 2]);
  xlim([0 18]);
  offsetAxes();
  spaceOutAxes();
  yl = ylim;
  if(itType == 1)
    pp = patch([0 2 2 0]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(1, :), 'EdgeColor', 'none');
    pp = patch([2 4 4 2]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(2, :), 'EdgeColor', 'none');
    pp = patch([4 7 7 4]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(3, :), 'EdgeColor', 'none');
    pp = patch([7 11 11 7]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(4, :), 'EdgeColor', 'none');
    pp = patch([11 13 13 11]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(5, :), 'EdgeColor', 'none');
    pp = patch([13 16 16 13]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(6, :), 'EdgeColor', 'none');
    pp = patch([16 18 18 16]+[0.1 -0.1 -0.1 0.1], [0.95 0.95 1 1]*yl(2), cmapOrig(7, :), 'EdgeColor', 'none');
  end
  

  p(itType, 4).select();
  cmap = parula(256);
  imagesc(corrdata{itType}, 'XData', linspace(matRange(1), matRange(2), length(corrdata{itType})), 'YData', linspace(matRange(1), matRange(2), length(corrdata{itType})));
  xlim(matRange);
  ylim(matRange);
  hold on;
  plot([1 1]*linePos, ylim, 'k--');
  plot(xlim, [1 1]*linePos, 'k--');
  xlabel(matLabel);
  ylabel(matLabel);
  axis square ij;
  caxis([0 1]);
  colormap(gca, cmap);
  cb = colorbar;
  cb.Label.String = 'correlation';
  %p(itType,2).marginright = 17;
end

%export_fig('2ext_AreaMasks.pdf','-nocrop');


p(1,1).select();
%figure;
areaGroupIDs = {[1], [2, 3], [4, 6, 8], [9], [5, 7], [10]};
areaGroupNames = {'V1', 'ventral', 'parietal', 'dorsal', 'somatosensory', 'retrosplenial'};

movData = ddata{5};
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
axp = get(gca, 'Position');
caxis(gca, [1 max(movData(:))]);
colormap(gca, cmapOrig(2:7,:));

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
%for it3 = 1:length(groupLabelsPosition)
%  text(groupLabelsPosition(it3, 1), groupLabelsPosition(it3, 2), areaGroupNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
%end
axis off;
xlim([1 67]);
ylim([1 67]);  
axp = gca;
axp.Position(1:2) = [0.05 0.8];
axis equal;

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);


