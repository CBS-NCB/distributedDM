clearvars;
load('figs9.mat');
%%
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
XPositionsOrig = [ 1     3     5     6     8     9    10    12    14    15    17];
regionColor = [2     3     3     4     4     4     5     6     6     7];
regionOrder = [1     2     3     4     6     8     9     5     7    10];
regionPosition = [3 5 6 8 9 10 12 14 15 17]; % In the new order
regionPosition
hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 -1], 1);

p(1,1).pack(1, [0.15 0.35 0.15 0.35]);


p(1, 1, 1, 1).select();

movData = dpre;
cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(gca, cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
caxis([0 max(movData(:))]);

tf = [];
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    tf =[tf; text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), regionNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0])];
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

p(1, 1, 1, 3).select();
movData = dpost;

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(gca, cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
caxis([0 max(movData(:))]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    tf =[tf; text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), regionNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0])];
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

p(1, 1, 1, 2).select();
cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig([1 regionColor], :);
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(univarpre, 'Width', 1, 'Whiskers', 'CI', 'Compression', 10, 'PointSize', 4, 'XPositions', XPositionsOrig, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
  
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
ylabel('d''');
offsetAxes();
spaceOutAxes();
h6 = gca;
 
p(1, 1, 1, 4).select();
cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig([1 regionColor], :);
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(univarpost, 'Width', 1, 'Whiskers', 'CI', 'Compression', 10, 'PointSize', 4, 'XPositions', XPositionsOrig, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
  
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
ylabel('d''');
offsetAxes();
spaceOutAxes();
h8 = gca;

h6.Position(1) = h6.Position(1)+0.025;
h8.Position(1) = h8.Position(1)+0.025;

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

