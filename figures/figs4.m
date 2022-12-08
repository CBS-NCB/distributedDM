%%
clearvars;
load('figs4.mat');
%%
%hFig = createCenteredFigure('width', 6, 'height', 2);
hFig = createCenteredFigure('width', 21, 'height', 29.7);
p = panel();
p.margin = 10;
p.pack([0.14 0.1 0.1 0.1 -1], [0.1 0.35 0.2 0.35]);

titleName = 'contralateral stimulus';
matRange = [-0.5 0.5];
linePos = 0;
matLabel = 'time from stimulus (s)';
eventTimeDataT = (1/30:1/30:2.5)';
vectorModes = {'horizontal','vertical'};
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
XPositionsOrig = [ 1     3     5     6     8     9    10    12    14    15    17];
regionColor = [2     3     3     4     4     4     5     6     6     7];
regionOrder = [1     2     3     4     6     8     9     5     7    10];

%%% The projection
p(1, 2).select();
cmap = cblindmap();
cmap = cmap(1:end,:);
hold on;
h1 = [];
h2 = [];

fill([t(:)' fliplr(t(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}))];

fill([t(:)' fliplr(t(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(t, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', vectorModes{2}))];

yl = ylim;
plot([1 1], ylim, 'k--');
set(gca, 'XTick', 0:2);

xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([0 2.5]);
%xlim([-1 1]);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

% Now the d prime spatial map
p(1, 3).select();
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
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
p(1, 4).select();
cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig([1 regionColor], :);
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(univardata, 'Width', 1, 'Whiskers', 'CI', 'Compression', 10, 'PointSize', 4, 'XPositions', XPositionsOrig, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
set(gca, 'XTick', XPositionsOrig);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
%ylim([0 1.2*max(data(:))]);
%ylim([0 2]);
xlim([0 18]);
offsetAxes();
spaceOutAxes();
yl = ylim;

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);


%exportgraphics(hFig, 'sup4_ori.pdf', 'ContentType', 'vector');

%%



