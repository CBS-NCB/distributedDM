% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs15.mat');

regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
vectorModes = {'left', 'right'};      

hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);

p(1,1).pack(1, [0.3 0.15 0.3 0.15 0.1]);

cmap = cblindmap();
cmap = cmap(1:end,:);
p(1,1,1,1).select();

hold on;
h1 = [];
h2 = [];

fill([t1(:)' fliplr(t1(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(t1, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}));


yl = ylim;
plot([0 0], yl, 'k--');

xlabel('time (s)');
ylabel('d''');
xlim([-1 1]);
set(gca, 'XTick', 0:2);
offsetAxes();
spaceOutAxes();
h3 = gca;

%%% The corr mat
p(1,1,1,2).select();


cmap = parula(256);
x = imagesc(movData_a);
x.AlphaData = ~isnan(movData_a);
setImageAxis();
colormap(gca, cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
caxis([0 max(movData_a(:))]);

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
axis equal;
%%%

%%% The d prime
cmap = cblindmap();
cmap = cmap(1:end,:);
p(1,1,1,3).select();

hold on;

fill([t2(:)' fliplr(t2(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(t2, y2m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}));
yl = ylim;
plot([0 0], yl, 'k--');

xlabel('time (s)');
ylabel('d''');
xlim([-1 1]);
set(gca, 'XTick', 0:2);
offsetAxes();
spaceOutAxes();
h3 = gca;

%%% The corr mat
p(1,1,1,4).select();

cmap = parula(256);
x = imagesc(movData_b);
x.AlphaData = ~isnan(movData_b);
setImageAxis();
colormap(gca, cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
caxis([0 max(movData_b(:))]);

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
axis equal;

p(1, 1, 1, 5).select();

XPositionsOrig = 1:3;

cmapOrig = cblindmap;
cmap = cmapOrig([3 3 3], :);
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_c, 'Width', 1, 'Whiskers', 'CI', 'Compression', 10, 'PointSize', 4, 'XPositions', XPositionsOrig, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
  
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'orig', 'no stim mov', 'no sacc'});
set(gca, 'XTickLabelRotation', 90);
ylabel('d''');
offsetAxes();
spaceOutAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

