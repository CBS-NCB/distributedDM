clearvars;
load('fig3defg.mat');

hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);

p(2,1).pack(1, [0.25 0.25 0.25 0.25]);
p(2,1,1,1).select();

cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig(regionColor, :);
hold on;
for it1 = 1:length(regionNames)
  h = plot(areaDt, areaD(:, it1), 'Color', cmap(it1, :));
end

xlabel('time from mov (s)');
ylabel('d''');
ylim([0 1.2]);
plot([0 0], ylim, 'k--');
spaceOutAxes(gca, 0.05);
offsetAxes(gca);

p(2,1,1,2).select();
hold on;
for it1 = 1:length(regionNames)
  h1 = plot(areaDfitX, areaDfit(:, it1), 'Color', cmap(it1, :), 'LineWidth', 0.5);
end

xlabel('time from mov (s)');
ylabel('d''');
xlim([-1 0.2]);
ylim([0 1.2]);
plot([0 0], ylim, 'k--');
spaceOutAxes(gca, 0.05);
offsetAxes(gca);

p(2,1,1,3).select();
hold on;
cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig([1 regionColor], :);
errorbar(slopeDx, slopeDm, slopeDl, slopeDu, 'o', 'Color', cmap(3,:), 'MarkerSize', 4);

xlim([0 max(regionPosition)+1]);
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
plot(xlim,[0 0], 'k--');
offsetAxes();
spaceOutAxes();
ylabel('Slope pre movement')

p(2,1,1,4).select();
hold on;
errorbar(slopeDx, slopeTm, slopeTl, slopeTu, 'o', 'Color', cmap(4,:), 'MarkerSize', 4);

xlim([0 max(regionPosition)+1]);
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
plot(xlim,[0 0], 'k--');
offsetAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);