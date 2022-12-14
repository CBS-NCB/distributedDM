% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%% PANEL A
clearvars;
load('data/figs1a.mat');

hFig = createCenteredFigure('width', 5, 'height', 5);
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'k');
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNamesOrig{it3}, 'HorizontalAlignment', 'center');
  end
end
setImageAxis();
xlim([1 67]);
ylim([1 67]);
box on;

for it1 = 1:7
  ax = gca;
  ax.ColorOrderIndex = 1;
  for it2 = 1:3
    hb = plot(area_x{it1}{it2}, area_y{it1}{it2}, 'DisplayName', fnames{it2});
  end
end

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);


%% PANEL B
clearvars;
load('data/figs1b.mat');

hFig = createCenteredFigure('width', 8, 'height', 6);
p = panel();
p.margin = 2;
p.pack(2, 5);
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'VL', 'RS'};
curComp = 0;

for it1 = 1:2
  for it2 = 1:5
    curComp = curComp + 1;
    p(it1,it2).select();
    hold off;
    imagesc(1-baseImgFull(:, :, curComp));
    %imagesc(baseImgFull(:, :, curComp));
    hold on;
    axis equal;
    for it3 = 1:length(contourListTransform)
      plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'k');
    end
    setImageAxis();
    xlim([1 67]);
    ylim([1 67]);
    colormap([1 1 1;parula]);
    title(regionNames{curComp});  
  end
end
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

%% PANEL C
clearvars;
load('data/figs1c.mat');

hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(2);

fill([t(:)' fliplr(t(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', 'locaNMF');
fill([t(:)' fliplr(t(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t2, y2m, '-', 'Color', cmap(2, :), 'DisplayName', 'SVD');

title('total variance explained');
xlabel('Num components');
l = legend([h1(end) h2(end)]);
l.Box = 'off';
l.Location = 'best';
xlim([0 100]);
set(gca,'XScale','log');
spaceOutAxes(gca);
xl = xlim;
xlim([0.9 xl(2)]);
set(gca, 'XTick', [1 10 100]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);
ylabel('Total variance explained (%)');


%% PANELS D E F G
clearvars;
load('data/figs1defg.mat');

cmap = cblindmap;
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 15;
p.pack([0.125 -1], 1);
p(1,1).pack(1, [0.25 0.25 0.25 0.25]);

p(1,1,1,1).select();
errorbar(1:10, d_y, d_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Rel. var explained');

ax = gca;
set(gca, 'YTick', [0:0.2:1]);
yyaxis right;
ax2 = gca;
plot(1:10, d2_y, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(2, :));
ylabel('Normalized area');
xlim([0.5 10.5]);
axes(ax);
set(gca, 'YTick', [0:0.2:1]);
spaceOutAxes(gca);
box off;
set(gca, 'XTickLabelRotation', 90);
ax.YAxis(2).Limits = ax.YAxis(1).Limits;

p(1,1,1,2).select();
errorbar(1:10, e_y, e_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Rel. var explained / pixel');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

%%% Num components per region
p(1,1,1,3).select();
errorbar(1:10, f_y, f_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Avg num components');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

%%% Num components per region normalized
p(1,1,1,4).select();
errorbar(1:10, g_y, g_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Norm avg num components');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
