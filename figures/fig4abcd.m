clearvars;
load('fig4abcd.mat');

%% Final figure

labels = {'attention', 'stimulus', 'wheel mov', 'saccades', 'choice early', 'choice late'};
perm = [1 3 4 6 5 2];

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 15;
p.pack(1, 4);

p(1, 1).select();
imagesc(a_angleMat);
cmap = parula(256);
colormap(cmap);
cb = colorbar;
cb.Label.String = 'angle (º)';
caxis([0 90]);
cb.Ticks = 0:30:90;
set(gca, 'XTick', 1:6);
set(gca, 'YTick', 1:6);
set(gca, 'XTickLabels', labels(perm));
set(gca, 'XTickLabelRotation', 90);
set(gca, 'YTickLabels', labels(perm));
title('Angle between state vectors');
axis ij equal square;
xlim([0.5 6.5]);
ylim([0.5 6.5]);

p(1, 2).select();
[~, ~, perm] = dendrogram(b_Y, 'Labels', labels);
set(gca, 'XTickLabelRotation', 90);
ylabel('Hierarchical clustering');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

p(1,3).select();
cmap = cblindmap();
cmap = cmap(1:end,:);
plot(t, ym, '-', 'Color', cmap(1, :));
hold on;
fill([t' fliplr(t')], [yu' fliplr(yl')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel('time from mov (s)');
ylabel('\theta_{C-M}');
xlim([-0.5 0.5]);
upper = [89 89];
lower = [79 79]; 
hold on;
set(gca,'XTick', -0.5:0.5:0.5);
set(gca,'YTick', 0:30:90);
ylim([60 90]);
plot([0 0], ylim, 'k--');
spaceOutAxes(gca, 0.05);
offsetAxes(gca);
x = xlim;
plot_handle = fill([x fliplr(x)],[upper fliplr(lower)],[0.2 0.2 0.2],'FaceAlpha',0.2 , 'EdgeColor', 'none');

%
p(1, 4).select();
cmap = cblindmap;
cmap = cmap([4 1 2 3 5 6 7], :);
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(sdi, 'Compression', 10, 'PointSize', 4, 'Color', cmap, 'MarkerFaceColor', cmap, 'MarkerEdgeColor', cmap);
set(gca, 'XTickLabels', labels);
set(gca, 'XTickLabelRotation', 90);
ylabel('Spatial-Distribution index (%)');
offsetAxes();
spaceOutAxes();
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
%exportgraphics(hFig, 'svCorrAreas.pdf');



