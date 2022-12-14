% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs10.mat');


hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);
p(1,1).pack(1, [0.06 0.32 0.22 0.32 0.08]);

vectorModes = {'left','right'}; % Vector will be V1 - V2


%%% The projection
p(1, 1, 1, 2).select();
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
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% The corr mat
p(1, 1, 1, 3).select();
imagesc(corrMat, 'XData', t, 'YData', t);
xlim([0 2.5]);
ylim([0 2.5]);
cmap = parula;
hold on;
plot([1 1], ylim, 'k--');
plot(xlim, [1 1], 'k--');
xlabel('time (s)');
ylabel('time (s)');
axis square ij tight;
caxis([0 1]);
colormap(gca,cmap);
hcb = colorbar;
hcb.Label.String = 'correlation';
h2 = gca;

%%% The d prime
cmap = cblindmap();
cmap = cmap(1:end,:);
p(1, 1, 1, 4).select();

hold on;
h1 = [];

fill([t(:)' fliplr(t(:)')], [du(:)' fliplr(dl(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, dm, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}))];

ylim([0 1.2])
set(gca, 'YTick', 0:0.5:1);
yl = ylim;
plot([1 1], yl, 'k--');

xlabel('time (s)');
ylabel('d''');
xlim([0 2.5]);
set(gca, 'XTick', 0:2);
offsetAxes();
spaceOutAxes();
h3 = gca;

% Now the d prime spatial map
p(1, 1, 1, 5).select();

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(univard, 'Compression', 10, 'PointSize', 16, 'Color', cmap(3,:), 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:));
set(gca, 'XTick', 1);
set(gca, 'XTickLabel', 'd'' (t = 2.5 s)');
xlim([0.5 1.5]);
ylabel('d''');
h4 = gca;

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
