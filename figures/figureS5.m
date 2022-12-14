% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs5.mat');

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.2 0.3 0.3 0.2]);

p(1, 1).select();

x = imagesc(data_a);
x.AlphaData = ~isnan(data_a);
axis square;
setImageAxis();
cmap = parula(256);
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'coeff';
set(gca, 'Position', axp-[0.025 0 0 0]);
axis off


%%% The d prime
cmap = cblindmap();
p(1, 2).select();
hold on;
vectorModes = {'before corr', 'after corr'};
xRange = [-0.5 1.5];
fill([t1(:)' fliplr(t1(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t1, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}));
fill([t2(:)' fliplr(t2(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t2, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', vectorModes{2}));

ylim([0 2.1])
set(gca, 'YTick', 0:0.5:2);
yl = ylim;

plot([evTime evTime], yl, 'k--');

l = legend(h1);
l.Box = 'off';
l.Location = 'best';
xlabel('time (s)');
ylabel('d''');
title('stimulus d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

p(1, 3).select();
hold on;

fill([t1b(:)' fliplr(t1b(:)')], [y1ub(:)' fliplr(y1lb(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t1b, y1mb, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}));
fill([t2b(:)' fliplr(t2b(:)')], [y2ub(:)' fliplr(y2lb(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t2b, y2mb, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', vectorModes{2}));

ylim([0 1.7])
set(gca, 'YTick', 0:0.5:2);
yl = ylim;

plot([evTimeb evTimeb], yl, 'k--');

l = legend(h1);
l.Box = 'off';
l.Location = 'best';
% title(datasetList.dataset{itAnimal});
xlabel('time (s)');
ylabel('d''');
title('wheel mov d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

p(1, 4).select();


h1 = errorbar(mxs, mys, sys, sys, sxs, sys,...
                                    '.', 'Color', cmap(1,:));
hold on;
h2 = errorbar(mxw, myw, syw, syw, sxw, syw,...
                                    '.', 'Color', cmap(2,:));
xl = xlim;
yl = ylim;
M = max([xl(2) yl(2)]);
m = min([xl(1) yl(1)]);
plot([m M], [m M], 'k--');
l = legend([h1 h2], 'stim', 'wheel');
xlabel('d'' before corr.');
ylabel('d'' after corr.');
l.Box = 'off';
l.Location = 'southeast';

