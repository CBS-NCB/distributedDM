% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs6.mat');

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.1 0.3 0.6]);
p(1, 2).select();
cmap = lines(8);
hold on;
h1 = [];
h2 = [];
vectorModes = {'corr. pre (t=-0.2s)',  'corr. post (t=-0.2s)'};
t = torig;
fill([t(:)' fliplr(t(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', vectorModes{1}))];
fill([t(:)' fliplr(t(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(t, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', vectorModes{2}))];

ylim([0 1]);
xlabel('time (s)');
ylabel('correlation');
legend([h1 h2]);
xlim([-0.8 0.8])
offsetAxes(); 
spaceOutAxes();
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

p(1, 3).select();

H = dendrogram(tree, 0, 'ColorThreshold','default');
tt = get(gca,'XTickLabels');
idx = str2num(tt);
set(gca,'XTickLabels',sprintf('%.2f\n',t(idx)))
%set(gca,'YTickLabels',sprintf('%2d\n',round(tim(idx)*1000)));
set(findall(gcf,'-property','FontName'),'FontName', 'Arial');
set(findall(gcf,'-property','FontSize'),'FontSize', 8);
ylabel('linkage dist.');
xlabel('time from mov.');
set(gca,'XTickLabelRotation', 90)

