% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs14.mat');

hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_a, 'Compression', 10, 'PointSize', 16);
hold on;

set(gca,'Xtick', 1:3);
set(gca,'XtickLabel', {'dorsal', 'posterior', 'ventral'});
xlim(0.5+[0 3]);

ylabel('d'' increase in choice axis');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();

hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_b, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Pre vs Post', 'Pre replicas', 'Post replicas'});
xlim([0.5 3.5]);
ylim([0 90]);
set(gca, 'YTick', 0:15:90);
ylabel('choice angle');
plot(xlim, [1 1]*85.5,'k--');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();

