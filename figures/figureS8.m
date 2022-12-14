% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs8.mat');


hFig = createCenteredFigure('width', 6, 'height', 5);
nu = 87.12;
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_a, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'C vs I', 'C replicas', 'I replicas'});
xlim([0.5 3.5]);
ylim([0 90]);
set(gca, 'YTick', 0:15:90);
ylabel('attention angle');
plot(xlim, [1 1]*nu,'k--');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();


hFig = createCenteredFigure('width', 6, 'height', 5);

nur =  0.0022;
[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_b, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'C vs I', 'C replicas', 'I replicas'});
xlim([0.5 3.5]);
ylabel('attention correlation');
ylim([0 1]);
plot(xlim, [1 1]*nur,'k--');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();

hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data_c, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'C with I', 'I with C'});
xlim([0.5 3.5]);
ylabel('d''');
ylim([0 1.8])
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();

