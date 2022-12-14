% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%% PANEL A
clearvars;
load('data/figs13a.mat');

hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Compression', 10, 'PointSize', 16);
hold on;
set(gca, 'XTick', [1 2 3 4 5]);
set(gca, 'XTickLabel', {'choice', '⟂ att', '⟂ mov', '⟂ sacc', '⟂ stim'});
xlim([0.5 5.5]);
ylabel('d''');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL B
clearvars;
load('data/figs13b.mat');
hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'mov', '⟂ sacc', '⟂ att'});
xlim([0.5 5.5]);
ylabel('d''');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();


%% PANEL C
clearvars;
load('data/figs13c.mat');
hFig = createCenteredFigure('width', 6, 'height', 5);

[xPositions, yPositions, Label, RangeCut, hList] = UnivarScatterV2(data, 'Compression', 10, 'PointSize', 16);
hold on;
% plot(data','-')
set(gca, 'XTick', [1 2 3 4]);
set(gca, 'XTickLabel', {'low', 'low on high', 'high', 'high on low'});
xlim([0.5 5.5]);
ylabel('d''');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
offsetAxes();
spaceOutAxes();