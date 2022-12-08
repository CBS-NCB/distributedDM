%%
clearvars;
load('figs12.mat');
%% The figure plot

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.06 0.2 0.2 0.11]);

p(1, 2).select();
hList = plotPairedData(data_a);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'low', 'low with high axis'});
xlim([0.5 2.5]);
ylabel('d''');
h4 = gca;
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
h4.Position(1) = h4.Position(1)+0.045;
ylim([0 2]);

% Now the d prime spatial map
p(1, 3).select();
hList = plotPairedData(data_b);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'high', 'high with low axis'});
xlim([0.5 2.5]);
ylabel('d''');
h4 = gca;
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
h4.Position(1) = h4.Position(1)+0.045;
ylim([0 2]);
