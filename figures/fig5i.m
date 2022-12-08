clearvars;
load('fig5i.mat');
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();
imagesc(attList, attList, diffMat);
set(gca, 'XTick', 0:0.25:1);
set(gca, 'YTick', 0:0.25:1);
xlabel('Attention level');
ylabel('Attention level');
axis equal square tight ij;
cb = colorbar;
caxis([0 30]);
cb.Label.String = 'Angle between choice vectors(°)';
xlim([0 1]);
ylim([0 1]);
%spaceOutAxes();
%offsetAxes(gca);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);