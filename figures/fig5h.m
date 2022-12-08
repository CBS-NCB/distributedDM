clearvars;
load('fig5h.mat');
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();

imagesc(diffList, diffList, diffMat);
set(gca, 'XTick', 0:30:90);
set(gca, 'YTick', 0:30:90);
xlabel('Difficulty ||\theta_L-\theta_R||');
ylabel('Difficulty ||\theta_L-\theta_R||');
axis equal square tight ij;
cb = colorbar;
caxis([0 30]);
cb.Label.String = 'Angle between choice vectors(°)';
xlim([0 90]);
ylim([0 90]);
%spaceOutAxes();
%offsetAxes(gca);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);