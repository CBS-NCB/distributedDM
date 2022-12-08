clearvars;
load('fig5j.mat');
cmap = parula(length(attList)+1);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();
hold on;
for it1 = 1:length(attList)
  plot(diffList, attListPsy(:, it1), '-', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(it1, :), 'Color', cmap(it1, :))
end

set(gca, 'XTick', -90:30:90);
set(gca, 'YTick', 0:0.2:1);
spaceOutAxes();
offsetAxes(gca);
xlabel('Angle difference |\theta_L|-|\theta_R|');
ylabel('Prob of choosing right');
cb = colorbar;
cb.Label.String = 'Attention';
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);