clearvars;
load('fig5a.mat');

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);

p(1, 1).select();

cmap = cblindmap();
hold on;
errorbar(angleList, ylowm, ylowl, ylowu, '-o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :), 'Color', cmap(1, :), 'MarkerSize', 4)
errorbar(angleList, yhighm, yhighl, yhighu', '-o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(2, :), 'Color', cmap(2, :), 'MarkerSize', 4)
l = legend('low attention', 'high attention');
l.Location = 'southeast';
l.Box = 'off';
l.Position(1) = l.Position(1) + 0.05;
set(gca, 'XTick', -90:30:90);
set(gca, 'YTick', 0:0.2:1);
spaceOutAxes();
offsetAxes(gca);
xlabel('Angle difference |\theta_L|-|\theta_R|');
ylabel('Prob of choosing right');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);