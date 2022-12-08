clearvars;
load('figs1defg.mat');

%%
cmap = cblindmap;
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 15;
p.pack([0.125 -1], 1);
p(1,1).pack(1, [0.25 0.25 0.25 0.25]);


p(1,1,1,1).select();
errorbar(1:10, d_y, d_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Rel. var explained');

ax = gca;
set(gca, 'YTick', [0:0.2:1]);
yyaxis right;
ax2 = gca;
plot(1:10, d2_y, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(2, :));
ylabel('Normalized area');
xlim([0.5 10.5]);
axes(ax);
set(gca, 'YTick', [0:0.2:1]);
spaceOutAxes(gca);
box off;
set(gca, 'XTickLabelRotation', 90);
ax.YAxis(2).Limits = ax.YAxis(1).Limits;

p(1,1,1,2).select();
errorbar(1:10, e_y, e_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Rel. var explained / pixel');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

%%% Num components per region
p(1,1,1,3).select();
errorbar(1:10, f_y, f_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Avg num components');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

%%% Num components per region normalized
p(1,1,1,4).select();
errorbar(1:10, g_y, g_e, 'o-', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :));
set(gca, 'XTick', 1:10);
set(gca, 'XTickLabels', labels);
ylabel('Norm avg num components');
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(gca, 'XTickLabelRotation', 90);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

