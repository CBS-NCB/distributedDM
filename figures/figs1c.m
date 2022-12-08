clearvars;
load('figs1c.mat');
%%
hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(2);

fill([t(:)' fliplr(t(:)')], [y1u(:)' fliplr(y1l(:)')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', 'locaNMF');
fill([t(:)' fliplr(t(:)')], [y2u(:)' fliplr(y2l(:)')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t2, y2m, '-', 'Color', cmap(2, :), 'DisplayName', 'SVD');

title('total variance explained');
xlabel('Num components');
l = legend([h1(end) h2(end)]);
l.Box = 'off';
l.Location = 'best';
xlim([0 100]);
set(gca,'XScale','log');
spaceOutAxes(gca);
xl = xlim;
xlim([0.9 xl(2)]);
set(gca, 'XTick', [1 10 100]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);
ylabel('Total variance explained (%)');

