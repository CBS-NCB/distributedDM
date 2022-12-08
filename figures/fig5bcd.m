clearvars;
load('fig5bcd.mat');
projLabels = {'left', 'right'};

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.1 0.35 0.35 0.2]);

%%% Projection trajs
p(1 ,2).select();
cmap = cblindmap();
hold on;
h1 = [];
h2 = [];
fill([t' fliplr(t')], [y1u' fliplr(y1l')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', projLabels{1}))];
fill([t' fliplr(t')], [y2u' fliplr(y2l')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(t, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', projLabels{2}))];
yl = ylim;
%ylim([0 yl(2)*1.2]);
plot([1 1], yl, 'k--');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([0.5 2.5]);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([t' fliplr(t')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(t, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([1 1], yl, 'k--');
xlabel('time (s)');
ylabel('d''');
xlim([0.5 2.5]);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(CC, 'XData', [0 2.5], 'YData', [0 2.5]);
x.AlphaData = ~isnan(CC);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'Correlation';
caxis([0 1]);
axis square;
hold on;
xl = xlim;
yl = ylim;
set(gca, 'XTick', 0:2);
set(gca, 'YTick', 0:2);
plot([1 1], yl, 'k--');
plot(xl, [1 1], 'k--');
xlabel('Time from mov (s)');
ylabel('Time from mov (s)');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);