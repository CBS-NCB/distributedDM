clearvars;
load('figs11ab.mat');

%%
cmap = lines(6);
hFig = createCenteredFigure('width', 21, 'height', 6);
p = panel();
p.margin = 10;
p.pack(1, [0.3, 0.3 0.4]);
p(1,1).select();
hold on;

y = a1;idx = 1;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));

y = a2;idx = 2;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));

y = a3;idx = 3;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h3 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));

y = a4;idx = 4;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h4 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));



xlabel('time from movement (s)');
ylabel('Wheel position (a.u.)');

xlim([-1 1]);
yl = ylim;
plot([0 0], yl, 'k--');
plot(xlim, [0 0], 'k-');
l = legend([h1 h2 h3 h4], 'easy left', 'hard left', 'easy right', 'hard right');
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();


p(1,2).select();
hold on;
t = bt;
y = b1;idx = 1;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));

y = b2;idx = 2;
fill([t(:)' fliplr(t(:)')], [y(:,3)' fliplr(y(:,2)')], cmap(idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = plot(t, y(:,1), '-', 'Color', cmap(idx, :));


xlabel('time from movement (s)');
ylabel('d''');
xlim([-1 1]);
plot([0 0], yl, 'k--');
plot(xlim, [0 0], 'k-');
l = legend([h1 h2], 'left easy-hard', 'right easy-hard');
l.Box = 'off';
l.Location = 'best';
ylim([-1 1])
yl = ylim;

offsetAxes();
spaceOutAxes();
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
